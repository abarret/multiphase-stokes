#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/app_namespaces.h>

#include <multiphase/MultiphaseStaggeredStokesBlockOperator.h>
#include <multiphase/MultiphaseStaggeredStokesBlockPreconditioner.h>
#include <multiphase/utility_functions.h>

#include <SAMRAIVectorReal.h>

namespace multiphase
{

static constexpr int UN_COMP_IDX = 0;
static constexpr int US_COMP_IDX = 1;
static constexpr int P_COMP_IDX = 2;

MultiphaseStaggeredStokesBlockPreconditioner::MultiphaseStaggeredStokesBlockPreconditioner(std::string object_name)
    : d_object_name(std::move(object_name)),
      d_thn_ths_sq_var(new SideVariable<NDIM, double>(d_object_name + "::ThnThsSq"))
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_thn_ths_sq_idx =
        var_db->registerVariableAndContext(d_thn_ths_sq_var, var_db->getContext(d_object_name + "::CTX"));
}

MultiphaseStaggeredStokesBlockPreconditioner::~MultiphaseStaggeredStokesBlockPreconditioner()
{
    // If we are allocated, deallocate before we delete the object
    if (d_is_initialized) deallocateSolverState();
}

bool
MultiphaseStaggeredStokesBlockPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& x,
                                                          SAMRAIVectorReal<NDIM, double>& b)
{
    // Pull out components of the vectors
    const int un_idx = x.getComponentDescriptorIndex(UN_COMP_IDX);
    const int us_idx = x.getComponentDescriptorIndex(US_COMP_IDX);
    const int p_idx = x.getComponentDescriptorIndex(P_COMP_IDX);

    const int bn_idx = b.getComponentDescriptorIndex(UN_COMP_IDX);
    const int bs_idx = b.getComponentDescriptorIndex(US_COMP_IDX);
    const int bp_idx = b.getComponentDescriptorIndex(P_COMP_IDX);

    // Create the individual vectors, add components to the vectors
    SAMRAIVectorReal<NDIM, double> u_vec, p_vec;
    u_vec.addComponent(
        x.getComponentVariable(0), x.getComponentDescriptorIndex(0), x.getControlVolumeIndex(0), d_hier_sc_data_ops);
    u_vec.addComponent(
        x.getComponentVariable(1), x.getComponentDescriptorIndex(1), x.getControlVolumeIndex(1), d_hier_sc_data_ops);
    p_vec.addComponent(
        x.getComponentVariable(2), x.getComponentDescriptorIndex(2), x.getControlVolumeIndex(2), d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM, double> bu_vec, bp_vec;
    bu_vec.addComponent(
        b.getComponentVariable(0), b.getComponentDescriptorIndex(0), b.getControlVolumeIndex(0), d_hier_sc_data_ops);
    bu_vec.addComponent(
        b.getComponentVariable(1), b.getComponentDescriptorIndex(1), b.getControlVolumeIndex(1), d_hier_sc_data_ops);
    bp_vec.addComponent(
        b.getComponentVariable(2), b.getComponentDescriptorIndex(2), b.getControlVolumeIndex(2), d_hier_cc_data_ops);

    /*
     * The preconditioner solves:
     * (A G )(u) = (f_u)
     * (0 -S)(p)   (f_p)
     *
     * in which S is the least squares commutator approximation to the Schur complement: S = (G^T G)*(G^T*A*G)^(-1)*(G^T
     * G)
     *
     * We solve this sytem by using the block decomposition:
     * (A G )^-1 = (A^-1 0)(I -G)(I  0   )
     * (0 -S)      (0    I)(0  I)(0 -S^-1)
     */

    // First solve the approximate Schur complement system
    d_pressure_solver->solveSystem(p_vec, bp_vec);

    // Now compute the appropriate RHS for the velocity sub problem.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const Box<NDIM>& box = patch->getBox();
            const double* const dx = pgeom->getDx();

            Pointer<SideData<NDIM, double>> bun_data = bu_vec.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double>> bus_data = bu_vec.getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM, double>> p_data = p_vec.getComponentPatchData(0, *patch);
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(d_thn_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(box, axis); si; si++)
                {
                    const SideIndex<NDIM>& idx = si();
                    // Note: This does not account for any synchronization that needs to occur at coarse-fine
                    // interfaces. Consider using HierarchyMathOps::grad() instead.
                    const double dp = ((*p_data)(idx.toCell(1)) - (*p_data)(idx.toCell(0))) / dx[axis];
                    const double thn = 0.5 * ((*thn_data)(idx.toCell(1)) + (*thn_data)(idx.toCell(0)));
                    (*bun_data)(idx) = (*bun_data)(idx) + thn * dp;
                    (*bus_data)(idx) = (*bun_data)(idx) + convertToThs(thn) * dp;
                }
            }
        }
    }

    // Now solve the velocity block
    d_stokes_solver->solveSystem(u_vec, bu_vec);
}

void
MultiphaseStaggeredStokesBlockPreconditioner::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                                    const SAMRAIVectorReal<NDIM, double>& b)
{
    LinearSolver::initializeSolverState(x, b);
    d_hierarchy = x.getPatchHierarchy();
#ifndef NDEBUG
    TBOX_ASSERT(d_hierarchy.getPointer() == b.getPatchHierarchy());
#endif
    d_coarsest_ln = 0;
    d_finest_ln = d_hierarchy->getFinestLevelNumber();

    d_hier_cc_data_ops = new HierarchyCellDataOpsReal<NDIM, double>(d_hierarchy, d_coarsest_ln, d_finest_ln);
    d_hier_sc_data_ops = new HierarchySideDataOpsReal<NDIM, double>(d_hierarchy, d_coarsest_ln, d_finest_ln);

    SAMRAIVectorReal<NDIM, double> u_vec, p_vec;
    u_vec.addComponent(
        x.getComponentVariable(0), x.getComponentDescriptorIndex(0), x.getControlVolumeIndex(0), d_hier_sc_data_ops);
    u_vec.addComponent(
        x.getComponentVariable(1), x.getComponentDescriptorIndex(1), x.getControlVolumeIndex(1), d_hier_sc_data_ops);
    p_vec.addComponent(
        x.getComponentVariable(2), x.getComponentDescriptorIndex(2), x.getControlVolumeIndex(2), d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM, double> bu_vec, bp_vec;
    bu_vec.addComponent(
        b.getComponentVariable(0), b.getComponentDescriptorIndex(0), b.getControlVolumeIndex(0), d_hier_sc_data_ops);
    bu_vec.addComponent(
        b.getComponentVariable(1), b.getComponentDescriptorIndex(1), b.getControlVolumeIndex(1), d_hier_sc_data_ops);
    bp_vec.addComponent(
        b.getComponentVariable(2), b.getComponentDescriptorIndex(2), b.getControlVolumeIndex(2), d_hier_cc_data_ops);

    // Allocate patch data and compute thn_ths_sq.
    allocate_patch_data(d_thn_ths_sq_idx, d_hierarchy, d_solution_time, d_coarsest_ln, d_finest_ln);

    // Note we need ghost cells of thn to compute thn_ths_sq.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp{ ITC(
        d_thn_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE", "LINEAR", false, d_thn_bc_coefs) };
    HierarchyGhostCellInterpolation hier_ghost_fill;
    hier_ghost_fill.initializeOperatorState(ghost_cell_comp, d_hierarchy, d_coarsest_ln, d_finest_ln);
    hier_ghost_fill.fillData(d_solution_time);

    // Now compute thn_ths_sq.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(d_thn_idx);
            Pointer<SideData<NDIM, double>> thn_ths_sq_data = patch->getPatchData(d_thn_ths_sq_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SideIndex<NDIM>& idx = si();

                    const double thn = 0.5 * ((*thn_data)(idx.toCell(0)) + (*thn_data)(idx.toCell(1)));
                    const double ths = convertToThs(thn);
                    (*thn_ths_sq_data)(idx) = thn * thn + ths * ths;
                }
            }
        }
    }

    // Create the pressure solver
    auto manager = CCPoissonSolverManager::getManager();
    d_pressure_solver = manager->allocateSolver(d_pressure_solver_type,
                                                d_object_name + "::pressure_solve",
                                                d_pressure_solver_db,
                                                "stokes_pressure_",
                                                d_pressure_precond_type,
                                                d_object_name + "::pressure_precond",
                                                d_pressure_precond_db,
                                                "stokes_pressure_precond_");
    d_pressure_coefs.setCConstant(0.0);
    d_pressure_coefs.setDPatchDataId(d_thn_ths_sq_idx);
    d_pressure_solver->setPoissonSpecifications(d_pressure_coefs);
    d_pressure_solver->initializeSolverState(p_vec, bp_vec);

    // Create the velocity solver
    Pointer<MultiphaseStaggeredStokesBlockOperator> stokes_op =
        new MultiphaseStaggeredStokesBlockOperator("stokes_op", true, d_params);
    d_stokes_solver =
        std::make_unique<PETScKrylovLinearSolver>("stokes_solver", d_stokes_solver_db, d_object_name + "_stokes_");
    d_stokes_solver->setOperator(stokes_op);
    d_stokes_precond_op =
        new MultiphaseStaggeredStokesBlockFACOperator(d_object_name + "::stokes_precond_op", "", d_params);
    d_stokes_precond = new FullFACPreconditioner(
        "stokes_precond", d_stokes_precond_op, d_stokes_precond_db, d_object_name + "_stokes_pc_");
    d_stokes_precond->transferToDense(d_thn_idx);
    d_stokes_precond_op->setThnIdx(d_thn_idx);
    d_stokes_solver->setPreconditioner(d_stokes_precond);
    d_stokes_solver->initializeSolverState(u_vec, bu_vec);

    // Setup any needed communication algorithms and scratch indices
}

void
MultiphaseStaggeredStokesBlockPreconditioner::deallocateSolverState()
{
    LinearSolver::deallocateSolverState();

    // Deallocate any scratch data.
}

} // namespace multiphase
