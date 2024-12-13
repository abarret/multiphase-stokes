#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/app_namespaces.h>

#include <multiphase/MultiphaseStaggeredStokesBlockOperator.h>
#include <multiphase/MultiphaseStaggeredStokesBlockPreconditioner.h>
#include <multiphase/fd_operators.h>
#include <multiphase/utility_functions.h>

#include <SAMRAIVectorReal.h>

namespace multiphase
{

static constexpr int UN_COMP_IDX = 0;
static constexpr int US_COMP_IDX = 1;
static constexpr int P_COMP_IDX = 2;

MultiphaseStaggeredStokesBlockPreconditioner::MultiphaseStaggeredStokesBlockPreconditioner(
    std::string object_name,
    const MultiphaseParameters& params,
    Pointer<Database> input_db)
    : d_params(params),
      d_thn_ths_sq_var(new SideVariable<NDIM, double>(d_object_name + "::ThnThsSq")),
      d_un_scr_var(new SideVariable<NDIM, double>(d_object_name + "::UN_SCR")),
      d_us_scr_var(new SideVariable<NDIM, double>(d_object_name + "::US_SCR")),
      d_fn_scr_var(new SideVariable<NDIM, double>(d_object_name + "::FN_SCR")),
      d_fs_scr_var(new SideVariable<NDIM, double>(d_object_name + "::FS_SCR")),
      d_p_scr_var(new CellVariable<NDIM, double>(d_object_name + "::P_SCR"))
{
    LinearSolver::init(std::move(object_name), d_homogeneous_bc);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_thn_ths_sq_idx =
        var_db->registerVariableAndContext(d_thn_ths_sq_var, var_db->getContext(d_object_name + "::CTX"));
    d_un_scr_idx = var_db->registerVariableAndContext(d_un_scr_var, var_db->getContext(d_object_name + "::CTX"), 1);
    d_us_scr_idx = var_db->registerVariableAndContext(d_us_scr_var, var_db->getContext(d_object_name + "::CTX"), 1);
    d_fn_scr_idx = var_db->registerVariableAndContext(d_fn_scr_var, var_db->getContext(d_object_name + "::CTX"), 1);
    d_fs_scr_idx = var_db->registerVariableAndContext(d_fs_scr_var, var_db->getContext(d_object_name + "::CTX"), 1);

    // Set up the solver databases
    d_stokes_solver_db = input_db->getDatabase("stokes_solver_db");

    d_stokes_precond_db = input_db->getDatabase("stokes_precond_db");

    d_pressure_solver_type = input_db->getStringWithDefault("pressure_solver_type", d_pressure_solver_type);
    d_pressure_solver_db = input_db->getDatabase("pressure_solver_db");
    d_pressure_precond_type = input_db->getStringWithDefault("pressure_precond_type", d_pressure_precond_type);
    d_pressure_precond_db = input_db->getDatabase("pressure_precond_db");
}

MultiphaseStaggeredStokesBlockPreconditioner::~MultiphaseStaggeredStokesBlockPreconditioner()
{
    // If we are allocated, deallocate before we delete the object
    if (d_is_initialized) deallocateSolverState();
}

void
MultiphaseStaggeredStokesBlockPreconditioner::setThnIdx(const int thn_idx)
{
    d_thn_idx = thn_idx;
}

void
MultiphaseStaggeredStokesBlockPreconditioner::setCAndDCoefficients(const double C, const double D)
{
    d_C = C;
    d_D = D;
}

bool
MultiphaseStaggeredStokesBlockPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& x,
                                                          SAMRAIVectorReal<NDIM, double>& b)
{
    // Pull out components of the vectors
    x.setToScalar(0.0); // Currently inputs are NaN's, so set them to 0 for now...
    const int un_idx = x.getComponentDescriptorIndex(UN_COMP_IDX);
    const int us_idx = x.getComponentDescriptorIndex(US_COMP_IDX);
    const int p_idx = x.getComponentDescriptorIndex(P_COMP_IDX);

    const int bn_idx = b.getComponentDescriptorIndex(UN_COMP_IDX);
    const int bs_idx = b.getComponentDescriptorIndex(US_COMP_IDX);
    const int bp_idx = b.getComponentDescriptorIndex(P_COMP_IDX);

    // Create the individual vectors, add components to the vectors
    SAMRAIVectorReal<NDIM, double> u_vec(d_object_name + "::U_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln),
        p_vec(d_object_name + "::P_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln);
    u_vec.addComponent(
        x.getComponentVariable(0), x.getComponentDescriptorIndex(0), x.getControlVolumeIndex(0), d_hier_sc_data_ops);
    u_vec.addComponent(
        x.getComponentVariable(1), x.getComponentDescriptorIndex(1), x.getControlVolumeIndex(1), d_hier_sc_data_ops);
    p_vec.addComponent(
        x.getComponentVariable(2), x.getComponentDescriptorIndex(2), x.getControlVolumeIndex(2), d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM, double> bu_vec(d_object_name + "::BU_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln),
        bp_vec(d_object_name + "::BP_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln);
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
     * in which S is the least squares commutator approximation to the Schur complement: S =
     * (G^T*G)*(G^T*A*G)^(-1)*(G^T*G)
     *
     * Note that using this preconditioner with an exact S gives the preconditioned operator
     * (I        0)
     * (D*A^(-1) I)
     *
     * We solve this system by using the block decomposition:
     * (A G )^-1 = (A^-1 0)(I -G)(I  0   )
     * (0 -S)      (0    I)(0  I)(0 -S^-1)
     */

    // First solve the approximate Schur complement system
    // The approximate Schur complement is S*p=db with S = (G^T*G)*(G^T*A*G)^(-1)*(G^T*G)
    // We first solve G^T*G*p_1 = b_p.
    d_pressure_solver->solveSystem(p_vec, bp_vec);

    // Fill ghost cells for p_vec.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_comps{ ITC(
        p_vec.getComponentDescriptorIndex(0), "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN") };
    HierarchyGhostCellInterpolation hier_ghost_fill;
    hier_ghost_fill.initializeOperatorState(ghost_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);
    hier_ghost_fill.fillData(d_solution_time);

    // Now compute the RHS for the next solve:
    // b_p = G^T*A*G*p_1
    multiphase_grad_on_hierarchy(*d_hierarchy,
                                 d_un_scr_idx,
                                 d_us_scr_idx,
                                 d_thn_idx,
                                 p_vec.getComponentDescriptorIndex(0),
                                 1.0,
                                 true,
                                 d_coarsest_ln,
                                 d_finest_ln);
    accumulateMomentumWithoutPressureConstantCoefficient(*d_hierarchy,
                                                         d_fn_scr_idx,
                                                         d_fs_scr_idx,
                                                         d_un_scr_idx,
                                                         d_us_scr_idx,
                                                         d_thn_idx,
                                                         d_params,
                                                         1.0,
                                                         1.0,
                                                         d_coarsest_ln,
                                                         d_finest_ln);
    applyCoincompressibility(*d_hierarchy,
                             bp_vec.getComponentDescriptorIndex(0),
                             d_fn_scr_idx,
                             d_fs_scr_idx,
                             d_thn_idx,
                             1.0,
                             d_coarsest_ln,
                             d_finest_ln);

    // Finally, solve for p:
    // G^T*G*p = b_p.
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

    return true;
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
    set_valid_level_numbers(*d_hierarchy, d_coarsest_ln, d_finest_ln);

    d_hier_cc_data_ops = new HierarchyCellDataOpsReal<NDIM, double>(d_hierarchy, d_coarsest_ln, d_finest_ln);
    d_hier_sc_data_ops = new HierarchySideDataOpsReal<NDIM, double>(d_hierarchy, d_coarsest_ln, d_finest_ln);

    SAMRAIVectorReal<NDIM, double> u_vec(d_object_name + "::U_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln),
        p_vec(d_object_name + "::P_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln);
    u_vec.addComponent(
        x.getComponentVariable(0), x.getComponentDescriptorIndex(0), x.getControlVolumeIndex(0), d_hier_sc_data_ops);
    u_vec.addComponent(
        x.getComponentVariable(1), x.getComponentDescriptorIndex(1), x.getControlVolumeIndex(1), d_hier_sc_data_ops);
    p_vec.addComponent(
        x.getComponentVariable(2), x.getComponentDescriptorIndex(2), x.getControlVolumeIndex(2), d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM, double> bu_vec(d_object_name + "::BU_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln),
        bp_vec(d_object_name + "::BP_VEC", d_hierarchy, d_coarsest_ln, d_finest_ln);
    bu_vec.addComponent(
        b.getComponentVariable(0), b.getComponentDescriptorIndex(0), b.getControlVolumeIndex(0), d_hier_sc_data_ops);
    bu_vec.addComponent(
        b.getComponentVariable(1), b.getComponentDescriptorIndex(1), b.getControlVolumeIndex(1), d_hier_sc_data_ops);
    bp_vec.addComponent(
        b.getComponentVariable(2), b.getComponentDescriptorIndex(2), b.getControlVolumeIndex(2), d_hier_cc_data_ops);

    // Allocate patch data and compute thn_ths_sq.
    allocate_patch_data({ d_thn_ths_sq_idx, d_un_scr_idx, d_us_scr_idx, d_fn_scr_idx, d_fs_scr_idx },
                        d_hierarchy,
                        d_solution_time,
                        d_coarsest_ln,
                        d_finest_ln);

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
    d_stokes_precond_op->setThnIdx(d_thn_idx);
    stokes_op->setThnIdx(d_thn_idx);
    stokes_op->setCandDCoefficients(d_C, d_D);
    d_stokes_precond_op->setCandDCoefficients(d_C, d_D);
    d_stokes_solver->setPreconditioner(d_stokes_precond);
    d_stokes_solver->initializeSolverState(u_vec, bu_vec);

    d_stokes_precond->transferToDense(d_thn_idx);

    // Setup any needed communication algorithms and scratch indices
}

void
MultiphaseStaggeredStokesBlockPreconditioner::deallocateSolverState()
{
    LinearSolver::deallocateSolverState();

    // Deallocate any scratch data.
    deallocate_patch_data({ d_thn_ths_sq_idx, d_un_scr_idx, d_us_scr_idx, d_fn_scr_idx, d_fs_scr_idx },
                          d_hierarchy,
                          d_coarsest_ln,
                          d_finest_ln);

    d_coarsest_ln = IBTK::invalid_level_number;
    d_finest_ln = IBTK::invalid_level_number;
}

} // namespace multiphase
