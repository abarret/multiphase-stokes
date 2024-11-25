#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/app_namespaces.h>

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
    SAMRAIVectorReal<NDIM, double> u_vec(d_object_name + "::u", d_hierarchy, d_coarsest_ln, d_finest_ln);

    SAMRAIVectorReal<NDIM, double> p_vec(d_object_name + "::p", d_hierarchy, d_coarsest_ln, d_finest_ln);

    SAMRAIVectorReal<NDIM, double> bu_vec(d_object_name + "::bu", d_hierarchy, d_coarsest_ln, d_finest_ln);

    SAMRAIVectorReal<NDIM, double> bp_vec(d_object_name + "::bp", d_hierarchy, d_coarsest_ln, d_finest_ln);

    // Now solve the pressure system
    d_pressure_solver->solveSystem(p_vec, bp_vec);

    // Now fill in the correct components for the velocity solve
    // bu_vec = bu_vec - grad(p_vec)

    // Now solve the velocity system
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

    // Setup any needed communication algorithms and scratch indices
}

void
MultiphaseStaggeredStokesBlockPreconditioner::deallocateSolverState()
{
    LinearSolver::deallocateSolverState();

    // Deallocate any scratch data.
}

} // namespace multiphase
