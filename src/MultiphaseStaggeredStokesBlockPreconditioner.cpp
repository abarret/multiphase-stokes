#include <ibtk/app_namespaces.h>

#include <multiphase/MultiphaseStaggeredStokesBlockPreconditioner.h>

#include <SAMRAIVectorReal.h>

namespace multiphase
{

static constexpr int UN_COMP_IDX = 0;
static constexpr int US_COMP_IDX = 1;
static constexpr int P_COMP_IDX = 2;

MultiphaseStaggeredStokesBlockPreconditioner::MultiphaseStaggeredStokesBlockPreconditioner(std::string object_name)
    : d_object_name(std::move(object_name))
{
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

    // Setup any needed communication algorithms and scratch indices
}

void
MultiphaseStaggeredStokesBlockPreconditioner::deallocateSolverState()
{
    LinearSolver::deallocateSolverState();

    // Deallocate any scratch data.
}

} // namespace multiphase
