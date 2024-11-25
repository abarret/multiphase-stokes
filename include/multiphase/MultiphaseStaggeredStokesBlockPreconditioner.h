#ifndef included_MultiphaseStaggeredStokesBlockPreconditioner
#define included_MultiphaseStaggeredStokesBlockPreconditioner

#include <ibtk/LinearSolver.h>

#include <multiphase/MultiphaseLSCSchurComplementSolver.h>
#include <multiphase/MultiphaseStokesBlockSolver.h>
#include <ibtk/PoissonSolver.h>

namespace multiphase
{
class MultiphaseStaggeredStokesBlockPreconditioner : IBTK::LinearSolver
{
public:
    MultiphaseStaggeredStokesBlockPreconditioner(std::string object_name);

    virtual ~MultiphaseStaggeredStokesBlockPreconditioner();

    /*!
     * Solve Mx = b with M the approximate schur complement preconditioner with
     * M = [F, G; 0, -S]
     *
     * First we solve S*xp = G^T*G(G^T*F*G)^-1(G^T*G)*xp = -bp.
     *
     * Then we solve
     * F*xu = bu - G*xp.
     *
     * Note that u consists of two velocity fields
     */
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * Initialize the solver state. Setup the sub solvers, allocate patch data, and set up communication algorithms.
     */
    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * Deallocate all temporary data used to solve the linear system.
     */
    void deallocateSolverState() override;

private:
    // Velocity solver settings.
    std::unique_ptr<MultiphaseStokesBlockSolver> d_stokes_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_solver_db;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_precond_db;

    // Pressure solver settings.
    std::unique_ptr<IBTK::PoissonSolver> d_pressure_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_pressure_solver_db;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_pressure_precond_db;

    // Hierarchy data
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_thn_ths_sq_var;
    int d_thn_ths_sq_idx = IBTK::invalid_index;

    int d_thn_idx = IBTK::invalid_index;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_thn_bc_coefs = nullptr;
};
} // namespace multiphase
#endif
