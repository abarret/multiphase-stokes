#ifndef included_MultiphaseStaggeredStokesBlockPreconditioner
#define included_MultiphaseStaggeredStokesBlockPreconditioner

#include <ibtk/LinearSolver.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PoissonSolver.h>

#include <multiphase/FullFACPreconditioner.h>
#include <multiphase/MultiphaseLSCSchurComplementSolver.h>
#include <multiphase/MultiphaseParameters.h>
#include <multiphase/MultiphaseStaggeredStokesBlockFACOperator.h>
#include <multiphase/MultiphaseStokesBlockSolver.h>

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
    // Velocity and Pressure vectors
    SAMRAI::solv::SAMRAIVectorReal<NDIM, double> d_u_vec, d_p_vec;
    SAMRAI::solv::SAMRAIVectorReal<NDIM, double> d_bu_vec, d_bp_vec;

    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double>> d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double>> d_hier_sc_data_ops;

    // Velocity solver settings.
    std::unique_ptr<IBTK::PETScKrylovLinearSolver> d_stokes_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_solver_db;
    SAMRAI::tbox::Pointer<FullFACPreconditioner> d_stokes_precond;
    SAMRAI::tbox::Pointer<MultiphaseStaggeredStokesBlockFACOperator> d_stokes_precond_op;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_precond_db;
    MultiphaseParameters d_params;

    // Pressure solver settings.
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_pressure_solver;
    std::string d_pressure_solver_type = "PETSC_KRYLOV_SOLVER";
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_pressure_solver_db;
    std::string d_pressure_precond_type = "POINT_RELAXATION_FAC_PRECONDITIONER";
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_pressure_precond_db;
    SAMRAI::solv::PoissonSpecifications d_pressure_coefs;

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
