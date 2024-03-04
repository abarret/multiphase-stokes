#ifndef included_multiphase_VCTwoFluidStaggeredStokesBoxRelaxationFACOperator
#define included_multiphase_VCTwoFluidStaggeredStokesBoxRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PoissonFACPreconditioner.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"

#include "IntVector.h"
#include "PoissonSpecifications.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include <map>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace multiphase
{
class VCTwoFluidStaggeredStokesBoxRelaxationFACOperator : public IBTK::FACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     * \param w under relaxation factor in box relaxation scheme
     * \param C scaler-valued C in C*u term used to add diagonal dominance
     */
    VCTwoFluidStaggeredStokesBoxRelaxationFACOperator(const std::string& object_name,
                                                      const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~VCTwoFluidStaggeredStokesBoxRelaxationFACOperator();

    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef);

    /*!
     * \name Functions for configuring the solver.
     */
    //\{

    /*!
     * \brief Zero-out the provided vector on the specified level of the patch
     * hierarchy.
     */
    void setToZero(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error, int level_num) override;

    /*!
     * \brief Restrict the residual from the source vector to the destination
     * vector on the specified level of the patch hierarchy.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    void restrictResidual(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                          int dest_level_num) override;

    /*!
     * \brief Prolong the error from the source vector to the destination
     * vector on the specified level of the patch hierarchy.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    void prolongError(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                      int dest_level_num) override;

    /*!
     * \brief Prolong the error from the source vector to the destination vector
     * on the specified level of the patch hierarchy and correct the fine-grid
     * error.
     *
     * \note Implementations must support the case in which source and dest are
     * the same vector.
     */
    void prolongErrorAndCorrect(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                                int dest_level_num) override;

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are being performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are being performed
     */

    void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps) override;

    /*! \brief This method sets up the patch data indices for \param Thn */
    void setThnIdx(int thn_idx);

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */

    bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                            int coarsest_ln) override;

    /*!
     * \brief Compute composite grid residual on a range of levels.
     *
     * \param residual residual vector
     * \param solution solution vector
     * \param rhs source (right hand side) vector
     * \param coarsest_level_num coarsest level number
     * \param finest_level_num finest level number
     */
    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num) override;

    //\}
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs) override;

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void deallocateOperatorState() override;

    /*!
     * \brief Set the under relaxation parameter.
     */
    void setUnderRelaxationParamater(double w);

    /*!
     * \brief Set the C and D coefficients.
     */
    void setCandDCoefficients(double C, double D);

    /*!
     * \brief Set the viscosity coefficients for the viscous stresses.
     */
    void setViscosityCoefficient(double eta_n, double eta_s);

    /*!
     * \brief Set the drag coefficients for each phase.
     */
    void setDragCoefficient(double xi, double nu_n, double nu_s);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCTwoFluidStaggeredStokesBoxRelaxationFACOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCTwoFluidStaggeredStokesBoxRelaxationFACOperator(const VCTwoFluidStaggeredStokesBoxRelaxationFACOperator& from) =
        delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCTwoFluidStaggeredStokesBoxRelaxationFACOperator&
    operator=(const VCTwoFluidStaggeredStokesBoxRelaxationFACOperator& that) = delete;

    /*!
     * \brief Perform prolongation or restriction on the provided indices.
     */
    //\{
    void performProlongation(const std::array<int, 3>& dst_idxs, const std::array<int, 3>& src_idxs, int dst_ln);
    void performRestriction(const std::array<int, 3>& dst_idxs, const std::array<int, 3>& src_idxs, int dst_ln);
    void performGhostFilling(const std::array<int, 3>& dst_idxs, int dst_ln);
    //\}

    /*
     * Level solvers and solver parameters.
     */
    std::string d_level_solver_type = IBTK::CCPoissonSolverManager::PETSC_LEVEL_SOLVER,
                d_level_solver_default_options_prefix;
    double d_level_solver_abs_residual_tol = 1.0e-50, d_level_solver_rel_residual_tol = 1.0e-5;
    int d_level_solver_max_iterations = 1;
    std::vector<SAMRAI::tbox::Pointer<IBTK::PoissonSolver>> d_level_solvers;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_level_solver_db;

    /*
     * Coarse level solvers and solver parameters.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_coarse_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_coarse_solver_db;

    /*
     * Patch overlap data.
     */
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM>>> d_patch_bc_box_overlap;

    /*
     * Coarse-fine interface interpolation objects.
     */
    SAMRAI::tbox::Pointer<IBTK::CoarseFineBoundaryRefinePatchStrategy> d_sc_bdry_op, d_cc_bdry_op;

    // Cache the prolongation and restriction schedules. Note we also cache the algorithms so that we can reset the
    // schedules to their previous state.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM>> d_un_prolong_op, d_us_prolong_op, d_p_prolong_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM>> d_un_restrict_op, d_us_restrict_op, d_p_restrict_op;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> d_prolong_scheds;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> d_prolong_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM>>> d_restrict_scheds;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM>> d_restrict_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> d_ghostfill_no_restrict_scheds;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> d_ghostfill_no_restrict_alg;

    int d_thn_idx = IBTK::invalid_index;
    int d_un_scr_idx = IBTK::invalid_index, d_us_scr_idx = IBTK::invalid_index, d_p_scr_idx = IBTK::invalid_index;
    double d_w = 0.75;                                         // under relaxation factor
    double d_C = std::numeric_limits<double>::quiet_NaN();     // C*u
    double d_D = std::numeric_limits<double>::quiet_NaN();     // D depends on time stepping scheme
    double d_eta_n = std::numeric_limits<double>::quiet_NaN(); // Network viscosity
    double d_eta_s = std::numeric_limits<double>::quiet_NaN(); // Solvent viscosity
    double d_xi = std::numeric_limits<double>::quiet_NaN();    // Drag coefficient
    double d_nu_n = std::numeric_limits<double>::quiet_NaN();
    double d_nu_s = std::numeric_limits<double>::quiet_NaN();

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy; // Reference patch hierarchy
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_un_bc_coefs;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_us_bc_coefs;
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>> d_default_un_bc_coef, d_default_us_bc_coef,
        d_default_P_bc_coef, d_default_thn_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef = nullptr;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_thn_bc_coef = nullptr;
    SAMRAI::tbox::Pointer<IBTK::CartSideRobinPhysBdryOp> d_un_bc_op, d_us_bc_op;
    SAMRAI::tbox::Pointer<IBTK::CartCellRobinPhysBdryOp> d_P_bc_op, d_thn_bc_op;
    std::unique_ptr<SAMRAI::xfer::RefinePatchStrategy<NDIM>> d_vel_P_bc_op;
    // Cached communications operators.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> d_un_fill_pattern, d_us_fill_pattern,
        d_P_fill_pattern;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill, d_no_fill;

    std::vector<std::vector<std::array<SAMRAI::hier::BoxList<NDIM>, NDIM>>> d_patch_side_bc_box_overlap;
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM>>> d_patch_cell_bc_box_overlap;

    SAMRAI::tbox::Pointer<IBTK::StaggeredPhysicalBoundaryHelper> d_bc_un_helper, d_bc_us_helper;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int>> d_mask_var;
    int d_mask_idx = IBTK::invalid_index;
};
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_VCTwoFluidStaggeredStokesBoxRelaxationFACOperator
