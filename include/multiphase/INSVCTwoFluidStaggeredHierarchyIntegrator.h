#ifndef included_multiphase_INSVCTwoFluidStaggeredHierarchyIntegrator
#define included_multiphase_INSVCTwoFluidStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "multiphase/FullFACPreconditioner.h"
#include "multiphase/INSVCTwoFluidConvectiveManager.h"
#include "multiphase/VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.h"
#include "multiphase/VCTwoFluidStaggeredStokesOperator.h"

#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/SideDataSynchronization.h"
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/muParserCartGridFunction.h>

#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace multiphase
{
/*!
 * \brief Class INSVCTwoFluidStaggeredHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations on an AMR grid hierarchy.
 */
class INSVCTwoFluidStaggeredHierarchyIntegrator : public IBAMR::INSHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSVCTwoFluidStaggeredHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSVCTwoFluidStaggeredHierarchyIntegrator(std::string object_name,
                                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                              std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> u_dummy_coefs,
                                              SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_dummy_coefs,
                                              bool register_for_restart = true);

    /*!
     * The destructor for class INSVCTwoFluidStaggeredHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSVCTwoFluidStaggeredHierarchyIntegrator();

    /*!
     * Non-zero Reynolds numbers are not implemented. Returns nullptr.
     */
    SAMRAI::tbox::Pointer<IBAMR::ConvectiveOperator> getConvectiveOperator() override;

    /*!
     * Not in use. Returns nullptr.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getVelocitySubdomainSolver() override;

    /*!
     * Not in use. Returns nullptr.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getPressureSubdomainSolver() override;

    /*!
     * Get solvent velocity variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> getSolventVariable() const;

    /*!
     * Get network velocity variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> getNetworkVariable() const;

    /*!
     * Get pressure velocity variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> getPressureVariable() const;

    inline SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> getNetworkVolumeFractionVariable() const
    {
        return d_thn_cc_var;
    }

    inline SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> getNetworkBodyForceVariable() const
    {
        return d_f_un_sc_var;
    }

    inline SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> getSolventBodyForceVariable() const
    {
        return d_f_us_sc_var;
    }

    /*!
     * \brief Set the viscosity coefficients for the viscous stresses.
     */
    void setViscosityCoefficient(double eta_n, double eta_s);

    /*!
     * \brief Set the drag coefficients for each phase.
     */
    void setDragCoefficient(double xi, double nu_n, double nu_s);

    /*!
     * Set initial conditions for the state variables.
     *
     * NOTE: These pointers are set to nullptr after they are used.
     */
    void setInitialData(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> un_fcn,
                        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> us_fcn,
                        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> p_fcn);

    /*!
     * Set the initial conditions for the network volume fraction.
     */
    void setInitialNetworkVolumeFraction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> thn_init_fcn);

    /*!
     * Register forcing functions. These are not scaled by the volume fraction. Any can be NULL.
     */
    void setForcingFunctions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> fn_fcn,
                             SAMRAI::tbox::Pointer<IBTK::CartGridFunction> fs_fcn,
                             SAMRAI::tbox::Pointer<IBTK::CartGridFunction> fp_fcn = nullptr);

    /*!
     * Regiser forcing functions that should be scaled by the volume fraction.
     */
    void setForcingFunctionsScaled(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> fn_fcn,
                                   SAMRAI::tbox::Pointer<IBTK::CartGridFunction> fs_fcn);

    /*!
     * Register forcing functions that should be scaled by both volume fractions.
     */
    void setForcingFunctionsScaledByBoth(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> fn_fcn,
                                         SAMRAI::tbox::Pointer<IBTK::CartGridFunction> fs_fcn);

    /*!
     * Register the volume fraction function. If a function is not registered, the volume fraction is advected with the
     * fluid.
     *
     * An optional argument allows this function to be used as the initial condition for the volume fraction. If this is
     * set to false, users are expected to call setInitialNetworkVolumeFraction.
     */
    void setNetworkVolumeFractionFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> thn_fcn,
                                          bool use_as_initial_condition = true);

    /*!
     * Sets that the network volume fraction should be advected with the specified advection diffusion integrator. Also
     * registers the integrator with the hierarchy.
     *
     * The optional input has_meaningful_mid_value should be set to false if there is no meaningful "new" time value
     * after a single evolution. In particular, this should be set to false if using semi-Lagrangian methods to evolve
     * the volume fraction.
     */
    void advectNetworkVolumeFraction(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_integrator,
                                     bool has_meaningful_mid_value = true);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override;

    /*!
     * Virtual method to perform implementation-specific data initialization on
     * a new level after it is inserted into an AMR patch hierarchy by the
     * gridding algorithm.
     *
     * An empty default implementation is provided.
     */
    void initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM>> old_level,
                                        bool allocate_data) override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchySpecialized(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    void regridProjection() override;

    double getStableTimestep(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch) const override;

    /*!
     * Reset cached hierarchy dependent data.
     */
    void resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                                int coarsest_level,
                                                int finest_level) override;

    void applyGradientDetectorSpecialized(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                          int level_num,
                                          double error_data_time,
                                          int tag_idx,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too);

    void regridHierarchyBeginSpecialized() override;
    void initializeCompositeHierarchyDataSpecialized(double init_data_time, bool initial_time) override;

protected:
    void setupPlotDataSpecialized() override;

    /*!
     * Return the maximum time step size as determined by https://doi.org/10.1006/jcph.1998.5967
     */
    double getMaximumTimeStepSizeSpecialized() override;

private:
    void approxConvecOp(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>& f_vec,
                        double current_time,
                        double new_time,
                        int un_cur_idx,
                        int us_cur_idx,
                        int thn_cur_idx,
                        int un_new_idx,
                        int us_new_idx,
                        int thn_new_idx);

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCTwoFluidStaggeredHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCTwoFluidStaggeredHierarchyIntegrator(const INSVCTwoFluidStaggeredHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCTwoFluidStaggeredHierarchyIntegrator&
    operator=(const INSVCTwoFluidStaggeredHierarchyIntegrator& that) = delete;

    /*!
     * \brief A helper function to set the volume fractions. Sets thn_scr_idx = 0.5 * (thn_cur_idx + thn_new_idx).
     *
     * If start_of_ts = true, then this function will recalculate thn at current and new times if a function is
     * prescribed. If thn is advected, then we set thn_new_idx to thn_cur_idx.
     */
    void setThnAtHalf(int& thn_cur_idx,
                      int& thn_new_idx,
                      int& thn_scr_idx,
                      double current_time,
                      double new_time,
                      bool start_of_ts = false);

    /*!
     * Adds all body forces to the provided SAMRAI vector. Assumes first component of f_vec is the network forces and
     * the second ocmponent is the solvent forces.
     *
     * The volume fraction index thn_idx must have at least one layer of ghost cells to account for scaling the forces.
     */
    void addBodyForces(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> f_vec,
                       double eval_time,
                       int thn_idx);

    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_f_un_fcn, d_f_us_fcn, d_f_p_fcn;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_f_un_thn_fcn, d_f_us_ths_fcn;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_f_un_thn_ths_fcn, d_f_us_thn_ths_fcn;

    // CartGridFunctions that set the initial values for state variables.
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_un_init_fcn, d_us_init_fcn, d_p_init_fcn, d_thn_init_fcn;

    /*!
     * Fluid solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_un_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_us_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_un_rhs_var, d_us_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_p_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_f_un_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_f_us_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_f_cc_var;

    /*!
     * Gradient detector variables
     */
    bool d_use_grad_tagging = false;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_grad_thn_var;
    SAMRAI::tbox::Array<double> d_abs_grad_thresh, d_rel_grad_thresh;
    double d_max_grad_thn = std::numeric_limits<double>::quiet_NaN();

    /*!
     * Vectors
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_sol_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_rhs_vec;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>> d_nul_vecs;
    bool d_has_vel_nullspace = false;

    // Physical parameters
    double d_rho = std::numeric_limits<double>::quiet_NaN();
    double d_xi = std::numeric_limits<double>::quiet_NaN(), d_nu_n = std::numeric_limits<double>::quiet_NaN(),
           d_nu_s = std::numeric_limits<double>::quiet_NaN();
    double d_eta_n = std::numeric_limits<double>::quiet_NaN(), d_eta_s = std::numeric_limits<double>::quiet_NaN();

    /*!
     * Solver information
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_solver_db;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> d_stokes_solver;
    SAMRAI::tbox::Pointer<VCTwoFluidStaggeredStokesOperator> d_stokes_op;

    /*!
     * Preconditioner information
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_precond_db;
    SAMRAI::tbox::Pointer<FullFACPreconditioner> d_stokes_precond;
    SAMRAI::tbox::Pointer<VCTwoFluidStaggeredStokesBoxRelaxationFACOperator> d_precond_op;
    double d_w = std::numeric_limits<double>::quiet_NaN();
    bool d_use_preconditioner = true;

    /*!
     * Drawing information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_un_draw_var, d_us_draw_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_div_draw_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_fn_draw_var, d_fs_draw_var;

    /*!
     * Objects that can do the operations we need
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double>> d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double>> d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double>> d_hier_fc_data_ops;

    /*!
     * Network volume fraction. We either specify the function thn_fcn, or we use the d_thn_integrator.
     *
     * IMPORTANT: d_thn_integrator should be listed in d_adv_diff_hier_integrators. We ONLY use d_thn_integrator for the
     * contexts.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_thn_cc_var;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_thn_fcn;
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_thn_integrator;

    bool d_make_div_rhs_sum_to_zero = true;

    /*!
     * Convective operator
     */
    std::unique_ptr<INSVCTwoFluidConvectiveManager> d_convec_op;
    IBAMR::LimiterType d_convec_limiter_type = IBAMR::LimiterType::UNKNOWN_LIMITER_TYPE;

    /*!
     * AB2 patch index for convective operator
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_Nn_old_var, d_Ns_old_var;

    // Scratch force index
    int d_fn_scr_idx = IBTK::invalid_index, d_fs_scr_idx = IBTK::invalid_index;

    // Flag for if volume fraction has a meaningful value in middle of time step.
    bool d_use_new_thn = true;

    bool d_use_accel_ts = false;
    double d_accel_ts_safety_fac = 0.1;
};
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSVCTwoFluidStaggeredHierarchyIntegrator
