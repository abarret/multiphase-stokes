#ifndef included_multiphase_MultiphaseStandardHierarchyIntegrator
#define included_multiphase_MultiphaseStandardHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "multiphase/FullFACPreconditioner.h"
#include "multiphase/MultiphaseConvectiveManager.h"
#include "multiphase/MultiphaseStaggeredHierarchyIntegrator.h"
#include "multiphase/MultiphaseStaggeredStokesBoxRelaxationFACOperator.h"
#include "multiphase/MultiphaseStaggeredStokesOperator.h"
#include "multiphase/utility_functions.h"

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
 * \brief Class MultiphaseStandardHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations on an AMR grid hierarchy.
 *
 * This class will regularize the volume fraction during the application of the Stokes operators. The regularization
 * parameter can be set in the input database. Note that no regularization occurs in the momentum transport component.
 */
class MultiphaseStandardHierarchyIntegrator : public MultiphaseStaggeredHierarchyIntegrator
{
public:
    /*!
     * The constructor for class MultiphaseStandardHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    MultiphaseStandardHierarchyIntegrator(std::string object_name,
                                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                          bool register_for_restart = true);

    /*!
     * The destructor for class MultiphaseStandardHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~MultiphaseStandardHierarchyIntegrator();

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

    void regridHierarchyBeginSpecialized() override;
    void initializeCompositeHierarchyDataSpecialized(double init_data_time, bool initial_time) override;

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
     * Adds all body forces to the provided SAMRAI vector. Assumes first component of f_vec is the network forces and
     * the second ocmponent is the solvent forces.
     *
     * The volume fraction index thn_idx must have at least one layer of ghost cells to account for scaling the forces.
     */
    void addBodyForces(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>& f_vec,
                       double current_time,
                       double new_time,
                       int thn_cur_idx,
                       int thn_half_idx,
                       int thn_new_idx);

    /*!
     * Vectors
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_sol_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_rhs_vec;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>> d_nul_vecs;
    bool d_has_vel_nullspace = false, d_has_pressure_nullspace = true;

    /*!
     * Solver information
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_solver_db;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> d_stokes_solver;
    SAMRAI::tbox::Pointer<MultiphaseStaggeredStokesOperator> d_stokes_op;

    /*!
     * Preconditioner information
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_precond_db;
    SAMRAI::tbox::Pointer<FullFACPreconditioner> d_stokes_precond;
    SAMRAI::tbox::Pointer<MultiphaseStaggeredStokesBoxRelaxationFACOperator> d_precond_op;
    double d_w = std::numeric_limits<double>::quiet_NaN();
    bool d_use_preconditioner = true;

    bool d_make_div_rhs_sum_to_zero = true;

    /*!
     * Convective operator
     */
    std::unique_ptr<MultiphaseConvectiveManager> d_convec_op;
    IBAMR::LimiterType d_convec_limiter_type = IBAMR::LimiterType::UNKNOWN_LIMITER_TYPE;

    /*!
     * AB2 patch index for convective operator
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_Nn_old_var, d_Ns_old_var;

    /*!
     * BDF2 patch data. Note that for second order accuracy, we need
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_un_old_var, d_us_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_fn_old_var, d_fs_old_var;

    TimeSteppingType d_viscous_ts_type = TimeSteppingType::TRAPEZOIDAL_RULE;
};
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_MultiphaseStandardHierarchyIntegrator
