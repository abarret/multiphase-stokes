/////////////////////////////// INCLUDES /////////////////////////////////////

#include "multiphase/FullFACPreconditioner.h"
#include "multiphase/MultiphaseStaggeredStokesBoxRelaxationFACOperator.h"
#include "multiphase/MultiphaseStaggeredStokesOperator.h"
#include "multiphase/MultiphaseStandardHierarchyIntegrator.h"
#include "multiphase/utility_functions.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/INSStaggeredConvectiveOperatorManager.h"
#include "ibamr/INSStaggeredPressureBcCoef.h"
#include "ibamr/INSStaggeredVelocityBcCoef.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include <ibamr/INSVCStaggeredPressureBcCoef.h>
#include <ibamr/INSVCStaggeredVelocityBcCoef.h>

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellDoubleBoundsPreservingConservativeLinearRefine.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/CartSideDoubleSpecializedLinearRefine.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/ibtk_utilities.h"
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/PETScKrylovLinearSolver.h>

#include "ArrayData.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace multiphase
{

namespace
{
static Timer* t_integrate_hierarchy = nullptr;
static Timer* t_preprocess_integrate_hierarchy = nullptr;
static Timer* t_postprocess_integrate_hierarchy = nullptr;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////
MultiphaseStandardHierarchyIntegrator::MultiphaseStandardHierarchyIntegrator(std::string object_name,
                                                                             Pointer<Database> input_db,
                                                                             bool register_for_restart)
    : MultiphaseStaggeredHierarchyIntegrator(std::move(object_name), input_db, register_for_restart)
{
    // We do not use d_viscous_time_stepping_type. Use d_viscous_ts_type instead (allows for BDF2).
    d_viscous_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
    if (input_db->keyExists("viscous_time_stepping_type"))
        d_viscous_ts_type = string_to_enum<TimeSteppingType>(input_db->getString("viscous_time_stepping_type"));
    if (input_db->keyExists("solver_db")) d_solver_db = input_db->getDatabase("solver_db");
    if (input_db->keyExists("precond_db")) d_precond_db = input_db->getDatabase("precond_db");
    if (input_db->keyExists("use_preconditioner")) d_use_preconditioner = input_db->getBool("use_preconditioner");
    if (input_db->keyExists("precond_type") && d_use_preconditioner) d_precond_type = multiphase::string_to_enum<PreconditionerType>(input_db->getString("precond_type"));
    if (input_db->keyExists("make_div_rhs_sum_to_zero"))
        d_make_div_rhs_sum_to_zero = input_db->getBool("make_div_rhs_sum_to_zero");
    d_has_vel_nullspace = input_db->getBoolWithDefault("has_vel_nullspace", d_has_vel_nullspace);
    d_has_pressure_nullspace = input_db->getBoolWithDefault("has_pressure_nullspace", d_has_pressure_nullspace);
    d_convec_limiter_type = IBAMR::string_to_enum<LimiterType>(
        input_db->getStringWithDefault("convec_limiter_type", IBAMR::enum_to_string(d_convec_limiter_type)));
    d_convective_time_stepping_type = IBAMR::string_to_enum<IBAMR::TimeSteppingType>(
        input_db->getStringWithDefault("convec_ts_type", "FORWARD_EULER"));

    // Make sure viscous time stepping type is valid for this class
    if (d_viscous_ts_type != TimeSteppingType::BACKWARD_EULER &&
        d_viscous_ts_type != TimeSteppingType::TRAPEZOIDAL_RULE && d_viscous_ts_type != TimeSteppingType::BDF2)
        TBOX_ERROR(d_object_name + ": Viscous time step type " +
                   IBAMR::enum_to_string<TimeSteppingType>(d_viscous_ts_type) +
                   " not valid. Must use BACKWARD_EULER, TRAPEZOIDAL_RULE, or BDF2");

    // Note we only support convective operating type of AB2 and forward Euler if we are using BDF for viscosity
    if (d_viscous_ts_type == TimeSteppingType::BDF2)
    {
        if (d_convective_time_stepping_type != ADAMS_BASHFORTH && d_convective_time_stepping_type != FORWARD_EULER)
        {
            TBOX_ERROR(d_object_name +
                           ": Viscous operator is set as BDF2. We only support convective operators of Adams Bashforth "
                           "or forward Euler. Choice of "
                       << IBAMR::enum_to_string(d_convective_time_stepping_type) << " is invalid.\n");
        }
    }

    IBTK_DO_ONCE(t_integrate_hierarchy = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStandardHierarchyIntegrator::integrateHierarchy()");
                 t_preprocess_integrate_hierarchy = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStandardHierarchyIntegrator::preprocess()");
                 t_postprocess_integrate_hierarchy = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStandardHierarchyIntegrator::postprocess()"););
    return;
} // MultiphaseStandardHierarchyIntegrator

MultiphaseStandardHierarchyIntegrator::~MultiphaseStandardHierarchyIntegrator()
{
    // intentionally blank
} // ~MultiphaseStandardHierarchyIntegrator

void
MultiphaseStandardHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                                     Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized) return;

    MultiphaseStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Here we do all we need to ensure that calls to advanceHierarchy() or integrateHierarchy() are valid.
    // NOTE: This function is called before the patch hierarchy has valid patch levels.
    // To set initial data, we should do this in initializeLevelDataSpecialized().

    // Note we only register AB2 variables if we are using that time stepping routine
    if (d_convective_time_stepping_type == ADAMS_BASHFORTH)
    {
        d_Nn_old_var = new SideVariable<NDIM, double>(d_object_name + "::Nn");
        d_Ns_old_var = new SideVariable<NDIM, double>(d_object_name + "::Ns");
        int Nn_idx, Ns_idx;
        registerVariable(Nn_idx, d_Nn_old_var, 0, getCurrentContext());
        registerVariable(Ns_idx, d_Ns_old_var, 0, getCurrentContext());
    }

    // We only register old velocity variables if we are using BDF2
    if (d_viscous_ts_type == TimeSteppingType::BDF2)
    {
        d_un_old_var = new SideVariable<NDIM, double>(d_object_name + "::un_old");
        d_us_old_var = new SideVariable<NDIM, double>(d_object_name + "::us_old");
        d_fn_old_var = new SideVariable<NDIM, double>(d_object_name + "::fn_old");
        d_fs_old_var = new SideVariable<NDIM, double>(d_object_name + "::fs_old");
        int un_old_idx, us_old_idx, fn_old_idx, fs_old_idx;
        registerVariable(un_old_idx, d_un_old_var, 0, getCurrentContext());
        registerVariable(us_old_idx, d_us_old_var, 0, getCurrentContext());
        registerVariable(fn_old_idx, d_fn_old_var, 0, getCurrentContext());
        registerVariable(fs_old_idx, d_fs_old_var, 0, getCurrentContext());
    }

    // Note that we MUST create the convective operators here because they require large ghost cell widths. The maximum
    // ghost cell width MUST be specified before the patch hierarchy is finished being created. Note that because we
    // also create the convective operator in ::initializeCompositeHierarchyDataSpecialized(), this object will be
    // reset.
    if (!d_creeping_flow)
        d_convec_op = std::make_unique<MultiphaseConvectiveManager>(d_object_name + "::ConvectiveOp",
                                                                    d_hierarchy,
                                                                    d_convec_limiter_type,
                                                                    d_un_bc_coefs,
                                                                    d_us_bc_coefs,
                                                                    d_thn_bc_coef);

    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
MultiphaseStandardHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
    if (!d_creeping_flow) d_convec_op.reset();
}

void
MultiphaseStandardHierarchyIntegrator::initializeCompositeHierarchyDataSpecialized(const double init_data_time,
                                                                                   const bool initial_time)
{
    // Set up the convective operator
    if (!d_creeping_flow)
        d_convec_op = std::make_unique<MultiphaseConvectiveManager>(d_object_name + "::ConvectiveOp",
                                                                    d_hierarchy,
                                                                    d_convec_limiter_type,
                                                                    d_un_bc_coefs,
                                                                    d_us_bc_coefs,
                                                                    d_thn_bc_coef);
}

void
MultiphaseStandardHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                    const double new_time,
                                                                    const int num_cycles)
{
    IBTK_TIMER_START(t_preprocess_integrate_hierarchy);
    // Do anything that needs to be done before we call integrateHierarchy().
    MultiphaseStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const double half_time = 0.5 * (new_time + current_time);

    // Pull out current solution components
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_cur_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getCurrentContext());
    const int us_cur_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getCurrentContext());
    const int p_cur_idx = var_db->mapVariableAndContextToIndex(d_P_var, getCurrentContext());
    const int un_new_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getNewContext());
    const int us_new_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getNewContext());
    const int p_new_idx = var_db->mapVariableAndContextToIndex(d_P_var, getNewContext());
    const int un_scr_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getScratchContext());
    const int us_scr_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getScratchContext());
    const int p_scr_idx = var_db->mapVariableAndContextToIndex(d_P_var, getScratchContext());

    // Allocate scratch and new data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Set initial guess
    d_hier_sc_data_ops->copyData(un_new_idx, un_cur_idx);
    d_hier_sc_data_ops->copyData(us_new_idx, us_cur_idx);
    d_hier_cc_data_ops->copyData(p_new_idx, p_cur_idx);
    d_hier_sc_data_ops->copyData(un_scr_idx, un_cur_idx);
    d_hier_sc_data_ops->copyData(us_scr_idx, us_cur_idx);
    d_hier_cc_data_ops->copyData(p_scr_idx, p_cur_idx);

    // Set up our solution vector
    d_sol_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::sol_vec", d_hierarchy, coarsest_ln, finest_ln);

    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    d_sol_vec->addComponent(d_un_sc_var, un_scr_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_sol_vec->addComponent(d_us_sc_var, us_scr_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_sol_vec->addComponent(d_P_var, p_scr_idx, wgt_cc_idx, d_hier_cc_data_ops);

    // Set up the RHS vector
    const int un_rhs_idx = var_db->mapVariableAndContextToIndex(d_un_rhs_var, getScratchContext());
    const int us_rhs_idx = var_db->mapVariableAndContextToIndex(d_us_rhs_var, getScratchContext());
    const int p_rhs_idx = var_db->mapVariableAndContextToIndex(d_p_rhs_var, getScratchContext());
    d_rhs_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_rhs_vec->addComponent(d_un_rhs_var, un_rhs_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_rhs_vec->addComponent(d_us_rhs_var, us_rhs_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_rhs_vec->addComponent(d_p_rhs_var, p_rhs_idx, wgt_cc_idx, d_hier_cc_data_ops);
    d_rhs_vec->setToScalar(0.0);

    // Grab the theta index
    int thn_cur_idx, thn_new_idx, thn_scr_idx;
    setThnAtHalf(thn_cur_idx, thn_new_idx, thn_scr_idx, current_time, new_time, /*start_of_ts*/ true);

    // Set up null vectors (if applicable). Note that "d_has_vel_nullspace" corresponds to cases when rho=0 and drag
    // coefficient != 0. If the drag coefficient is zero, there are additional elements in the nullspace.
    int num_null_vecs = (d_has_pressure_nullspace ? 1 : 0) + (d_has_vel_nullspace ? NDIM : 0);
    d_nul_vecs.resize(num_null_vecs);
    for (size_t i = 0; i < d_nul_vecs.size(); ++i)
    {
        d_nul_vecs[i] = d_sol_vec->cloneVector("NullVec_" + std::to_string(i)); // should delete the vector at the end
        d_nul_vecs[i]->allocateVectorData();
        d_nul_vecs[i]->setToScalar(0.0);
    }
    // Pull out pressure component and set to constant
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            if (d_has_pressure_nullspace)
            {
                Pointer<CellData<NDIM, double>> p_data = d_nul_vecs[0]->getComponentPatchData(2, *patch);
                p_data->fillAll(1.0);
            }

            if (d_has_vel_nullspace)
            {
                int shft = d_has_pressure_nullspace ? 1 : 0;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    Pointer<SideData<NDIM, double>> un_data = d_nul_vecs[axis + shft]->getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double>> us_data = d_nul_vecs[axis + shft]->getComponentPatchData(1, *patch);
                    un_data->getArrayData(axis).fillAll(1.0);
                    us_data->getArrayData(axis).fillAll(1.0);
                }
            }
        }
    }

    // set-up RHS to treat viscosity and drag with backward Euler or Implicit Trapezoidal Rule:
    // RHS = f(n) + C*theta_i(n)*u_i(n) + D1*(viscous + drag) for  i = n, s
    double D1 = std::numeric_limits<double>::signaling_NaN();
    double D2 = std::numeric_limits<double>::signaling_NaN();
    double C = d_params.rho / dt;

    switch (d_viscous_ts_type)
    {
    case TimeSteppingType::TRAPEZOIDAL_RULE:
        D1 = 0.5;
        D2 = -0.5;
        break;

    case TimeSteppingType::BACKWARD_EULER:
    case TimeSteppingType::BDF2:
        D1 = 0.0;
        D2 = -1.0;
        break;

    default:
        TBOX_ERROR("Unknown time stepping type " + IBAMR::enum_to_string<TimeSteppingType>(d_viscous_ts_type) +
                   ". Valid options are BACKWARD_EULER, TRAPEZOIDAL_RULE, or BDF2.");
    }

    // Set drag coefficient if necessary
    if (isVariableDrag())
    {
        const double eval_time = half_time;
        const int xi_idx = var_db->mapVariableAndContextToIndex(d_xi_var, getScratchContext());
        d_xi_fcn->setDataOnPatchHierarchy(xi_idx, d_xi_var, d_hierarchy, eval_time, false, coarsest_ln, finest_ln);
    }
    // Boundary Condition helper.
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_un_helper = new StaggeredStokesPhysicalBoundaryHelper();
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_us_helper = new StaggeredStokesPhysicalBoundaryHelper();

    if (getIntegratorStep() != 0 && d_viscous_ts_type == TimeSteppingType::BDF2)
    {
        // If the integrator step is 0 (so initial time step), we reduce to backward Euler, and this block is skipped.
        // Note that the old velocities are already multiplied by volume fraction.
        const int un_old_idx = var_db->mapVariableAndContextToIndex(d_un_old_var, getCurrentContext());
        const int us_old_idx = var_db->mapVariableAndContextToIndex(d_us_old_var, getCurrentContext());

        const int rhs_un_idx = d_rhs_vec->getComponentDescriptorIndex(0);
        const int rhs_us_idx = d_rhs_vec->getComponentDescriptorIndex(1);

        double alpha = d_dt_previous[0] / dt;
        C = d_params.rho * (2.0 + alpha) / (dt * (1.0 + alpha));

        d_hier_sc_data_ops->scale(rhs_un_idx, -d_params.rho / (dt * alpha * (alpha + 1.0)), un_old_idx);
        d_hier_sc_data_ops->scale(rhs_us_idx, -d_params.rho / (dt * alpha * (alpha + 1.0)), us_old_idx);

        multiply_sc_and_thn(d_fn_scr_idx, un_cur_idx, thn_cur_idx, d_hierarchy);
        multiply_sc_and_ths(d_fs_scr_idx, us_cur_idx, thn_cur_idx, d_hierarchy);

        d_hier_sc_data_ops->axpy(rhs_un_idx, d_params.rho * (alpha + 1.0) / (dt * alpha), d_fn_scr_idx, rhs_un_idx);
        d_hier_sc_data_ops->axpy(rhs_us_idx, d_params.rho * (alpha + 1.0) / (dt * alpha), d_fs_scr_idx, rhs_us_idx);
    }
    else
    {
        MultiphaseStaggeredStokesOperator RHS_op("RHS_op", true, d_params, d_thn_cur_manager);
        RHS_op.setPhysicalBoundaryHelper(bc_un_helper, bc_us_helper);
        RHS_op.setPhysicalBcCoefs(d_un_bc_coefs, d_us_bc_coefs, d_p_bc_coef);
        RHS_op.setSolutionTime(current_time);
        RHS_op.setTimeInterval(current_time, new_time);
        // Divergence free condition and pressure are not time stepped. We do not need to account for the contributions
        // in the RHS.
        RHS_op.setCandDCoefficients(C, D1, 0.0, 0.0);

        // Store results of applying stokes operator in rhs_vec
        RHS_op.initializeOperatorState(*d_sol_vec, *d_rhs_vec);
        RHS_op.apply(*d_sol_vec, *d_rhs_vec);
    }

    // Set up the operators and solvers needed to solve the linear system.
    d_stokes_op = new MultiphaseStaggeredStokesOperator("stokes_op", false, d_params, d_thn_new_manager);
    d_stokes_op->setPhysicalBcCoefs(d_un_bc_coefs, d_us_bc_coefs, d_p_bc_coef);
    d_stokes_op->setCandDCoefficients(C, D2);
    d_stokes_op->setPhysicalBoundaryHelper(bc_un_helper, bc_us_helper);
    d_stokes_op->setSolutionTime(new_time);
    d_stokes_op->setTimeInterval(current_time, new_time);

    d_stokes_solver = new PETScKrylovLinearSolver("solver", d_solver_db, "solver_");
    d_stokes_solver->setOperator(d_stokes_op);

    // Now create a preconditioner
    if (d_use_preconditioner)
    {
        d_precond = std::make_unique<MultiphasePreconditioner>(d_precond_type,
                                                               d_hierarchy,
                                                               d_precond_db,
                                                               d_params,
                                                               C,
                                                               D2,
                                                               d_thn_new_manager,
                                                               d_nul_vecs,
                                                               d_un_bc_coefs,
                                                               d_us_bc_coefs,
                                                               d_p_bc_coef);
        d_stokes_solver->setPreconditioner(d_precond->getPreconditioner());
    }

    d_stokes_solver->setSolutionTime(new_time);
    d_stokes_solver->setTimeInterval(current_time, new_time);
    d_stokes_solver->setNullSpace(false, d_nul_vecs);
    d_stokes_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec);

    // Set thn_cc_idx on the dense hierarchy.
    if (d_use_preconditioner)
    {
        if (d_precond->updateVolumeFraction(d_thn_new_manager, new_time))
        {
            d_stokes_solver->deallocateSolverState();
            d_stokes_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec);
        }
    }

    // Set the forcing data if applicable. We only do this for pressure. The momentum forces are applied in
    // integrateHierarchy().
    const int f_p_idx = var_db->mapVariableAndContextToIndex(d_f_cc_var, getScratchContext());
    // Divergence is not time-stepped. We need to prescribe the correct divergence at the time point for which we are
    // solving.
    if (d_f_p_fcn)
    {
        d_f_p_fcn->setDataOnPatchHierarchy(f_p_idx, d_f_cc_var, d_hierarchy, new_time);
        d_hier_cc_data_ops->add(
            d_rhs_vec->getComponentDescriptorIndex(2), f_p_idx, d_rhs_vec->getComponentDescriptorIndex(2));
    }

    // Set up the advection diffusion integrator
    for (const auto& adv_diff_integrator : d_adv_diff_hier_integrators)
    {
        const int adv_diff_num_cycles = adv_diff_integrator->getNumberOfCycles();
        // Network advection velocity
        const int U_adv_diff_cur_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_integrator->getCurrentContext());
        const int U_adv_diff_scr_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_integrator->getScratchContext());
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_cur_idx)) copy_side_to_face(U_adv_diff_cur_idx, un_cur_idx, d_hierarchy);
        const int us_adv_diff_cur_idx =
            var_db->mapVariableAndContextToIndex(d_us_adv_diff_var, adv_diff_integrator->getCurrentContext());
        const int us_adv_diff_scr_idx =
            var_db->mapVariableAndContextToIndex(d_us_adv_diff_var, adv_diff_integrator->getScratchContext());
        const int us_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_us_adv_diff_var, adv_diff_integrator->getNewContext());
        if (isAllocatedPatchData(us_adv_diff_cur_idx)) copy_side_to_face(us_adv_diff_cur_idx, us_cur_idx, d_hierarchy);

        adv_diff_integrator->preprocessIntegrateHierarchy(current_time, new_time, adv_diff_num_cycles);

        if (isAllocatedPatchData(U_adv_diff_scr_idx))
            d_hier_fc_data_ops->copyData(U_adv_diff_scr_idx, U_adv_diff_cur_idx);
        if (isAllocatedPatchData(U_adv_diff_new_idx))
            d_hier_fc_data_ops->copyData(U_adv_diff_new_idx, U_adv_diff_cur_idx);
        if (isAllocatedPatchData(us_adv_diff_scr_idx))
            d_hier_fc_data_ops->copyData(us_adv_diff_scr_idx, us_adv_diff_cur_idx);
        if (isAllocatedPatchData(us_adv_diff_new_idx))
            d_hier_fc_data_ops->copyData(us_adv_diff_new_idx, us_adv_diff_cur_idx);
    }

    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    IBTK_TIMER_STOP(t_preprocess_integrate_hierarchy);
    return;
} // preprocessIntegrateHierarchy

void
MultiphaseStandardHierarchyIntegrator::integrateHierarchySpecialized(const double current_time,
                                                                     const double new_time,
                                                                     const int cycle_num)
{
    IBTK_TIMER_START(t_integrate_hierarchy);
    MultiphaseStaggeredHierarchyIntegrator::integrateHierarchySpecialized(current_time, new_time, cycle_num);
    double half_time = 0.5 * (current_time + new_time);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_cur_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getCurrentContext());
    const int us_cur_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getCurrentContext());
    const int un_new_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getNewContext());
    const int us_new_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getNewContext());
    const int p_new_idx = var_db->mapVariableAndContextToIndex(d_P_var, getNewContext());
    // Update the state of the advection diffusion integrator
    for (const auto& adv_diff_integrator : d_adv_diff_hier_integrators)
        adv_diff_integrator->integrateHierarchy(current_time, new_time, cycle_num);

    // Determine new values of volume fraction.
    int thn_cur_idx, thn_new_idx, thn_scr_idx;
    setThnAtHalf(thn_cur_idx, thn_new_idx, thn_scr_idx, current_time, new_time, /*start_of_ts*/ !d_use_new_thn);

    // Update the preconditioner with new volume fraction
    if (d_precond)
    {
        bool needs_reallocation = d_precond->updateVolumeFraction(d_thn_new_manager, new_time);
        if (isVariableDrag())
        {
            const int xi_idx = var_db->mapVariableAndContextToIndex(d_xi_var, getScratchContext());
            needs_reallocation =
                needs_reallocation || d_precond->updateDragCoefficient(xi_idx, d_xi_var, new_time, d_xi_fcn);
        }
        if (needs_reallocation)
        {
            d_stokes_solver->deallocateSolverState();
            d_stokes_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec);
        }
    }
    const double dt = new_time - current_time;

    // Compute forces
    Pointer<SAMRAIVectorReal<NDIM, double>> f_vec = d_rhs_vec->cloneVector(d_object_name + "::F_temp");
    f_vec->allocateVectorData(current_time);
    f_vec->setToScalar(0.0);

    addBodyForces(f_vec, current_time, new_time, thn_cur_idx, thn_scr_idx, thn_new_idx);
    if (!d_creeping_flow)
    {
        approxConvecOp(
            f_vec, current_time, new_time, un_cur_idx, us_cur_idx, thn_cur_idx, un_new_idx, us_new_idx, thn_new_idx);
    }

    d_rhs_vec->add(d_rhs_vec, f_vec);

    // Compute weighted sum for rhs divergence.
    // TODO: This needs a permanent home, rather than an if statement.
    if (d_make_div_rhs_sum_to_zero)
    {
        const int rhs_p_idx = d_rhs_vec->getComponentDescriptorIndex(2);
        HierarchyCellDataOpsReal<NDIM, double> hier_sc_data_ops(d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
        double integral = hier_sc_data_ops.integral(rhs_p_idx, d_hier_math_ops->getCellWeightPatchDescriptorIndex());
        hier_sc_data_ops.addScalar(rhs_p_idx, rhs_p_idx, -1.0 * integral);
        integral = hier_sc_data_ops.integral(rhs_p_idx, d_hier_math_ops->getCellWeightPatchDescriptorIndex());
    }

    // Set the initial guess for the system to be the most recent approximation to t^{n+1}
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(0), un_new_idx);
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(1), us_new_idx);
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(2), p_new_idx);

    // Synchronize the rhs before we solve
    {
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp(2);
        ghost_cell_comp[0] = ITC(d_rhs_vec->getComponentDescriptorIndex(0), "NONE", false, "CONSERVATIVE_COARSEN");
        ghost_cell_comp[1] = ITC(d_rhs_vec->getComponentDescriptorIndex(1), "NONE", false, "CONSERVATIVE_COARSEN");
        HierarchyGhostCellInterpolation hier_ghost_fill;
        hier_ghost_fill.initializeOperatorState(ghost_cell_comp, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
        hier_ghost_fill.fillData(current_time);
    }

    // Solve for un(n+1), us(n+1), p(n+1).
    bool converged = d_stokes_solver->solveSystem(*d_sol_vec, *d_rhs_vec);
    if (d_enable_logging)
        pout << "Stokes solver " << (converged ? "converged" : "failed to converge") << " after "
             << d_stokes_solver->getNumIterations() << " iterations\n";

    if (d_normalize_pressure)
    {
        for (const auto& nul_vec : d_nul_vecs)
        {
            const double sol_dot_nul = d_sol_vec->dot(nul_vec);
            const double nul_L2_norm = std::sqrt(nul_vec->dot(nul_vec));
            d_sol_vec->axpy(-sol_dot_nul / nul_L2_norm, nul_vec, d_sol_vec);
        }
    }

    // Reset the solve vector to copy the "scratch" data into the "new" data
    d_hier_sc_data_ops->copyData(un_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_sc_data_ops->copyData(us_new_idx, d_sol_vec->getComponentDescriptorIndex(1));
    d_hier_cc_data_ops->copyData(p_new_idx, d_sol_vec->getComponentDescriptorIndex(2));

    // Reset the RHS
    d_rhs_vec->subtract(d_rhs_vec, f_vec);

    // Destroy forces
    f_vec->deallocateVectorData();
    f_vec->freeVectorComponents();

    // Do any more updates to the advection diffusion variables
    for (const auto& adv_diff_hier_integrator : d_adv_diff_hier_integrators)
    {
        // Reset the velocities
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_new_idx)) copy_side_to_face(U_adv_diff_new_idx, un_new_idx, d_hierarchy);
        const int U_adv_diff_cur_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_hier_integrator->getCurrentContext());
        const int U_adv_diff_scr_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_hier_integrator->getScratchContext());
        if (isAllocatedPatchData(U_adv_diff_scr_idx))
            d_hier_fc_data_ops->linearSum(U_adv_diff_scr_idx, 0.5, U_adv_diff_new_idx, 0.5, U_adv_diff_cur_idx);

        const int us_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_us_adv_diff_var, adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(us_adv_diff_new_idx)) copy_side_to_face(us_adv_diff_new_idx, us_new_idx, d_hierarchy);
        const int us_adv_diff_cur_idx =
            var_db->mapVariableAndContextToIndex(d_us_adv_diff_var, adv_diff_hier_integrator->getCurrentContext());
        const int us_adv_diff_scr_idx =
            var_db->mapVariableAndContextToIndex(d_us_adv_diff_var, adv_diff_hier_integrator->getScratchContext());
        if (isAllocatedPatchData(us_adv_diff_scr_idx))
            d_hier_fc_data_ops->linearSum(us_adv_diff_scr_idx, 0.5, us_adv_diff_new_idx, 0.5, us_adv_diff_cur_idx);

        // Now update the state variables. Note that cycle 0 has already been performed
        const int adv_diff_num_cycles = adv_diff_hier_integrator->getNumberOfCycles();
        if (d_current_num_cycles != adv_diff_num_cycles)
        {
            for (int adv_diff_cycle_num = 1; adv_diff_cycle_num < adv_diff_num_cycles; ++adv_diff_cycle_num)
                adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, adv_diff_cycle_num);
        }
    }
    IBTK_TIMER_STOP(t_integrate_hierarchy);
    return;
} // integrateHierarchy

void
MultiphaseStandardHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                     const double new_time,
                                                                     const bool skip_synchronize_new_state_data,
                                                                     const int num_cycles)
{
    IBTK_TIMER_START(t_postprocess_integrate_hierarchy);

    // Replace N_old data if necessary
    if (!d_creeping_flow && is_multistep_time_stepping_type(d_convective_time_stepping_type))
    {
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int Nn_old_idx = var_db->mapVariableAndContextToIndex(d_Nn_old_var, getCurrentContext());
        const int Ns_old_idx = var_db->mapVariableAndContextToIndex(d_Ns_old_var, getCurrentContext());
        d_convec_op->fillWithConvectiveOperator(Nn_old_idx, Ns_old_idx);
    }

    if (d_viscous_ts_type == TimeSteppingType::BDF2)
    {
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int un_old_idx = var_db->mapVariableAndContextToIndex(d_un_old_var, getCurrentContext());
        const int us_old_idx = var_db->mapVariableAndContextToIndex(d_us_old_var, getCurrentContext());
        const int fn_old_idx = var_db->mapVariableAndContextToIndex(d_fn_old_var, getCurrentContext());
        const int fs_old_idx = var_db->mapVariableAndContextToIndex(d_fs_old_var, getCurrentContext());
        const int un_cur_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getCurrentContext());
        const int us_cur_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getCurrentContext());
        const int thn_cur_idx = d_thn_cur_manager->getCellIndex();
        const int fn_cur_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, getCurrentContext());
        const int fs_cur_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, getCurrentContext());

        multiply_sc_and_thn(un_old_idx, un_cur_idx, thn_cur_idx, d_hierarchy);
        multiply_sc_and_ths(us_old_idx, us_cur_idx, thn_cur_idx, d_hierarchy);
        d_hier_sc_data_ops->copyData(fn_old_idx, fn_cur_idx);
        d_hier_sc_data_ops->copyData(fs_old_idx, fs_cur_idx);
    }

    // Do anything that needs to be done after integrateHierarchy().
    MultiphaseStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Note: Preconditioner should deallocate data as necessary.
    // Deallocate scratch data
    for (const auto& adv_diff_hier_integrator : d_adv_diff_hier_integrators)
    {
        const int adv_diff_num_cycles = adv_diff_hier_integrator->getNumberOfCycles();
        adv_diff_hier_integrator->postprocessIntegrateHierarchy(
            current_time, new_time, skip_synchronize_new_state_data, adv_diff_num_cycles);
    }

    for (auto& nul_vec : d_nul_vecs)
    {
        if (nul_vec)
        {
            nul_vec->deallocateVectorData();
            nul_vec->freeVectorComponents();
        }
    }
    d_nul_vecs.clear();

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    d_stokes_op = nullptr;
    d_stokes_solver = nullptr;
    d_precond = nullptr;
    IBTK_TIMER_STOP(t_postprocess_integrate_hierarchy);
    return;
} // postprocessIntegrateHierarchy

void
MultiphaseStandardHierarchyIntegrator::approxConvecOp(Pointer<SAMRAIVectorReal<NDIM, double>>& f_vec,
                                                      double current_time,
                                                      double new_time,
                                                      int un_cur_idx,
                                                      int us_cur_idx,
                                                      int thn_cur_idx,
                                                      int un_new_idx,
                                                      int us_new_idx,
                                                      int thn_new_idx)
{
    // Note that MultiphaseConvectiveOperator is not smart enough to handle multistep time steppers. Therefore, we
    // have to do something special if we are using a multistep algorithm.
    IBAMR::TimeSteppingType ts_type = d_convective_time_stepping_type;
    if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(ts_type))
        ts_type = d_init_convective_time_stepping_type;
    if (ts_type != ADAMS_BASHFORTH)
    {
        d_convec_op->approximateConvectiveOperator(d_fn_scr_idx,
                                                   d_fs_scr_idx,
                                                   ts_type,
                                                   current_time,
                                                   new_time,
                                                   un_cur_idx,
                                                   us_cur_idx,
                                                   thn_cur_idx,
                                                   un_new_idx,
                                                   us_new_idx,
                                                   thn_new_idx);
    }
    else
    {
        d_convec_op->approximateConvectiveOperator(d_fn_scr_idx,
                                                   d_fs_scr_idx,
                                                   FORWARD_EULER,
                                                   current_time,
                                                   new_time,
                                                   un_cur_idx,
                                                   us_cur_idx,
                                                   thn_cur_idx,
                                                   un_new_idx,
                                                   us_new_idx,
                                                   thn_new_idx);
        // Now add previous time step. Note if we are using BDF2 for viscosity, we require different coefficients.
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int Nn_old_idx = var_db->mapVariableAndContextToIndex(d_Nn_old_var, getCurrentContext());
        const int Ns_old_idx = var_db->mapVariableAndContextToIndex(d_Ns_old_var, getCurrentContext());
        const double dt = new_time - current_time;
        const double alpha = d_dt_previous[0] / dt;
        if (d_viscous_ts_type == TimeSteppingType::BDF2)
        {
            d_hier_sc_data_ops->linearSum(d_fn_scr_idx, 1.0 + 1.0 / alpha, d_fn_scr_idx, -1.0 / alpha, Nn_old_idx);
            d_hier_sc_data_ops->linearSum(d_fs_scr_idx, 1.0 + 1.0 / alpha, d_fs_scr_idx, -1.0 / alpha, Ns_old_idx);
        }
        else
        {
            d_hier_sc_data_ops->linearSum(d_fn_scr_idx, 1.0 + 0.5 / alpha, d_fn_scr_idx, -0.5 / alpha, Nn_old_idx);
            d_hier_sc_data_ops->linearSum(d_fs_scr_idx, 1.0 + 0.5 / alpha, d_fs_scr_idx, -0.5 / alpha, Ns_old_idx);
        }
    }
    d_hier_sc_data_ops->linearSum(
        f_vec->getComponentDescriptorIndex(0), 1.0, f_vec->getComponentDescriptorIndex(0), -d_params.rho, d_fn_scr_idx);
    d_hier_sc_data_ops->linearSum(
        f_vec->getComponentDescriptorIndex(1), 1.0, f_vec->getComponentDescriptorIndex(1), -d_params.rho, d_fs_scr_idx);
}

void
MultiphaseStandardHierarchyIntegrator::addBodyForces(Pointer<SAMRAIVectorReal<NDIM, double>>& f_vec,
                                                     const double current_time,
                                                     const double new_time,
                                                     const int thn_cur_idx,
                                                     const int thn_half_idx,
                                                     const int thn_new_idx)
{
    // If we are using BDF2, then we evaluate the forces using AB2
    TimeSteppingType ts_type = d_viscous_ts_type;
    if (getIntegratorStep() == 0) ts_type = TimeSteppingType::BACKWARD_EULER;
    double eval_time = ts_type == TimeSteppingType::BDF2 ? current_time : 0.5 * (current_time + new_time);
    double dt = new_time - current_time;
    double alpha = ts_type == TimeSteppingType::BDF2 ? d_dt_previous[0] / dt : 1.0;
    const int thn_idx = ts_type == TimeSteppingType::BDF2 ? thn_cur_idx : thn_half_idx;

    // Create some scratch indices for accumulation
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int fn_cloned_idx = var_db->registerClonedPatchDataIndex(d_f_un_sc_var, d_fn_scr_idx);
    const int fs_cloned_idx = var_db->registerClonedPatchDataIndex(d_f_us_sc_var, d_fs_scr_idx);
    allocate_patch_data(
        { fn_cloned_idx, fs_cloned_idx }, d_hierarchy, eval_time, 0, d_hierarchy->getFinestLevelNumber());
    d_hier_sc_data_ops->setToScalar(d_fn_scr_idx, 0.0);
    d_hier_sc_data_ops->setToScalar(d_fs_scr_idx, 0.0);

    // Start with forces that do not need to be scaled.
    if (d_f_un_fcn)
    {
        d_f_un_fcn->setDataOnPatchHierarchy(fn_cloned_idx, d_f_un_sc_var, d_hierarchy, eval_time);
        d_hier_sc_data_ops->add(d_fn_scr_idx, fn_cloned_idx, d_fn_scr_idx);
    }
    if (d_f_us_fcn)
    {
        d_f_us_fcn->setDataOnPatchHierarchy(fs_cloned_idx, d_f_us_sc_var, d_hierarchy, eval_time);
        d_hier_sc_data_ops->add(d_fs_scr_idx, fs_cloned_idx, d_fs_scr_idx);
    }

    // Now account for scaled forces. Note volume fractions should always have ghost cells filled
    if (d_f_un_thn_fcn)
    {
        d_f_un_thn_fcn->setDataOnPatchHierarchy(fn_cloned_idx, d_f_un_sc_var, d_hierarchy, eval_time);
        multiply_sc_and_thn(fn_cloned_idx, fn_cloned_idx, thn_idx, d_hierarchy);
        d_hier_sc_data_ops->add(d_fn_scr_idx, fn_cloned_idx, d_fn_scr_idx);
    }
    if (d_f_us_ths_fcn)
    {
        d_f_us_ths_fcn->setDataOnPatchHierarchy(fs_cloned_idx, d_f_us_sc_var, d_hierarchy, eval_time);
        multiply_sc_and_ths(fs_cloned_idx, fs_cloned_idx, thn_idx, d_hierarchy);
        d_hier_sc_data_ops->add(d_fs_scr_idx, fs_cloned_idx, d_fs_scr_idx);
    }

    if (d_f_un_thn_ths_fcn)
    {
        d_f_un_thn_ths_fcn->setDataOnPatchHierarchy(fn_cloned_idx, d_f_un_sc_var, d_hierarchy, eval_time);
        multiply_sc_and_thn_and_ths(fn_cloned_idx, fn_cloned_idx, thn_idx, d_hierarchy);
        d_hier_sc_data_ops->add(d_fn_scr_idx, d_fn_scr_idx, fn_cloned_idx);
    }
    if (d_f_us_thn_ths_fcn)
    {
        d_f_us_thn_ths_fcn->setDataOnPatchHierarchy(fs_cloned_idx, d_f_us_sc_var, d_hierarchy, eval_time);
        multiply_sc_and_thn_and_ths(fs_cloned_idx, fs_cloned_idx, thn_idx, d_hierarchy);
        d_hier_sc_data_ops->add(d_fs_scr_idx, d_fs_scr_idx, fs_cloned_idx);
    }

    // Copy our accumulated forces into a permanent location
    const int fn_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, getCurrentContext());
    const int fs_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, getCurrentContext());
    d_hier_sc_data_ops->copyData(fn_idx, d_fn_scr_idx);
    d_hier_sc_data_ops->copyData(fs_idx, d_fs_scr_idx);

    if (ts_type == TimeSteppingType::BDF2)
    {
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int fn_old_idx = var_db->mapVariableAndContextToIndex(d_fn_old_var, getCurrentContext());
        const int fs_old_idx = var_db->mapVariableAndContextToIndex(d_fs_old_var, getCurrentContext());
        // Note old forces should already by scaled by volume fraction.
        d_hier_sc_data_ops->axpy(
            f_vec->getComponentDescriptorIndex(0), -1.0 / alpha, fn_old_idx, f_vec->getComponentDescriptorIndex(0));
        d_hier_sc_data_ops->axpy(
            f_vec->getComponentDescriptorIndex(1), -1.0 / alpha, fs_old_idx, f_vec->getComponentDescriptorIndex(1));
        d_hier_sc_data_ops->axpy(f_vec->getComponentDescriptorIndex(0),
                                 1.0 + 1.0 / alpha,
                                 d_fn_scr_idx,
                                 f_vec->getComponentDescriptorIndex(0));
        d_hier_sc_data_ops->axpy(f_vec->getComponentDescriptorIndex(1),
                                 1.0 + 1.0 / alpha,
                                 d_fs_scr_idx,
                                 f_vec->getComponentDescriptorIndex(1));
    }
    else
    {
        d_hier_sc_data_ops->add(
            f_vec->getComponentDescriptorIndex(0), d_fn_scr_idx, f_vec->getComponentDescriptorIndex(0));
        d_hier_sc_data_ops->add(
            f_vec->getComponentDescriptorIndex(1), d_fs_scr_idx, f_vec->getComponentDescriptorIndex(1));
    }

    // Remove cloned patch indices
    deallocate_patch_data({ fn_cloned_idx, fs_cloned_idx }, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    var_db->removePatchDataIndex(fn_cloned_idx);
    var_db->removePatchDataIndex(fs_cloned_idx);
}

MultiphaseStandardHierarchyIntegrator::MultiphasePreconditioner::MultiphasePreconditioner(
    PreconditionerType precond_type,
    Pointer<PatchHierarchy<NDIM>> hierarchy,
    Pointer<Database> input_db,
    const MultiphaseParameters& params,
    const double C,
    const double D,
    const std::unique_ptr<VolumeFractionDataManager>& thn_manager,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>& null_vecs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
    RobinBcCoefStrategy<NDIM>* p_bc_coef)
    : d_hierarchy(hierarchy), d_precond_type(precond_type)
{
    switch (d_precond_type)
    {
    case PreconditionerType::BLOCK:
    {
        d_block_precond =
            new MultiphaseStaggeredStokesBlockPreconditioner("KrylovPrecondStrategy", params, input_db, thn_manager);
        d_block_precond->setCAndDCoefficients(C, D);
        d_block_precond->setNullSpace(false, null_vecs);
        break;
    }
    case PreconditionerType::MULTIGRID:
    {
        d_fac_op = new MultiphaseStaggeredStokesBoxRelaxationFACOperator(
            "KrylovPrecondStrategy", "Krylov_precond_", params, thn_manager);
        d_fac_op->setUnderRelaxationParamater(input_db->getDouble("relaxation_parameter"));
        d_fac_op->setCandDCoefficients(C, D);
        d_fac_op->setPhysicalBcCoefs(un_bc_coefs, us_bc_coefs, p_bc_coef);
        d_fac_precond = new FullFACPreconditioner("KrylovPrecond", d_fac_op, input_db, "Krylov_precond_");
        d_fac_precond->setNullSpace(false, null_vecs);
        break;
    }
    default:
    {
        TBOX_ERROR("Unknown preconditioner type!\n");
        break;
    }
    }
}

bool
MultiphaseStandardHierarchyIntegrator::MultiphasePreconditioner::updateVolumeFraction(
    const std::unique_ptr<VolumeFractionDataManager>& thn_manager,
    const double time)
{
    switch (d_precond_type)
    {
    case PreconditionerType::BLOCK:
    {
        if (d_block_precond->getVolumeFractionManager() == thn_manager)
        {
            d_block_precond->deallocateSolverState();
        }
        else
        {
            TBOX_ERROR(
                "MultiphaseStandardHierarchyIntegrator::MultiphasePreconditioner::updateVolumeFraction(): Can't hot "
                "swap volume fraction managers!\n");
        }

        return true;
        break;
    }
    case PreconditionerType::MULTIGRID:
    {
        Pointer<FACPreconditionerStrategy> fac_strat = d_fac_precond->getFACPreconditionerStrategy();
        Pointer<MultiphaseStaggeredStokesBoxRelaxationFACOperator> box_relax_op = fac_strat;
        if (box_relax_op)
        {
            if (box_relax_op->getVolumeFractionManager() == thn_manager)
            {
                Pointer<PatchHierarchy<NDIM>> dense_hierarchy = d_fac_precond->getDenseHierarchy();
                d_fac_precond->transferToDense(thn_manager->getCellIndex());
                thn_manager->updateVolumeFraction(thn_manager->getCellIndex(), dense_hierarchy, time);
            }
            else
            {
                TBOX_ERROR(
                    "MultiphaseStandardHierarchyIntegrator::MultiphasePreconditioner::updateVolumeFraction(): Can't "
                    "hot swap volume fraction managers!\n");
            }
        }
        return false;
        break;
    }
    default:
    {
        TBOX_ERROR("Unknown preconditioner type!\n");
        return false;
    }
    }
}

bool
MultiphaseStandardHierarchyIntegrator::MultiphasePreconditioner::updateDragCoefficient(
    const int drag_idx,
    Pointer<SideVariable<NDIM, double>>& drag_var,
    const double time,
    IBTK::CartGridFunction* drag_fcn)
{
    switch (d_precond_type)
    {
    case PreconditionerType::BLOCK:
    {
        TBOX_ERROR("BLOCK preconditioner not set up for variable drag!\n");
        return false;
        break;
    }
    case PreconditionerType::MULTIGRID:
    {
        if (drag_fcn) drag_fcn->setDataOnPatchHierarchy(drag_idx, drag_var, d_hierarchy, time);
        d_fac_precond->transferToDense(drag_idx, true);
        return false;
        break;
    }
    default:
    {
        TBOX_ERROR("Unknown preconditioner type!\n");
        return true;
    }
    }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
