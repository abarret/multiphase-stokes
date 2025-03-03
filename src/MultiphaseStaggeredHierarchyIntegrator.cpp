/////////////////////////////// INCLUDES /////////////////////////////////////

#include "multiphase/FullFACPreconditioner.h"
#include "multiphase/MultiphaseStaggeredHierarchyIntegrator.h"
#include "multiphase/MultiphaseStaggeredStokesBlockPreconditioner.h"
#include "multiphase/MultiphaseStaggeredStokesBoxRelaxationFACOperator.h"
#include "multiphase/MultiphaseStaggeredStokesOperator.h"
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

/////////////////////////////// PUBLIC ///////////////////////////////////////
MultiphaseStaggeredHierarchyIntegrator::MultiphaseStaggeredHierarchyIntegrator(
    std::string object_name,
    Pointer<Database> input_db,
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry,
    bool register_for_restart)
    : INSHierarchyIntegrator(std::move(object_name),
                             input_db,
                             new SideVariable<NDIM, double>(object_name + "::Us"),
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             new CellVariable<NDIM, double>(object_name + "::P"),
                             "CONSERVATIVE_COARSEN",
                             "LINEAR_REFINE",
                             new SideVariable<NDIM, double>(object_name + "::F"),
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             new CellVariable<NDIM, double>(object_name + "::Q"),
                             "CONSERVATIVE_COARSEN",
                             "CONSTANT_REFINE",
                             register_for_restart)
{
    // We do not use d_viscous_time_stepping_type. Use d_viscous_ts_type instead (allows for BDF2).
    d_viscous_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
    if (input_db->keyExists("viscous_time_stepping_type"))
        d_viscous_ts_type = string_to_enum<TimeSteppingType>(input_db->getString("viscous_time_stepping_type"));
    if (input_db->keyExists("rho")) d_params.rho = input_db->getDouble("rho");
    if (input_db->keyExists("eta_n")) d_params.eta_n = input_db->getDouble("eta_n");
    if (input_db->keyExists("eta_s")) d_params.eta_s = input_db->getDouble("eta_s");
    if (input_db->keyExists("solver_db")) d_solver_db = input_db->getDatabase("solver_db");
    if (input_db->keyExists("precond_db")) d_precond_db = input_db->getDatabase("precond_db");
    if (input_db->keyExists("use_preconditioner")) d_use_preconditioner = input_db->getBool("use_preconditioner");
    if (input_db->keyExists("precond_type") && d_use_preconditioner) d_precond_type = multiphase::string_to_enum<PreconditionerType>(input_db->getString("precond_type"));
    if (input_db->keyExists("use_grad_tagging")) d_use_grad_tagging = input_db->getBool("use_grad_tagging");
    if (input_db->keyExists("grad_rel_thresh")) input_db->getArray("grad_rel_thresh", d_rel_grad_thresh);
    if (input_db->keyExists("grad_abs_thresh")) input_db->getArray("grad_abs_thresh", d_abs_grad_thresh);
    if (input_db->keyExists("make_div_rhs_sum_to_zero"))
        d_make_div_rhs_sum_to_zero = input_db->getBool("make_div_rhs_sum_to_zero");
    d_has_vel_nullspace = input_db->getBoolWithDefault("has_vel_nullspace", d_has_vel_nullspace);
    d_convec_limiter_type = IBAMR::string_to_enum<LimiterType>(
        input_db->getStringWithDefault("convec_limiter_type", IBAMR::enum_to_string(d_convec_limiter_type)));
    d_convective_time_stepping_type = IBAMR::string_to_enum<IBAMR::TimeSteppingType>(
        input_db->getStringWithDefault("convec_ts_type", "FORWARD_EULER"));

    d_input_db = input_db;
    // TODO: The default here should really be "false", but for now, this will not change the default behavior.
    d_creeping_flow = input_db->getBoolWithDefault("creeping_flow", true);
    d_un_sc_var = new SideVariable<NDIM, double>(d_object_name + "::un_sc");
    d_us_sc_var = new SideVariable<NDIM, double>(d_object_name + "::us_sc");

    // Arbitrarily set the base class velocity variable to the solvent variable. Note this is required for the base
    // class to compute the CFL number.
    d_U_var = d_us_sc_var;

    d_un_bc_coefs.resize(NDIM, nullptr);
    d_us_bc_coefs.resize(NDIM, nullptr);

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

    return;
} // MultiphaseStaggeredHierarchyIntegrator

MultiphaseStaggeredHierarchyIntegrator::~MultiphaseStaggeredHierarchyIntegrator()
{
    // intentionally blank
} // ~MultiphaseStaggeredHierarchyIntegrator

Pointer<ConvectiveOperator>
MultiphaseStaggeredHierarchyIntegrator::getConvectiveOperator()
{
    return nullptr;
} // getConvectiveOperator

Pointer<PoissonSolver>
MultiphaseStaggeredHierarchyIntegrator::getVelocitySubdomainSolver()
{
    return nullptr;
} // getVelocitySubdomainSolver

Pointer<PoissonSolver>
MultiphaseStaggeredHierarchyIntegrator::getPressureSubdomainSolver()
{
    return nullptr;
} // getPressureSubdomainSolver

Pointer<SideVariable<NDIM, double>>
MultiphaseStaggeredHierarchyIntegrator::getSolventVariable() const
{
    return d_us_sc_var;
}

Pointer<SideVariable<NDIM, double>>
MultiphaseStaggeredHierarchyIntegrator::getNetworkVariable() const
{
    return d_un_sc_var;
}

Pointer<CellVariable<NDIM, double>>
MultiphaseStaggeredHierarchyIntegrator::getPressureVariable() const
{
    return d_P_var;
}

void
MultiphaseStaggeredHierarchyIntegrator::registerPhysicalBoundaryConditions(
    std::vector<RobinBcCoefStrategy<NDIM>*> un_bc_coefs,
    std::vector<RobinBcCoefStrategy<NDIM>*> us_bc_coefs)
{
    d_un_bc_coefs = std::move(un_bc_coefs);
    d_us_bc_coefs = std::move(us_bc_coefs);
}

const std::vector<RobinBcCoefStrategy<NDIM>*>&
MultiphaseStaggeredHierarchyIntegrator::getNetworkBoundaryConditions() const
{
    return d_un_bc_coefs;
}

const std::vector<RobinBcCoefStrategy<NDIM>*>&
MultiphaseStaggeredHierarchyIntegrator::getSolventBoundaryConditions() const
{
    return d_us_bc_coefs;
}

void
MultiphaseStaggeredHierarchyIntegrator::registerVolumeFractionBoundaryConditions(RobinBcCoefStrategy<NDIM>* thn_bc_coef)
{
    d_thn_bc_coef = thn_bc_coef;
}

void
MultiphaseStaggeredHierarchyIntegrator::setViscosityCoefficient(const double eta_n, const double eta_s)
{
    d_params.eta_n = eta_n;
    d_params.eta_s = eta_s;
}

void
MultiphaseStaggeredHierarchyIntegrator::setDragCoefficient(const double xi, const double nu_n, const double nu_s)
{
    d_params.xi = xi;
    d_params.nu_n = nu_n;
    d_params.nu_s = nu_s;
}

void
MultiphaseStaggeredHierarchyIntegrator::setDragCoefficientFunction(Pointer<CartGridFunction> xi_fcn)
{
    d_xi_fcn = xi_fcn;
    if (!d_xi_var) d_xi_var = new SideVariable<NDIM, double>(d_object_name + "::Xi");
}

void
MultiphaseStaggeredHierarchyIntegrator::setInitialData(Pointer<CartGridFunction> un_fcn,
                                                       Pointer<CartGridFunction> us_fcn,
                                                       Pointer<CartGridFunction> p_fcn)
{
    d_un_init_fcn = un_fcn;
    d_us_init_fcn = us_fcn;
    d_p_init_fcn = p_fcn;
    return;
}

void
MultiphaseStaggeredHierarchyIntegrator::setForcingFunctions(Pointer<CartGridFunction> fn_fcn,
                                                            Pointer<CartGridFunction> fs_fcn,
                                                            Pointer<CartGridFunction> fp_fcn)
{
    if (fn_fcn)
    {
        if (d_f_un_fcn)
        {
            Pointer<CartGridFunctionSet> p_f_un_fcn = d_f_un_fcn;
            if (p_f_un_fcn)
            {
                p_f_un_fcn->addFunction(fn_fcn);
            }
            else
            {
                p_f_un_fcn = new CartGridFunctionSet(d_object_name + "::fn_fcn_set");
                p_f_un_fcn->addFunction(d_f_un_fcn);
                p_f_un_fcn->addFunction(fn_fcn);
                d_f_un_fcn = p_f_un_fcn;
            }
        }
        else
        {
            d_f_un_fcn = fn_fcn;
        }
    }

    if (fs_fcn)
    {
        if (d_f_us_fcn)
        {
            Pointer<CartGridFunctionSet> p_f_us_fcn = d_f_us_fcn;
            if (p_f_us_fcn)
            {
                p_f_us_fcn->addFunction(fs_fcn);
            }
            else
            {
                p_f_us_fcn = new CartGridFunctionSet(d_object_name + "::fs_fcn_set");
                p_f_us_fcn->addFunction(d_f_us_fcn);
                p_f_us_fcn->addFunction(fs_fcn);
                d_f_us_fcn = p_f_us_fcn;
            }
        }
        else
        {
            d_f_us_fcn = fs_fcn;
        }
    }

    if (fp_fcn)
    {
        if (d_f_p_fcn)
        {
            Pointer<CartGridFunctionSet> p_f_p_fcn = d_f_p_fcn;
            if (p_f_p_fcn)
            {
                p_f_p_fcn->addFunction(fp_fcn);
            }
            else
            {
                p_f_p_fcn = new CartGridFunctionSet(d_object_name + "::fp_fcn_set");
                p_f_p_fcn->addFunction(d_f_p_fcn);
                p_f_p_fcn->addFunction(fp_fcn);
                d_f_p_fcn = p_f_p_fcn;
            }
        }
        else
        {
            d_f_p_fcn = fp_fcn;
        }
    }
    return;
}

void
MultiphaseStaggeredHierarchyIntegrator::setForcingFunctionsScaled(Pointer<CartGridFunction> fn_fcn,
                                                                  Pointer<CartGridFunction> fs_fcn)
{
    if (fn_fcn)
    {
        if (d_f_un_thn_fcn)
        {
            Pointer<CartGridFunctionSet> p_f_un_fcn = d_f_un_thn_fcn;
            if (p_f_un_fcn)
            {
                p_f_un_fcn->addFunction(fn_fcn);
            }
            else
            {
                p_f_un_fcn = new CartGridFunctionSet(d_object_name + "::fn_thn_fcn_set");
                p_f_un_fcn->addFunction(d_f_un_thn_fcn);
                p_f_un_fcn->addFunction(fn_fcn);
                d_f_un_thn_fcn = p_f_un_fcn;
            }
        }
        else
        {
            d_f_un_thn_fcn = fn_fcn;
        }
    }

    if (fs_fcn)
    {
        if (d_f_us_ths_fcn)
        {
            Pointer<CartGridFunctionSet> p_f_us_fcn = d_f_us_ths_fcn;
            if (p_f_us_fcn)
            {
                p_f_us_fcn->addFunction(fs_fcn);
            }
            else
            {
                p_f_us_fcn = new CartGridFunctionSet(d_object_name + "::fs_ths_fcn_set");
                p_f_us_fcn->addFunction(d_f_us_ths_fcn);
                p_f_us_fcn->addFunction(fs_fcn);
                d_f_us_ths_fcn = p_f_us_fcn;
            }
        }
        else
        {
            d_f_us_ths_fcn = fs_fcn;
        }
    }
}

void
MultiphaseStaggeredHierarchyIntegrator::setInitialNetworkVolumeFraction(Pointer<CartGridFunction> thn_init_fcn)
{
    d_thn_init_fcn = thn_init_fcn;
}

void
MultiphaseStaggeredHierarchyIntegrator::setNetworkVolumeFractionFunction(Pointer<CartGridFunction> thn_fcn,
                                                                         bool use_as_initial_data)
{
    // Make sure we have a valid pointer.
    TBOX_ASSERT(thn_fcn);
    d_thn_fcn = thn_fcn;
    if (use_as_initial_data) d_thn_init_fcn = thn_fcn;
}

void
MultiphaseStaggeredHierarchyIntegrator::advectNetworkVolumeFraction(
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
    const bool has_meaningful_mid_value)
{
    registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
    d_thn_integrator = adv_diff_integrator;

    // Set up thn to be advected
    d_thn_cc_var = new CellVariable<NDIM, double>(d_object_name + "::thn_cc");
    d_thn_integrator->registerTransportedQuantity(d_thn_cc_var, true /*output_Q*/);
    d_thn_integrator->setAdvectionVelocity(d_thn_cc_var, d_U_adv_diff_var);
    d_thn_integrator->setAdvectionVelocityIsDivergenceFree(d_U_adv_diff_var, false); // Not divergence free in general.
    d_thn_integrator->setInitialConditions(d_thn_cc_var, d_thn_init_fcn);

    d_use_new_thn = has_meaningful_mid_value;
}

void
MultiphaseStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                                      Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Here we do all we need to ensure that calls to advanceHierarchy() or integrateHierarchy() are valid.
    // NOTE: This function is called before the patch hierarchy has valid patch levels.
    // To set initial data, we should do this in initializeLevelDataSpecialized().

    // First create the variables we need.
    // NOTE: d_P_var is a member variable of the base class.
    // NOTE: This may have been created above with the advection diffusion integrator.
    if (!d_thn_cc_var) d_thn_cc_var = new CellVariable<NDIM, double>(d_object_name + "::thn_cc");
    d_grad_thn_var = new CellVariable<NDIM, double>(d_object_name + "grad_thn_cc", NDIM);
    d_f_un_sc_var = new SideVariable<NDIM, double>(d_object_name + "::f_un_sc");
    d_f_us_sc_var = new SideVariable<NDIM, double>(d_object_name + "::f_us_sc");
    d_f_cc_var = new CellVariable<NDIM, double>(d_object_name + "::f_cc");
    d_un_rhs_var = new SideVariable<NDIM, double>(d_object_name + "::un_rhs");
    d_us_rhs_var = new SideVariable<NDIM, double>(d_object_name + "::us_rhs");
    d_p_rhs_var = new CellVariable<NDIM, double>(d_object_name + "::p_rhs");

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

    // Register variables with the integrator. Those with states (Velocities) have associated current, new, and scratch
    // indices. NOTE: Pressure is NOT a state variable, but we keep track of current and new values for initial guesses
    // for the solver.
    int un_cur_idx, un_scr_idx, un_new_idx;
    int us_cur_idx, us_scr_idx, us_new_idx;
    int p_cur_idx, p_scr_idx, p_new_idx;
    registerVariable(un_cur_idx,
                     un_new_idx,
                     un_scr_idx,
                     d_un_sc_var,
                     IntVector<NDIM>(1),
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_un_init_fcn);
    registerVariable(us_cur_idx,
                     us_new_idx,
                     us_scr_idx,
                     d_us_sc_var,
                     IntVector<NDIM>(1),
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_us_init_fcn);
    registerVariable(p_cur_idx,
                     p_new_idx,
                     p_scr_idx,
                     d_P_var,
                     IntVector<NDIM>(1),
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_p_init_fcn);

    // Everything else only gets a scratch context, which is deallocated at the end of each time step.
    // Note the forces need ghost cells for modifying the RHS to account for non-homogenous boundary conditions.
    int thn_cur_idx, thn_scr_idx, thn_new_idx, f_p_idx, f_un_idx, f_us_idx, un_rhs_idx, us_rhs_idx, p_rhs_idx, xi_idx;
    registerVariable(f_p_idx, d_f_cc_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(f_un_idx, d_f_un_sc_var, IntVector<NDIM>(1), getCurrentContext());
    registerVariable(f_us_idx, d_f_us_sc_var, IntVector<NDIM>(1), getCurrentContext());
    registerVariable(thn_cur_idx, d_thn_cc_var, IntVector<NDIM>(1), getCurrentContext());
    registerVariable(thn_scr_idx, d_thn_cc_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(thn_new_idx, d_thn_cc_var, IntVector<NDIM>(1), getNewContext());
    registerVariable(un_rhs_idx, d_un_rhs_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(us_rhs_idx, d_us_rhs_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(p_rhs_idx, d_p_rhs_var, IntVector<NDIM>(1), getScratchContext());
    if (d_xi_var)
    {
        registerVariable(xi_idx, d_xi_var, IntVector<NDIM>(0), getScratchContext());
        d_params.xi_idx = xi_idx;
    }

    // Register a scratch force object
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_fn_scr_idx = var_db->registerClonedPatchDataIndex(d_f_un_sc_var, f_un_idx);
    d_fs_scr_idx = var_db->registerClonedPatchDataIndex(d_f_us_sc_var, f_us_idx);
    d_scratch_data.setFlag(d_fn_scr_idx);
    d_scratch_data.setFlag(d_fs_scr_idx);

    // Register gradient of thn with the current context
    int grad_thn_idx;
    registerVariable(grad_thn_idx, d_grad_thn_var, IntVector<NDIM>(0), getCurrentContext());

    // Drawing variables. Velocities are node centered that get allocated with the current context.
    // Pressure is cell centered, so we can just use the current context for that.
    // Note: Using the current context for the velocities means drawing variables are allocated with the state
    // variables. This is important because the interface has no way of separately allocating and deallocating "drawing"
    // data.
    if (d_visit_writer)
    {
        d_un_draw_var = new NodeVariable<NDIM, double>(d_object_name + "::Un_draw", NDIM);
        d_us_draw_var = new NodeVariable<NDIM, double>(d_object_name + "::Us_draw", NDIM);
        d_div_draw_var = new CellVariable<NDIM, double>(d_object_name + "::Div_Draw");
        d_fn_draw_var = new CellVariable<NDIM, double>(d_object_name + "::Fn_draw", NDIM);
        d_fs_draw_var = new CellVariable<NDIM, double>(d_object_name + "::Fs_draw", NDIM);
        int un_draw_idx, us_draw_idx, div_draw_idx, fn_draw_idx, fs_draw_idx;
        registerVariable(un_draw_idx, d_un_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(us_draw_idx, d_us_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(div_draw_idx, d_div_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(fn_draw_idx, d_fn_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(fs_draw_idx, d_fs_draw_var, IntVector<NDIM>(0), getCurrentContext());

        d_visit_writer->registerPlotQuantity("Un", "VECTOR", un_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Un_" + std::to_string(d), "SCALAR", un_draw_idx, d, 1.0, "NODE");

        d_visit_writer->registerPlotQuantity("Us", "VECTOR", us_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Us_" + std::to_string(d), "SCALAR", us_draw_idx, d, 1.0, "NODE");

        d_visit_writer->registerPlotQuantity("P", "SCALAR", p_cur_idx, 0, 1.0, "CELL");

        d_visit_writer->registerPlotQuantity("Div", "SCALAR", div_draw_idx, 0, 1.0, "CELL");

        d_visit_writer->registerPlotQuantity("Fn", "VECTOR", fn_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Fn_" + std::to_string(d), "SCALAR", fn_draw_idx, d, 1.0, "CELL");

        d_visit_writer->registerPlotQuantity("Fs", "VECTOR", fs_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Fs_" + std::to_string(d), "SCALAR", fs_draw_idx, d, 1.0, "CELL");

        // Only need to plot this variable if we aren't advecting it.
        // If we do advect theta, the advection integrator will plot it.
        if (d_thn_fcn) d_visit_writer->registerPlotQuantity("Thn", "SCALAR", thn_cur_idx, 0, 1.0, "CELL");
    }

    // Set the correct boundary conditions. If we advect the volume fraction, then we overwrite this class's bc's with
    // the integrators.
    if (d_thn_integrator) d_thn_bc_coef = d_thn_integrator->getPhysicalBcCoefs(d_thn_cc_var)[0];

    // Create the hierarchy data operations
    auto hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops =
        hier_ops_manager->getOperationsDouble(new CellVariable<NDIM, double>("cc_var"), hierarchy, true /*get_unique*/);
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), hierarchy, true /*get_unique*/);
    d_hier_fc_data_ops =
        hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM, double>("fc_var"), hierarchy, true /*get_unique*/);

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
MultiphaseStaggeredHierarchyIntegrator::initializeLevelDataSpecialized(Pointer<BasePatchHierarchy<NDIM>> hierarchy,
                                                                       int level_number,
                                                                       double init_data_time,
                                                                       bool can_be_refined,
                                                                       bool initial_time,
                                                                       Pointer<BasePatchLevel<NDIM>> old_level,
                                                                       bool allocate_data)
{
    // Do any kind of Hierarchy specific initialization on a given patch level.
    // Note: All initialization and regridding of state variables is managed by HierarchyIntegrator.
    // Example uses of this include mantaining a vorticity value to see where things should be refined or checking norms
    // of quantities before and after regridding.
    Pointer<PatchLevel<NDIM>> new_level = hierarchy->getPatchLevel(level_number);

    // Compute gradient of theta for initial tagging
    if (initial_time)
    {
        HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_cc_data_ops =
            hier_ops_manager->getOperationsDouble(d_grad_thn_var, d_hierarchy, true);

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int grad_thn_idx = var_db->mapVariableAndContextToIndex(d_grad_thn_var, getCurrentContext());
        const int thn_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getCurrentContext());

        // Fill in with initial conditions.
        d_thn_init_fcn->setDataOnPatchLevel(thn_idx, d_thn_cc_var, new_level, init_data_time, initial_time);

        // Now fill in ghost cells
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp(1);
        ghost_cell_comp[0] = ITC(thn_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, d_thn_bc_coef);

        Pointer<HierarchyGhostCellInterpolation> ghost_cell_fill = new HierarchyGhostCellInterpolation();
        // Note: when we create a new level, the coarser levels should already be created. So we can use that data to
        // fill ghost cells.
        ghost_cell_fill->initializeOperatorState(ghost_cell_comp, hierarchy, 0, level_number);
        HierarchyMathOps hier_math_ops("HierMathOps", hierarchy, 0, level_number);
        hier_math_ops.grad(grad_thn_idx, d_grad_thn_var, 1.0, thn_idx, d_thn_cc_var, ghost_cell_fill, init_data_time);
        // TODO: Replace this with a max over the L2 norm.
        d_max_grad_thn = hier_cc_data_ops->maxNorm(grad_thn_idx, IBTK::invalid_index);

        // Fill in initial values for drawing variables
        if (d_visit_writer)
        {
            const int fn_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, getCurrentContext());
            const int fs_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, getCurrentContext());
            for (PatchLevel<NDIM>::Iterator p(new_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = new_level->getPatch(p());
                Pointer<SideData<NDIM, double>> fn_data = patch->getPatchData(fn_idx);
                fn_data->fillAll(0.0);
                Pointer<SideData<NDIM, double>> fs_data = patch->getPatchData(fs_idx);
                fs_data->fillAll(0.0);
            }
        }
    }
}

void
MultiphaseStaggeredHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM>> hierarchy,
    const int level_num,
    const double error_data_time,
    const int tag_idx,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    // Fill in the tag_idx with 1 if the cell index on a given level should be refined.
    // tag_idx corresponds to patch data CellData<NDIM, int>.
    // Note: d_grad_thn_var is filled during postprocessHierarchyIntegrator.
    if (d_use_grad_tagging)
    {
        // Determine thresholds for the level
        double grad_rel_thresh = std::numeric_limits<double>::max();
        double grad_abs_thresh = std::numeric_limits<double>::max();
        if (d_abs_grad_thresh.size() > 0)
            grad_abs_thresh = d_abs_grad_thresh[std::max(std::min(level_num, d_abs_grad_thresh.size() - 1), 0)];
        if (d_rel_grad_thresh.size() > 0)
            grad_rel_thresh = d_rel_grad_thresh[std::max(std::min(level_num, d_rel_grad_thresh.size() - 1), 0)];
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            Pointer<CellData<NDIM, double>> grad_thn_data = patch->getPatchData(d_grad_thn_var, getCurrentContext());
            Pointer<CellData<NDIM, int>> tagged_data = patch->getPatchData(tag_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++) // cell-centers
            {
                const CellIndex<NDIM>& idx = ci();

                double grad_thn_norm = 0.0;
                for (int d = 0; d < NDIM; ++d) grad_thn_norm += (*grad_thn_data)(idx, d) * (*grad_thn_data)(idx, d);
                grad_thn_norm = std::sqrt(grad_thn_norm);

                if (grad_thn_norm > grad_abs_thresh || (grad_thn_norm / d_max_grad_thn) > grad_rel_thresh)
                    (*tagged_data)(idx) = 1;
            }
        }
    }
    return;
} // applyGradientDetectorSpecialized

void
MultiphaseStaggeredHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
    int coarsest_level,
    int finest_level)
{
    // Do any kind of Hierarchy specific configuration on a new patch hierarchy.
    // If there are synchronization objects that need to be created, do that here.

    // Reset the hierarchy operations objects.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, hierarchy->getFinestLevelNumber());

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, hierarchy->getFinestLevelNumber());

    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, hierarchy->getFinestLevelNumber());
}

void
MultiphaseStaggeredHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
    if (!d_creeping_flow) d_convec_op.reset();
}

void
MultiphaseStaggeredHierarchyIntegrator::initializeCompositeHierarchyDataSpecialized(const double init_data_time,
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
MultiphaseStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                     const double new_time,
                                                                     const int num_cycles)
{
    // Do anything that needs to be done before we call integrateHierarchy().
    INSHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

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
    d_nul_vecs.resize(1 + (d_has_vel_nullspace ? NDIM : 0));
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
            Pointer<CellData<NDIM, double>> p_data = d_nul_vecs[0]->getComponentPatchData(2, *patch);
            p_data->fillAll(1.0);

            if (d_has_vel_nullspace)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    Pointer<SideData<NDIM, double>> un_data = d_nul_vecs[axis + 1]->getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double>> us_data = d_nul_vecs[axis + 1]->getComponentPatchData(1, *patch);
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
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();

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

        // Set Thn boundary conditions
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp(1);
        ghost_cell_comp[0] =
            ITC(thn_cur_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, d_thn_bc_coef);
        HierarchyGhostCellInterpolation hier_ghost_fill;
        hier_ghost_fill.initializeOperatorState(
            ghost_cell_comp, d_hierarchy, 0 /*coarsest_ln*/, d_hierarchy->getFinestLevelNumber());
        hier_ghost_fill.fillData(current_time);

        multiply_sc_and_thn(d_fn_scr_idx, un_cur_idx, thn_cur_idx, d_hierarchy);
        multiply_sc_and_ths(d_fs_scr_idx, us_cur_idx, thn_cur_idx, d_hierarchy);

        d_hier_sc_data_ops->axpy(rhs_un_idx, d_params.rho * (alpha + 1.0) / (dt * alpha), d_fn_scr_idx, rhs_un_idx);
        d_hier_sc_data_ops->axpy(rhs_us_idx, d_params.rho * (alpha + 1.0) / (dt * alpha), d_fs_scr_idx, rhs_us_idx);
    }
    else
    {
        MultiphaseStaggeredStokesOperator RHS_op("RHS_op", true, d_params);
        RHS_op.setPhysicalBoundaryHelper(bc_helper);
        RHS_op.setPhysicalBcCoefs(d_un_bc_coefs, d_us_bc_coefs, nullptr, d_thn_bc_coef);
        // Divergence free condition and pressure are not time stepped. We do not need to account for the contributions
        // in the RHS.
        RHS_op.setCandDCoefficients(C, D1, 0.0, 0.0);
        RHS_op.setThnIdx(thn_cur_idx); // Values at time t_n

        // Store results of applying stokes operator in rhs_vec
        RHS_op.initializeOperatorState(*d_sol_vec, *d_rhs_vec);
        RHS_op.apply(*d_sol_vec, *d_rhs_vec);
    }

    // Set up the operators and solvers needed to solve the linear system.
    d_stokes_op = new MultiphaseStaggeredStokesOperator("stokes_op", false, d_params);
    d_stokes_op->setPhysicalBcCoefs(d_un_bc_coefs, d_us_bc_coefs, nullptr, d_thn_bc_coef);
    d_stokes_op->setCandDCoefficients(C, D2);
    d_stokes_op->setThnIdx(thn_new_idx); // Approximation at time t_{n+1}
    d_stokes_op->setPhysicalBoundaryHelper(bc_helper);

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
                                                               thn_new_idx,
                                                               d_nul_vecs,
                                                               d_un_bc_coefs,
                                                               d_us_bc_coefs,
                                                               d_thn_bc_coef);
        d_stokes_solver->setPreconditioner(d_precond->getPreconditioner());
    }

    d_stokes_solver->setNullspace(false, d_nul_vecs);
    d_stokes_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec); // Error after 1 timestep with BLOCK

    // Set thn_cc_idx on the dense hierarchy.
    if (d_precond)
    {
        if (d_precond->updateVolumeFraction(thn_new_idx, d_thn_cc_var, new_time, d_thn_fcn, d_thn_bc_coef))
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
        const int U_adv_diff_cur_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_integrator->getCurrentContext());
        const int U_adv_diff_scr_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_integrator->getScratchContext());
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, adv_diff_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_cur_idx)) copy_side_to_face(U_adv_diff_cur_idx, un_cur_idx, d_hierarchy);
        adv_diff_integrator->preprocessIntegrateHierarchy(current_time, new_time, adv_diff_num_cycles);
        if (isAllocatedPatchData(U_adv_diff_scr_idx))
            d_hier_fc_data_ops->copyData(U_adv_diff_scr_idx, U_adv_diff_cur_idx);
        if (isAllocatedPatchData(U_adv_diff_new_idx))
            d_hier_fc_data_ops->copyData(U_adv_diff_new_idx, U_adv_diff_cur_idx);
    }

    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
MultiphaseStaggeredHierarchyIntegrator::integrateHierarchySpecialized(const double current_time,
                                                                      const double new_time,
                                                                      const int cycle_num)
{
    INSHierarchyIntegrator::integrateHierarchySpecialized(current_time, new_time, cycle_num);
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

    // Update the Multigrid preconditioner with new volume fraction.
    if (d_precond)
    {
        bool needs_reallocation =
            d_precond->updateVolumeFraction(thn_new_idx, d_thn_cc_var, new_time, nullptr, d_thn_bc_coef);
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

    // Now compute momentum forces
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

    // Reset the solve vector to copy the "scratch" data into the "new" data
    d_hier_sc_data_ops->copyData(un_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_sc_data_ops->copyData(us_new_idx, d_sol_vec->getComponentDescriptorIndex(1));
    d_hier_cc_data_ops->copyData(p_new_idx, d_sol_vec->getComponentDescriptorIndex(2));

    if (d_normalize_pressure)
    {
        const double int_p =
            d_hier_cc_data_ops->integral(p_new_idx, d_hier_math_ops->getCellWeightPatchDescriptorIndex());
        d_hier_cc_data_ops->addScalar(p_new_idx, p_new_idx, -int_p);
    }
    // Reset the RHS
    d_rhs_vec->subtract(d_rhs_vec, f_vec);

    // Destroy the forces
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

        // Now update the state variables. Note that cycle 0 has already been performed
        const int adv_diff_num_cycles = adv_diff_hier_integrator->getNumberOfCycles();
        if (d_current_num_cycles != adv_diff_num_cycles)
        {
            for (int adv_diff_cycle_num = 1; adv_diff_cycle_num < adv_diff_num_cycles; ++adv_diff_cycle_num)
                adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, adv_diff_cycle_num);
        }
    }

    return;
} // integrateHierarchy

void
MultiphaseStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                      const double new_time,
                                                                      const bool skip_synchronize_new_state_data,
                                                                      const int num_cycles)
{
    // Do anything that needs to be done after integrateHierarchy().
    INSHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

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
        const int thn_cur_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getCurrentContext());
        const int fn_cur_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, getCurrentContext());
        const int fs_cur_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, getCurrentContext());

        multiply_sc_and_thn(un_old_idx, un_cur_idx, thn_cur_idx, d_hierarchy);
        multiply_sc_and_ths(us_old_idx, us_cur_idx, thn_cur_idx, d_hierarchy);
        d_hier_sc_data_ops->copyData(fn_old_idx, fn_cur_idx);
        d_hier_sc_data_ops->copyData(fs_old_idx, fs_cur_idx);
    }

    // Calculate gradient of thn for grid cell tagging.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int grad_thn_idx = var_db->mapVariableAndContextToIndex(d_grad_thn_var, getCurrentContext());
    const int thn_new_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getNewContext());
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(1);
    ghost_cell_comp[0] = ITC(thn_new_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, d_thn_bc_coef);
    HierarchyGhostCellInterpolation hier_ghost_fill;
    hier_ghost_fill.initializeOperatorState(
        ghost_cell_comp, d_hierarchy, 0 /*coarsest_ln*/, d_hierarchy->getFinestLevelNumber());
    d_hier_math_ops->grad(grad_thn_idx,
                          d_grad_thn_var,
                          1.0 /*alpha*/,
                          thn_new_idx,
                          d_thn_cc_var,
                          Pointer<HierarchyGhostCellInterpolation>(&hier_ghost_fill, false),
                          new_time);
    // TODO: Replace this with a max over the L2 norm.
    d_max_grad_thn = d_hier_cc_data_ops->maxNorm(grad_thn_idx, IBTK::invalid_index);

    // Synchronize new state
    if (!skip_synchronize_new_state_data) synchronizeHierarchyData(NEW_DATA);

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
    return;
} // postprocessIntegrateHierarchy

void
MultiphaseStaggeredHierarchyIntegrator::regridProjection()
{
    // During regridding, the coarsening and refining operations can introduce errors. Here, we project the velocity
    // onto a divergence free field.
    // TODO: Do we need to implement this? How do we project something that has a co-incompressibility condition.
    // TODO: Check divergence norms to ensure our regridding process does not introduce too large of errors.
}

double
MultiphaseStaggeredHierarchyIntegrator::getStableTimestep(Pointer<Patch<NDIM>> /*patch*/) const
{
    return d_dt_init;
}

void
MultiphaseStaggeredHierarchyIntegrator::setupPlotDataSpecialized()
{
    // Interpolate velocities to node centered.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_draw_idx = var_db->mapVariableAndContextToIndex(d_un_draw_var, getCurrentContext());
    const int us_draw_idx = var_db->mapVariableAndContextToIndex(d_us_draw_var, getCurrentContext());
    const int fn_draw_idx = var_db->mapVariableAndContextToIndex(d_fn_draw_var, getCurrentContext());
    const int fs_draw_idx = var_db->mapVariableAndContextToIndex(d_fs_draw_var, getCurrentContext());
    const int div_draw_idx = var_db->mapVariableAndContextToIndex(d_div_draw_var, getCurrentContext());
    const int un_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getCurrentContext());
    const int us_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getCurrentContext());
    const int fn_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, getCurrentContext());
    const int fs_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, getCurrentContext());
    const int un_scr_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getScratchContext());
    const int us_scr_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getScratchContext());
    const int thn_cur_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getCurrentContext());

    static const bool synch_cf_interface = true;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate scratch data if necessary
    allocatePatchData(un_scr_idx, 0.0, coarsest_ln, finest_ln);
    allocatePatchData(us_scr_idx, 0.0, coarsest_ln, finest_ln);

    // Set thn if necessary
    if (d_thn_fcn)
        d_thn_fcn->setDataOnPatchHierarchy(
            thn_cur_idx, d_thn_cc_var, d_hierarchy, d_integrator_time, false, coarsest_ln, finest_ln);
    if (d_thn_integrator)
    {
        // Copy thn if necessary.
        const int thn_adv_cur_idx =
            var_db->mapVariableAndContextToIndex(d_thn_cc_var, d_thn_integrator->getCurrentContext());
        d_hier_cc_data_ops->copyData(thn_cur_idx, thn_adv_cur_idx);
    }

    // We need ghost cells to interpolate to nodes.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_comp = {
        ITC(un_scr_idx,
            un_idx,
            "CONSERVATIVE_LINEAR_REFINE",
            true,
            "CONSERVATIVE_COARSEN",
            "LINEAR",
            false,
            d_un_bc_coefs),
        ITC(us_scr_idx,
            us_idx,
            "CONSERVATIVE_LINEAR_REFINE",
            true,
            "CONSERVATIVE_COARSEN",
            "LINEAR",
            false,
            d_us_bc_coefs),
        ITC(thn_cur_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, d_thn_bc_coef)
    };
    HierarchyGhostCellInterpolation hier_bdry_fill;
    hier_bdry_fill.initializeOperatorState(ghost_comp, d_hierarchy);
    hier_bdry_fill.fillData(d_integrator_time);

    // interpolate
    d_hier_math_ops->interp(un_draw_idx,
                            d_un_draw_var,
                            synch_cf_interface,
                            un_scr_idx,
                            d_un_sc_var,
                            nullptr,
                            d_integrator_time,
                            synch_cf_interface);
    d_hier_math_ops->interp(us_draw_idx,
                            d_us_draw_var,
                            synch_cf_interface,
                            us_scr_idx,
                            d_us_sc_var,
                            nullptr,
                            d_integrator_time,
                            synch_cf_interface);
    d_hier_math_ops->interp(
        fn_draw_idx, d_fn_draw_var, fn_idx, d_f_un_sc_var, nullptr, d_integrator_time, synch_cf_interface);
    d_hier_math_ops->interp(
        fs_draw_idx, d_fs_draw_var, fs_idx, d_f_us_sc_var, nullptr, d_integrator_time, synch_cf_interface);

    // Compute divergence of velocity fields.
    pre_div_interp(un_scr_idx, thn_cur_idx, un_scr_idx, us_scr_idx, d_hierarchy);
    d_hier_math_ops->div(div_draw_idx, d_div_draw_var, 1.0, un_scr_idx, d_un_sc_var, nullptr, 0.0, true);
    ghost_comp = { ITC(div_draw_idx, "NONE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr) };
    hier_bdry_fill.deallocateOperatorState();
    hier_bdry_fill.initializeOperatorState(ghost_comp, d_hierarchy);
    hier_bdry_fill.fillData(d_integrator_time);

    // Deallocate scratch data
    deallocatePatchData(un_scr_idx, coarsest_ln, finest_ln);
    deallocatePatchData(us_scr_idx, coarsest_ln, finest_ln);
}

void
MultiphaseStaggeredHierarchyIntegrator::setThnAtHalf(int& thn_cur_idx,
                                                     int& thn_new_idx,
                                                     int& thn_scr_idx,
                                                     const double current_time,
                                                     const double new_time,
                                                     const bool start_of_ts)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    thn_cur_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getCurrentContext());
    thn_new_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getNewContext());
    thn_scr_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getScratchContext());

    if (d_thn_integrator)
    {
        // We are evolving thn, so use the integrator to copy data.
        const int thn_evolve_cur_idx =
            var_db->mapVariableAndContextToIndex(d_thn_cc_var, d_thn_integrator->getCurrentContext());
        const int thn_evolve_new_idx =
            var_db->mapVariableAndContextToIndex(d_thn_cc_var, d_thn_integrator->getNewContext());
        d_hier_cc_data_ops->copyData(thn_cur_idx, thn_evolve_cur_idx);
        // We set the "new" data to the best approximation we have, which if we are at the beginning of the time step,
        // is the current data.
        d_hier_cc_data_ops->copyData(thn_new_idx, start_of_ts ? thn_evolve_cur_idx : thn_evolve_new_idx);
        d_hier_cc_data_ops->linearSum(thn_scr_idx, 0.5, thn_cur_idx, 0.5, thn_new_idx);
    }
    else
    {
        // Otherwise set the values with the function
#ifndef NDEBUG
        TBOX_ASSERT(d_thn_fcn);
#endif
        if (start_of_ts)
        {
            d_thn_fcn->setDataOnPatchHierarchy(thn_cur_idx, d_thn_cc_var, d_hierarchy, current_time, false);
            d_thn_fcn->setDataOnPatchHierarchy(thn_new_idx, d_thn_cc_var, d_hierarchy, new_time, false);
        }
        double half_time = 0.5 * (current_time + new_time);
        d_thn_fcn->setDataOnPatchHierarchy(thn_scr_idx, d_thn_cc_var, d_hierarchy, half_time, false);
    }
    return;
}

void
MultiphaseStaggeredHierarchyIntegrator::approxConvecOp(Pointer<SAMRAIVectorReal<NDIM, double>>& f_vec,
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

bool
MultiphaseStaggeredHierarchyIntegrator::isVariableDrag() const
{
    return d_params.isVariableDrag();
}

void
MultiphaseStaggeredHierarchyIntegrator::addBodyForces(Pointer<SAMRAIVectorReal<NDIM, double>>& f_vec,
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

    // Now account for scaled forces. We need to compute a volume fraction with ghost cells.
    {
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp = { ITC(
            thn_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", true, nullptr) };
        HierarchyGhostCellInterpolation ghost_fill;
        ghost_fill.initializeOperatorState(ghost_cell_comp, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
        ghost_fill.fillData(eval_time);
    }
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

    deallocate_patch_data({ fn_cloned_idx, fs_cloned_idx }, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    var_db->removePatchDataIndex(fn_cloned_idx);
    var_db->removePatchDataIndex(fs_cloned_idx);
}

MultiphaseStaggeredHierarchyIntegrator::MultiphasePreconditioner::MultiphasePreconditioner(
    PreconditionerType precond_type,
    Pointer<PatchHierarchy<NDIM>> hierarchy,
    Pointer<Database> input_db,
    const MultiphaseParameters& params,
    const double C,
    const double D,
    const int thn_idx,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>& null_vecs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
    RobinBcCoefStrategy<NDIM>* thn_bc_coef)
    : d_hierarchy(hierarchy), d_precond_type(precond_type)
{
    switch (d_precond_type)
    {
    case PreconditionerType::BLOCK:
    {
        d_block_precond = new MultiphaseStaggeredStokesBlockPreconditioner("KrylovPrecondStrategy", params, input_db);
        d_block_precond->setThnIdx(thn_idx);
        d_block_precond->setCAndDCoefficients(C, D);
        d_block_precond->setNullspace(false, null_vecs);
        break;
    }
    case PreconditionerType::MULTIGRID:
    {
        d_fac_op =
            new MultiphaseStaggeredStokesBoxRelaxationFACOperator("KrylovPrecondStrategy", "Krylov_precond_", params);
        d_fac_op->setThnIdx(thn_idx); // Approximation at time t_{n+1}
        d_fac_op->setUnderRelaxationParamater(input_db->getDouble("relaxation_parameter"));
        d_fac_op->setCandDCoefficients(C, D);
        d_fac_op->setPhysicalBcCoefs(un_bc_coefs, us_bc_coefs, nullptr, thn_bc_coef);
        d_fac_precond = new FullFACPreconditioner("KrylovPrecond", d_fac_op, input_db, "Krylov_precond_");
        d_fac_precond->setNullspace(false, null_vecs);
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
MultiphaseStaggeredHierarchyIntegrator::MultiphasePreconditioner::updateVolumeFraction(
    const int thn_idx,
    Pointer<CellVariable<NDIM, double>>& thn_var,
    const double time,
    IBTK::CartGridFunction* thn_fcn,
    RobinBcCoefStrategy<NDIM>* thn_bc_coef)
{
    switch (d_precond_type)
    {
    case PreconditionerType::BLOCK:
    {
        d_block_precond->deallocateSolverState();
        d_block_precond->setThnIdx(thn_idx);
        return true;
        break;
    }
    case PreconditionerType::MULTIGRID:
    {
        Pointer<PatchHierarchy<NDIM>> dense_hierarchy = d_fac_precond->getDenseHierarchy();
        if (thn_fcn)
        {
            // Allocate data
            for (int ln = 0; ln <= dense_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = dense_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(thn_idx)) level->allocatePatchData(thn_idx, time);
            }
            thn_fcn->setDataOnPatchHierarchy(
                thn_idx, thn_var, dense_hierarchy, time, false, 0, dense_hierarchy->getFinestLevelNumber());
        }
        else
        {
            d_fac_precond->transferToDense(thn_idx, true);
        }

        // Also fill in ghost cells on the dense hierarchy
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp(1);
        ghost_cell_comp[0] = ITC(thn_idx,
                                 "CONSERVATIVE_LINEAR_REFINE",
                                 false,
                                 "NONE",
                                 "LINEAR",
                                 true,
                                 thn_bc_coef); // defaults to fill corner
        HierarchyGhostCellInterpolation ghost_cell_fill;
        ghost_cell_fill.initializeOperatorState(
            ghost_cell_comp, dense_hierarchy, 0, dense_hierarchy->getFinestLevelNumber());
        ghost_cell_fill.fillData(time);
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
MultiphaseStaggeredHierarchyIntegrator::MultiphasePreconditioner::updateDragCoefficient(
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
