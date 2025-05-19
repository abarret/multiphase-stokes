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
MultiphaseStaggeredHierarchyIntegrator::MultiphaseStaggeredHierarchyIntegrator(std::string object_name,
                                                                               Pointer<Database> input_db,
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
                             register_for_restart),
      d_un_sc_var(new SideVariable<NDIM, double>(d_object_name + "::un_sc")),
      d_us_sc_var(new SideVariable<NDIM, double>(d_object_name + "::us_sc")),
      d_f_un_sc_var(new SideVariable<NDIM, double>(d_object_name + "::f_un_sc")),
      d_f_us_sc_var(new SideVariable<NDIM, double>(d_object_name + "::f_us_sc")),
      d_f_cc_var(new CellVariable<NDIM, double>(d_object_name + "::f_cc")),
      d_thn_cc_var(new CellVariable<NDIM, double>(d_object_name + "::thn_cc"))
{
    if (input_db->keyExists("rho")) d_params.rho = input_db->getDouble("rho");
    if (input_db->keyExists("eta_n")) d_params.eta_n = input_db->getDouble("eta_n");
    if (input_db->keyExists("eta_s")) d_params.eta_s = input_db->getDouble("eta_s");
    if (input_db->keyExists("l_n"))
        d_params.lambda_n = input_db->getDouble("l_n");
    else
        d_params.lambda_n = d_params.eta_n;
    if (input_db->keyExists("l_s"))
        d_params.lambda_s = input_db->getDouble("l_s");
    else
        d_params.lambda_s = d_params.eta_s;
    if (input_db->keyExists("use_grad_tagging")) d_use_grad_tagging = input_db->getBool("use_grad_tagging");
    if (input_db->keyExists("grad_rel_thresh")) input_db->getArray("grad_rel_thresh", d_rel_grad_thresh);
    if (input_db->keyExists("grad_abs_thresh")) input_db->getArray("grad_abs_thresh", d_abs_grad_thresh);
    d_use_accel_ts = input_db->getBoolWithDefault("use_accel_ts", d_use_accel_ts);
    d_accel_ts_safety_fac = input_db->getDoubleWithDefault("accel_ts_safety_factor", d_accel_ts_safety_fac);
    d_regularize_thn = input_db->getDoubleWithDefault("regualarize_thn", d_regularize_thn);

    // TODO: The default here should really be "false", but for now, this will not change the default behavior.
    d_creeping_flow = input_db->getBoolWithDefault("creeping_flow", true);

    // Arbitrarily set the base class velocity variable to the solvent variable. Note this is required for the base
    // class to compute the CFL number.
    d_U_var = d_us_sc_var;

    d_un_bc_coefs.resize(NDIM, nullptr);
    d_us_bc_coefs.resize(NDIM, nullptr);

    d_U_var = d_us_sc_var;

    // Create some dummy boundary condition objects. These can be overwritten by calling
    // registerPhysicalBoundaryConditions().
    d_un_bc_coefs.resize(NDIM);
    d_us_bc_coefs.resize(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        d_dummy_bcs.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("u_dummy", nullptr));
        d_un_bc_coefs[d] = d_dummy_bcs[d].get();
        d_us_bc_coefs[d] = d_dummy_bcs[d].get();
    }
    d_U_bc_coefs = d_us_bc_coefs;

    // Advection velocity variable for the solvent.
    d_us_adv_diff_var = new FaceVariable<NDIM, double>(d_object_name + "::us_adv_var");

    IBTK_DO_ONCE(t_integrate_hierarchy = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStaggeredHierarchyIntegrator::integrateHierarchy()");
                 t_preprocess_integrate_hierarchy = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStaggeredHierarchyIntegrator::preprocess()");
                 t_postprocess_integrate_hierarchy = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStaggeredHierarchyIntegrator::postprocess()"););
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
    std::vector<RobinBcCoefStrategy<NDIM>*> us_bc_coefs,
    RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
    for (int d = 0; d < NDIM; ++d)
    {
        if (un_bc_coefs[d]) d_un_bc_coefs[d] = un_bc_coefs[d];
        if (us_bc_coefs[d]) d_us_bc_coefs[d] = us_bc_coefs[d];
    }
    // Treat solvent as U
    d_U_bc_coefs = d_us_bc_coefs;
    d_p_bc_coef = p_bc_coef;
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
MultiphaseStaggeredHierarchyIntegrator::registerAdvDiffHierarchyIntegrator(
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(adv_diff_hier_integrator);
    // Bad things happen if the same integrator is registered twice.
    auto pointer_compare = [adv_diff_hier_integrator](Pointer<AdvDiffHierarchyIntegrator> integrator) -> bool
    { return adv_diff_hier_integrator.getPointer() == integrator.getPointer(); };
    TBOX_ASSERT(std::find_if(d_adv_diff_hier_integrators.begin(), d_adv_diff_hier_integrators.end(), pointer_compare) ==
                d_adv_diff_hier_integrators.end());
#endif
    d_adv_diff_hier_integrators.push_back(adv_diff_hier_integrator);
    registerChildHierarchyIntegrator(adv_diff_hier_integrator);
    adv_diff_hier_integrator->registerAdvectionVelocity(d_U_adv_diff_var);
    adv_diff_hier_integrator->setAdvectionVelocityIsDivergenceFree(d_U_adv_diff_var, false);
    adv_diff_hier_integrator->registerAdvectionVelocity(d_us_adv_diff_var);
    adv_diff_hier_integrator->setAdvectionVelocityIsDivergenceFree(d_us_adv_diff_var, false);
    return;
} // registerAdvDiffHierarchyIntegrator

void
MultiphaseStaggeredHierarchyIntegrator::setViscosityCoefficient(const double eta_n, const double eta_s)
{
    setViscosityCoefficient(eta_n, eta_s, eta_n, eta_s);
}

void
MultiphaseStaggeredHierarchyIntegrator::setViscosityCoefficient(const double eta_n,
                                                                const double eta_s,
                                                                const double lambda_n,
                                                                const double lambda_s)
{
    d_params.eta_n = eta_n;
    d_params.eta_s = eta_s;
    d_params.lambda_n = lambda_n;
    d_params.lambda_s = lambda_s;
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
MultiphaseStaggeredHierarchyIntegrator::setForcingFunctionsScaledByBoth(Pointer<CartGridFunction> fn_fcn,
                                                                        Pointer<CartGridFunction> fs_fcn)
{
    if (fn_fcn)
    {
        if (d_f_un_thn_ths_fcn)
        {
            Pointer<CartGridFunctionSet> p_f_un_fcn = d_f_un_thn_ths_fcn;
            if (p_f_un_fcn)
            {
                p_f_un_fcn->addFunction(fn_fcn);
            }
            else
            {
                p_f_un_fcn = new CartGridFunctionSet(d_object_name + "::fn_thn_ths_fcn_set");
                p_f_un_fcn->addFunction(d_f_un_thn_ths_fcn);
                p_f_un_fcn->addFunction(fn_fcn);
                d_f_un_thn_ths_fcn = p_f_un_fcn;
            }
        }
        else
        {
            d_f_un_thn_ths_fcn = fn_fcn;
        }
    }

    if (fs_fcn)
    {
        if (d_f_us_thn_ths_fcn)
        {
            Pointer<CartGridFunctionSet> p_f_us_fcn = d_f_us_thn_ths_fcn;
            if (p_f_us_fcn)
            {
                p_f_us_fcn->addFunction(fs_fcn);
            }
            else
            {
                p_f_us_fcn = new CartGridFunctionSet(d_object_name + "::fs_thn_ths_fcn_set");
                p_f_us_fcn->addFunction(d_f_us_thn_ths_fcn);
                p_f_us_fcn->addFunction(fs_fcn);
                d_f_us_thn_ths_fcn = p_f_us_fcn;
            }
        }
        else
        {
            d_f_us_thn_ths_fcn = fs_fcn;
        }
    }
}

void
MultiphaseStaggeredHierarchyIntegrator::setInitialNetworkVolumeFraction(Pointer<CartGridFunction> thn_init_fcn)
{
    d_thn_init_fcn = thn_init_fcn;
    // If we have already registered the volume fraction to be advected, overwrite the initial conditions
    if (d_thn_integrator) d_thn_integrator->setInitialConditions(d_thn_cc_var, d_thn_init_fcn);
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
    RobinBcCoefStrategy<NDIM>* thn_bc_coef,
    const bool has_meaningful_mid_value)
{
    registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
    d_thn_integrator = adv_diff_integrator;

    // Set up thn to be advected
    d_thn_integrator->registerTransportedQuantity(d_thn_cc_var, true /*output_Q*/);
    d_thn_integrator->setAdvectionVelocity(d_thn_cc_var, d_U_adv_diff_var);
    d_thn_integrator->setAdvectionVelocityIsDivergenceFree(d_U_adv_diff_var, false); // Not divergence free in general.
    d_thn_integrator->setInitialConditions(d_thn_cc_var, d_thn_init_fcn);
    d_thn_integrator->setPhysicalBcCoef(d_thn_cc_var, thn_bc_coef);
    d_thn_bc_coef = thn_bc_coef;

    d_use_new_thn = has_meaningful_mid_value;
}

void
MultiphaseStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                                      Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Set the correct boundary conditions. If we advect the volume fraction, then we overwrite this class's bc's with
    // the integrators.
    if (d_thn_integrator) d_thn_bc_coef = d_thn_integrator->getPhysicalBcCoefs(d_thn_cc_var)[0];

    // Here we do all we need to ensure that calls to advanceHierarchy() or integrateHierarchy() are valid.
    // NOTE: This function is called before the patch hierarchy has valid patch levels.
    // To set initial data, we should do this in initializeLevelDataSpecialized().

    // First create the variables we need.
    // NOTE: d_P_var is a member variable of the base class.
    d_grad_thn_var = new CellVariable<NDIM, double>(d_object_name + "grad_thn_cc", NDIM);
    d_un_rhs_var = new SideVariable<NDIM, double>(d_object_name + "::un_rhs");
    d_us_rhs_var = new SideVariable<NDIM, double>(d_object_name + "::us_rhs");
    d_p_rhs_var = new CellVariable<NDIM, double>(d_object_name + "::p_rhs");

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

    // Setup the volume fraction managers
    d_thn_cur_manager = std::make_unique<VolumeFractionDataManager>(
        d_object_name + "::ThnCurManager", d_thn_cc_var, getCurrentContext(), d_thn_bc_coef, d_regularize_thn);
    d_thn_scr_manager = std::make_unique<VolumeFractionDataManager>(
        d_object_name + "::ThnScrManager", d_thn_cc_var, getScratchContext(), d_thn_bc_coef, d_regularize_thn);
    d_thn_new_manager = std::make_unique<VolumeFractionDataManager>(
        d_object_name + "::ThnNewManager", d_thn_cc_var, getNewContext(), d_thn_bc_coef, d_regularize_thn);

    // Make sure the current data is maintained after regridding
    const int thn_cur_idx = d_thn_cur_manager->getCellIndex();
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineOperator<NDIM>> thn_refine_op =
        grid_geom->lookupRefineOperator(d_thn_cur_manager->getCellVariable(), "CONSERVATIVE_LINEAR_REFINE");
    d_fill_after_regrid_bc_idxs.setFlag(thn_cur_idx);
    d_fill_after_regrid_prolong_alg.registerRefine(thn_cur_idx, thn_cur_idx, thn_cur_idx, thn_refine_op);

    // Set the scratch and new indices according to be allocated with the corresponding context
    d_scratch_data.setFlag(d_thn_scr_manager->getCellIndex());
    d_new_data.setFlag(d_thn_new_manager->getCellIndex());

    // Everything else only gets a scratch context, which is deallocated at the end of each time step.
    // Note the forces need ghost cells for modifying the RHS to account for non-homogenous boundary conditions.
    int f_p_idx, f_un_idx, f_us_idx, un_rhs_idx, us_rhs_idx, p_rhs_idx, xi_idx;
    registerVariable(f_p_idx, d_f_cc_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(f_un_idx, d_f_un_sc_var, IntVector<NDIM>(1), getCurrentContext());
    registerVariable(f_us_idx, d_f_us_sc_var, IntVector<NDIM>(1), getCurrentContext());
    registerVariable(un_rhs_idx, d_un_rhs_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(us_rhs_idx, d_us_rhs_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(p_rhs_idx, d_p_rhs_var, IntVector<NDIM>(1), getScratchContext());
    if (d_xi_var)
    {
        registerVariable(xi_idx, d_xi_var, IntVector<NDIM>(0), getScratchContext());
        d_params.xi_idx = xi_idx;
        pout << "Setting xi idx in parameter list\n";
    }

    // Register a scratch force object.
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
        d_div_un_var = new CellVariable<NDIM, double>(d_object_name + "::Div_Un");
        d_div_us_var = new CellVariable<NDIM, double>(d_object_name + "::Div_Us");
        int un_draw_idx, us_draw_idx, div_draw_idx, fn_draw_idx, fs_draw_idx, div_un_idx, div_us_idx;
        registerVariable(un_draw_idx, d_un_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(us_draw_idx, d_us_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(div_draw_idx, d_div_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(fn_draw_idx, d_fn_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(fs_draw_idx, d_fs_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(div_un_idx, d_div_un_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(div_us_idx, d_div_us_var, IntVector<NDIM>(0), getCurrentContext());

        d_visit_writer->registerPlotQuantity("Un", "VECTOR", un_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Un_" + std::to_string(d), "SCALAR", un_draw_idx, d, 1.0, "NODE");

        d_visit_writer->registerPlotQuantity("Us", "VECTOR", us_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Us_" + std::to_string(d), "SCALAR", us_draw_idx, d, 1.0, "NODE");

        d_visit_writer->registerPlotQuantity("Fn", "VECTOR", fn_draw_idx, 0, 1.0, "CELL");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Fn_" + std::to_string(d), "SCALAR", fn_draw_idx, d, 1.0, "CELL");

        d_visit_writer->registerPlotQuantity("Fs", "VECTOR", fs_draw_idx, 0, 1.0, "CELL");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Fs_" + std::to_string(d), "SCALAR", fs_draw_idx, d, 1.0, "CELL");

        d_visit_writer->registerPlotQuantity("P", "SCALAR", p_cur_idx, 0, 1.0, "CELL");

        d_visit_writer->registerPlotQuantity("Div", "SCALAR", div_draw_idx, 0, 1.0, "CELL");
        d_visit_writer->registerPlotQuantity("Div_Un", "SCALAR", div_un_idx, 0, 1.0, "CELL");
        d_visit_writer->registerPlotQuantity("Div_Us", "SCALAR", div_us_idx, 0, 1.0, "CELL");

        // Only need to plot this variable if we aren't advecting it.
        // If we do advect theta, the advection integrator will plot it.
        if (d_thn_fcn)
            d_visit_writer->registerPlotQuantity("Thn", "SCALAR", d_thn_cur_manager->getCellIndex(), 0, 1.0, "CELL");
    }

    // Create the hierarchy data operations
    auto hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops =
        hier_ops_manager->getOperationsDouble(new CellVariable<NDIM, double>("cc_var"), hierarchy, true /*get_unique*/);
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), hierarchy, true /*get_unique*/);
    d_hier_fc_data_ops =
        hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM, double>("fc_var"), hierarchy, true /*get_unique*/);

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
        const int thn_cur_idx = d_thn_cur_manager->getCellIndex();
        new_level->allocatePatchData(thn_cur_idx, init_data_time);
        Pointer<CellVariable<NDIM, double>> thn_cur_var = d_thn_cur_manager->getCellVariable();

        // Fill in with initial conditions.
        // NOTE: We cannot use the VolumeFractionDataManager here because it can't set values on a given level.
        d_thn_init_fcn->setDataOnPatchLevel(thn_cur_idx, thn_cur_var, new_level, init_data_time, initial_time);

        // Now fill in ghost cells
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp(1);
        ghost_cell_comp[0] =
            ITC(thn_cur_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, d_thn_bc_coef);

        Pointer<HierarchyGhostCellInterpolation> ghost_cell_fill = new HierarchyGhostCellInterpolation();
        // Note: when we create a new level, the coarser levels should already be created. So we can use that data to
        // fill ghost cells.
        ghost_cell_fill->initializeOperatorState(ghost_cell_comp, hierarchy, 0, level_number);
        HierarchyMathOps hier_math_ops("HierMathOps", hierarchy, 0, level_number);
        hier_math_ops.grad(
            grad_thn_idx, d_grad_thn_var, 1.0, thn_cur_idx, thn_cur_var, ghost_cell_fill, init_data_time);
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
MultiphaseStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                      const double new_time,
                                                                      const bool skip_synchronize_new_state_data,
                                                                      const int num_cycles)
{
    IBTK_TIMER_START(t_postprocess_integrate_hierarchy);
    // Compute the network CFL number. Note the solvent CFL number is computed by the base class and is available in
    // d_cfl_current.
    PatchSideDataOpsReal<NDIM, double> patch_sc_ops;
    double cfl_max = 0.0;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double dx_min = *(std::min_element(dx, dx + NDIM));
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(d_un_sc_var, getNewContext());
#ifndef NDEBUG
            TBOX_ASSERT(un_data);
#endif
            double u_max = patch_sc_ops.maxNorm(un_data, patch_box);
            const double dt = new_time - current_time;
            cfl_max = std::max(cfl_max, u_max * dt / dx_min);
        }
    }
    d_cfl_un_current = IBTK_MPI::maxReduction(cfl_max);
    plog << d_object_name << "::postprocessIntegrateHierarchy(): Network CFL number = " << d_cfl_un_current << "\n";

    // Do anything that needs to be done after integrateHierarchy().
    INSHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Calculate gradient of thn for grid cell tagging.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int grad_thn_idx = var_db->mapVariableAndContextToIndex(d_grad_thn_var, getCurrentContext());
    const int thn_new_idx = d_thn_new_manager->getCellIndex();
    Pointer<CellVariable<NDIM, double>> thn_new_var = d_thn_new_manager->getCellVariable();
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
                          thn_new_var,
                          Pointer<HierarchyGhostCellInterpolation>(&hier_ghost_fill, false),
                          new_time);
    // TODO: Replace this with a max over the L2 norm.
    d_max_grad_thn = d_hier_cc_data_ops->maxNorm(grad_thn_idx, IBTK::invalid_index);

    IBTK_TIMER_STOP(t_postprocess_integrate_hierarchy);
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
    const int div_un_idx = var_db->mapVariableAndContextToIndex(d_div_un_var, getCurrentContext());
    const int div_us_idx = var_db->mapVariableAndContextToIndex(d_div_us_var, getCurrentContext());
    const int un_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getCurrentContext());
    const int us_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getCurrentContext());
    const int fn_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, getCurrentContext());
    const int fs_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, getCurrentContext());
    const int un_scr_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getScratchContext());
    const int us_scr_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getScratchContext());
    const int thn_cur_idx = d_thn_cur_manager->getCellIndex();

    static const bool synch_cf_interface = true;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate scratch data if necessary
    allocatePatchData(un_scr_idx, 0.0, coarsest_ln, finest_ln);
    allocatePatchData(us_scr_idx, 0.0, coarsest_ln, finest_ln);

    // Set thn if necessary
    if (d_thn_fcn)
    {
        d_thn_cur_manager->updateVolumeFraction(
            *d_thn_fcn, *d_hierarchy, d_integrator_time, TimePoint::CURRENT_TIME, false);
    }
    else
    {
        const int thn_adv_cur_idx =
            var_db->mapVariableAndContextToIndex(d_thn_cc_var, d_thn_integrator->getCurrentContext());
        d_thn_cur_manager->updateVolumeFraction(thn_adv_cur_idx, *d_hierarchy, d_integrator_time);
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
    d_hier_math_ops->div(div_un_idx, d_div_un_var, 1.0, un_scr_idx, d_un_sc_var, nullptr, 0.0, true);
    d_hier_math_ops->div(div_us_idx, d_div_us_var, 1.0, us_scr_idx, d_us_sc_var, nullptr, 0.0, true);
    pre_div_interp(un_scr_idx, thn_cur_idx, un_scr_idx, us_scr_idx, d_hierarchy);
    d_hier_math_ops->div(div_draw_idx, d_div_draw_var, 1.0, un_scr_idx, d_un_sc_var, nullptr, 0.0, true);
    ghost_comp = { ITC(div_draw_idx, "NONE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr),
                   ITC(div_un_idx, "NONE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr),
                   ITC(div_us_idx, "NONE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr) };
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
    thn_cur_idx = d_thn_cur_manager->getCellIndex();
    thn_new_idx = d_thn_new_manager->getCellIndex();
    thn_scr_idx = d_thn_scr_manager->getCellIndex();

    const double half_time = 0.5 * (current_time + new_time);

    if (d_thn_integrator)
    {
        // We are evolving thn, so use the integrator to copy data.
        const int thn_evolve_cur_idx =
            var_db->mapVariableAndContextToIndex(d_thn_cc_var, d_thn_integrator->getCurrentContext());
        const int thn_evolve_new_idx =
            var_db->mapVariableAndContextToIndex(d_thn_cc_var, d_thn_integrator->getNewContext());
        d_thn_cur_manager->updateVolumeFraction(thn_evolve_cur_idx, *d_hierarchy, current_time);
        d_thn_new_manager->updateVolumeFraction(
            start_of_ts ? thn_evolve_cur_idx : thn_evolve_new_idx, *d_hierarchy, new_time);
        d_hier_cc_data_ops->linearSum(thn_scr_idx, 0.5, thn_cur_idx, 0.5, thn_new_idx);
        d_thn_scr_manager->updateVolumeFraction(thn_scr_idx, *d_hierarchy, half_time);
    }
    else
    {
        // Otherwise set the values with the function
        d_thn_cur_manager->updateVolumeFraction(*d_thn_fcn, *d_hierarchy, current_time, TimePoint::CURRENT_TIME);
        d_thn_scr_manager->updateVolumeFraction(*d_thn_fcn, *d_hierarchy, half_time, TimePoint::HALF_TIME);
        d_thn_new_manager->updateVolumeFraction(*d_thn_fcn, *d_hierarchy, new_time, TimePoint::NEW_TIME);
    }
    return;
}

bool
MultiphaseStaggeredHierarchyIntegrator::isVariableDrag() const
{
    return d_params.isVariableDrag();
}

double
MultiphaseStaggeredHierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    double dt = INSHierarchyIntegrator::getMaximumTimeStepSizeSpecialized();
    if (!d_use_accel_ts) return dt;
    double dt_loc = std::numeric_limits<double>::max();
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> fn_data = patch->getPatchData(d_f_un_sc_var, getCurrentContext());
            Pointer<SideData<NDIM, double>> fs_data = patch->getPatchData(d_f_us_sc_var, getCurrentContext());

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double dx_min = *std::min_element(dx, dx + NDIM);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                VectorNd fn, fs;
                for (int d = 0; d < NDIM; ++d)
                {
                    fn[d] = 0.5 * ((*fn_data)(SideIndex<NDIM>(idx, d, 0)) + (*fn_data)(SideIndex<NDIM>(idx, d, 1)));
                    fs[d] = 0.5 * ((*fs_data)(SideIndex<NDIM>(idx, d, 0)) + (*fs_data)(SideIndex<NDIM>(idx, d, 1)));
                }
                double f_max = std::max(fn.norm(), fs.norm());
                dt_loc = std::min(dt_loc, d_accel_ts_safety_fac * std::sqrt(dx_min / f_max));
            }
        }
    }
    return std::max(dt, dt_loc);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
