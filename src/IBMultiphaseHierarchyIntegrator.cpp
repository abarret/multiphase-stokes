/////////////////////////////// INCLUDES /////////////////////////////////////

#include "multiphase/IBMultiphaseHierarchyIntegrator.h"
#include "multiphase/MultiphaseStaggeredHierarchyIntegrator.h"

#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/app_namespaces.h" // IWYU pragma: keep
#include "ibamr/ibamr_enums.h"
#include <ibamr/IBMethod.h>

#include "ibtk/CartGridFunction.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_enums.h"
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>

#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsManager.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <ostream>
#include <string>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace multiphase
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBMultiphaseHierarchyIntegrator restart file data.
static const int IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION = 2;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMultiphaseHierarchyIntegrator::IBMultiphaseHierarchyIntegrator(
    std::string object_name,
    Pointer<Database> input_db,
    Pointer<IBStrategy> ib_method_ops,
    Pointer<IBStrategy> ibn_method_ops,
    Pointer<MultiphaseStaggeredHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : IBHierarchyIntegrator(std::move(object_name), input_db, ib_method_ops, ins_hier_integrator, register_for_restart),
      d_ibn_method_ops(ibn_method_ops)
{
    // Set default configuration options.
    d_use_structure_predictor = true;

    // Set options from input.
    if (input_db)
    {
        if (input_db->keyExists("use_structure_predictor"))
            d_use_structure_predictor = input_db->getBool("use_structure_predictor");
    }

    d_ibn_method_ops->registerIBHierarchyIntegrator(this);
    d_u_var = ins_hier_integrator->getSolventVariable();
    d_p_var = ins_hier_integrator->getPressureVariable();
    d_f_var = ins_hier_integrator->getSolventBodyForceVariable();
    d_un_var = ins_hier_integrator->getNetworkVariable();
    d_fn_var = ins_hier_integrator->getNetworkBodyForceVariable();
    d_cross_links_un_var = new SideVariable<NDIM, double>(d_object_name + "::cross_links_un");
    d_cross_links_us_var = new SideVariable<NDIM, double>(d_object_name + "::cross_links_us");

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBMultiphaseHierarchyIntegrator

void
IBMultiphaseHierarchyIntegrator::registerCrossLinkStrategy(Pointer<MultiphaseCrossLinksStrategy> cross_link_strategy)
{
    TBOX_ASSERT(cross_link_strategy);
    d_cross_links_strategy = cross_link_strategy;
}

void
IBMultiphaseHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                              const double new_time,
                                                              const int num_cycles)
{
    // preprocess our dependencies...
    d_ibn_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
    case BACKWARD_EULER:
    case TRAPEZOIDAL_RULE:
        if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(current_time);
        d_ibn_method_ops->computeLagrangianForce(current_time);
        if (d_cross_links_strategy)
            d_cross_links_strategy->computeLagrangianForce(current_time, IBTK::TimePoint::CURRENT_TIME);
        if (d_enable_logging)
            plog << d_object_name
                 << "::preprocessIntegrateHierarchy(): spreading Lagrangian force "
                    "to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), current_time);

        d_hier_velocity_data_ops->setToScalar(d_fn_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_fn_idx);
        d_ib_method_ops->spreadForce(
            d_fn_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::fn"), current_time);

        d_hier_velocity_data_ops->setToScalar(d_cross_links_un_idx, 0.0);
        d_hier_velocity_data_ops->setToScalar(d_cross_links_us_idx, 0.0);

        if (d_cross_links_strategy)
        {
            d_u_phys_bdry_op->setPatchDataIndex(d_cross_links_un_idx);
            d_cross_links_strategy->spreadForce(d_cross_links_un_idx,
                                                d_u_phys_bdry_op,
                                                getProlongRefineSchedules(d_object_name + "::fn"),
                                                current_time,
                                                true /*spread_network*/,
                                                IBTK::TimePoint::CURRENT_TIME);
            d_u_phys_bdry_op->setPatchDataIndex(d_cross_links_us_idx);
            d_cross_links_strategy->spreadForce(d_cross_links_us_idx,
                                                d_u_phys_bdry_op,
                                                getProlongRefineSchedules(d_object_name + "::f"),
                                                current_time,
                                                false /*spread_network*/,
                                                IBTK::TimePoint::CURRENT_TIME);
        }

        d_u_phys_bdry_op->setHomogeneousBc(false);
        if (d_f_current_idx != invalid_index) d_hier_velocity_data_ops->copyData(d_f_current_idx, d_f_idx);
        if (d_fn_current_idx != invalid_index) d_hier_velocity_data_ops->copyData(d_fn_current_idx, d_fn_idx);
        if (d_cross_links_current_un_idx != invalid_index)
            d_hier_velocity_data_ops->copyData(d_cross_links_current_un_idx, d_cross_links_un_idx);
        if (d_cross_links_current_us_idx != invalid_index)
            d_hier_velocity_data_ops->copyData(d_cross_links_current_us_idx, d_cross_links_us_idx);
        break;
    case MIDPOINT_RULE:
        // intentionally blank
        break;
    default:
        TBOX_ERROR(
            d_object_name
            << "::preprocessIntegrateHierarchy():\n"
            << "  unsupported time stepping type: "
            << IBAMR::enum_to_string<IBAMR::TimeSteppingType>(d_time_stepping_type) << "\n"
            << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an initial prediction of the updated positions of the Lagrangian
    // structure.
    //
    // NOTE: The velocity should already have been interpolated to the
    // curvilinear mesh and should not need to be re-interpolated.
    if (d_use_structure_predictor)
    {
        if (d_enable_logging)
            plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian forward Euler step\n";
        d_ib_method_ops->forwardEulerStep(current_time, new_time);
        d_ibn_method_ops->forwardEulerStep(current_time, new_time);
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
IBMultiphaseHierarchyIntegrator::integrateHierarchySpecialized(const double current_time,
                                                               const double new_time,
                                                               const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchySpecialized(current_time, new_time, cycle_num);
    const double half_time = current_time + 0.5 * (new_time - current_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<MultiphaseStaggeredHierarchyIntegrator> ins_hier_integrator = d_ins_hier_integrator;
    const int un_cur_idx = var_db->mapVariableAndContextToIndex(d_un_var, ins_hier_integrator->getCurrentContext());
    const int un_new_idx = var_db->mapVariableAndContextToIndex(d_un_var, ins_hier_integrator->getNewContext());
    const int us_cur_idx = var_db->mapVariableAndContextToIndex(d_u_var, ins_hier_integrator->getCurrentContext());
    const int us_new_idx = var_db->mapVariableAndContextToIndex(d_u_var, ins_hier_integrator->getNewContext());
    const int p_new_idx = var_db->mapVariableAndContextToIndex(d_p_var, ins_hier_integrator->getNewContext());

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
    case BACKWARD_EULER:
        // intentionally blank
        break;
    case MIDPOINT_RULE:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(half_time);
        d_ibn_method_ops->computeLagrangianForce(half_time);
        if (d_cross_links_strategy)
            d_cross_links_strategy->computeLagrangianForce(half_time, IBTK::TimePoint::HALF_TIME);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);

        d_hier_velocity_data_ops->setToScalar(d_fn_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_fn_idx);
        d_ibn_method_ops->spreadForce(
            d_fn_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::fn"), half_time);

        d_hier_velocity_data_ops->setToScalar(d_cross_links_un_idx, 0.0);
        d_hier_velocity_data_ops->setToScalar(d_cross_links_us_idx, 0.0);
        if (d_cross_links_strategy)
        {
            d_u_phys_bdry_op->setPatchDataIndex(d_cross_links_un_idx);
            d_cross_links_strategy->spreadForce(d_cross_links_un_idx,
                                                d_u_phys_bdry_op,
                                                getProlongRefineSchedules(d_object_name + "::fn"),
                                                half_time,
                                                true /*spread_network*/,
                                                IBTK::TimePoint::HALF_TIME);
            d_u_phys_bdry_op->setPatchDataIndex(d_cross_links_us_idx);
            d_cross_links_strategy->spreadForce(d_cross_links_us_idx,
                                                d_u_phys_bdry_op,
                                                getProlongRefineSchedules(d_object_name + "::f"),
                                                half_time,
                                                false /*spread_network*/,
                                                IBTK::TimePoint::HALF_TIME);
        }

        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case TRAPEZOIDAL_RULE:
        if (d_use_structure_predictor || cycle_num > 0)
        {
            // NOTE: We do not re-compute the force unless it could have changed.
            if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
            d_ib_method_ops->computeLagrangianForce(new_time);
            d_ibn_method_ops->computeLagrangianForce(new_time);
            if (d_cross_links_strategy)
                d_cross_links_strategy->computeLagrangianForce(new_time, IBTK::TimePoint::NEW_TIME);
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
            d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
            d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
            d_u_phys_bdry_op->setHomogeneousBc(true);
            d_ib_method_ops->spreadForce(
                d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);

            d_hier_velocity_data_ops->setToScalar(d_fn_idx, 0.0);
            d_u_phys_bdry_op->setPatchDataIndex(d_fn_idx);
            d_ibn_method_ops->spreadForce(
                d_fn_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::fn"), new_time);

            d_hier_velocity_data_ops->setToScalar(d_cross_links_un_idx, 0.0);
            d_hier_velocity_data_ops->setToScalar(d_cross_links_us_idx, 0.0);
            if (d_cross_links_strategy)
            {
                d_u_phys_bdry_op->setPatchDataIndex(d_cross_links_un_idx);
                d_cross_links_strategy->spreadForce(d_cross_links_un_idx,
                                                    d_u_phys_bdry_op,
                                                    getProlongRefineSchedules(d_object_name + "::fn"),
                                                    new_time,
                                                    true /*spread_network*/,
                                                    IBTK::TimePoint::NEW_TIME);
                d_u_phys_bdry_op->setPatchDataIndex(d_cross_links_us_idx);
                d_cross_links_strategy->spreadForce(d_cross_links_us_idx,
                                                    d_u_phys_bdry_op,
                                                    getProlongRefineSchedules(d_object_name + "::f"),
                                                    new_time,
                                                    false /*spread_network*/,
                                                    IBTK::TimePoint::NEW_TIME);
            }

            d_u_phys_bdry_op->setHomogeneousBc(false);
            d_hier_velocity_data_ops->linearSum(d_f_idx, 0.5, d_f_current_idx, 0.5, d_f_idx);
            d_hier_velocity_data_ops->linearSum(d_fn_idx, 0.5, d_fn_current_idx, 0.5, d_fn_idx);
        }
        break;
    default:
        TBOX_ERROR(
            d_object_name
            << "::integrateHierarchy():\n"
            << "  unsupported time stepping type: "
            << IBAMR::enum_to_string<IBAMR::TimeSteppingType>(d_time_stepping_type) << "\n"
            << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Solve the incompressible Navier-Stokes equations.
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): solving the incompressible Navier-Stokes equations\n";
    if (d_current_num_cycles > 1)
    {
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_current_num_cycles == 1);
#endif
        const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
        for (int ins_cycle_num = 0; ins_cycle_num < ins_num_cycles; ++ins_cycle_num)
        {
            d_ins_hier_integrator->integrateHierarchy(current_time, new_time, ins_cycle_num);
        }
    }
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    d_ibn_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
    case BACKWARD_EULER:
        d_hier_velocity_data_ops->copyData(d_u_idx, us_new_idx);
        d_hier_velocity_data_ops->copyData(d_un_idx, un_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        d_u_phys_bdry_op->setPatchDataIndex(d_un_idx);
        d_ibn_method_ops->interpolateVelocity(d_un_idx,
                                              getCoarsenSchedules(d_object_name + "::un::CONSERVATIVE_COARSEN"),
                                              getGhostfillRefineSchedules(d_object_name + "::un"),
                                              new_time);
        break;
    case MIDPOINT_RULE:
    {
        d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, us_cur_idx, 0.5, us_new_idx);
        d_hier_velocity_data_ops->linearSum(d_un_idx, 0.5, un_cur_idx, 0.5, un_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             half_time);
        Pointer<IBMethod> ops = d_ib_method_ops;
        d_u_phys_bdry_op->setPatchDataIndex(d_un_idx);
        d_ibn_method_ops->interpolateVelocity(d_un_idx,
                                              getCoarsenSchedules(d_object_name + "::un::CONSERVATIVE_COARSEN"),
                                              getGhostfillRefineSchedules(d_object_name + "::un"),
                                              half_time);
    }
    break;
    case TRAPEZOIDAL_RULE:
    {
        d_hier_velocity_data_ops->copyData(d_u_idx, us_new_idx);
        d_hier_velocity_data_ops->copyData(d_un_idx, un_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        d_u_phys_bdry_op->setPatchDataIndex(d_un_idx);
        d_ibn_method_ops->interpolateVelocity(d_un_idx,
                                              getCoarsenSchedules(d_object_name + "::un::CONSERVATIVE_COARSEN"),
                                              getGhostfillRefineSchedules(d_object_name + "::un"),
                                              new_time);
    }
    break;
    default:
        TBOX_ERROR(
            d_object_name
            << "::integrateHierarchy():\n"
            << "  unsupported time stepping type: "
            << IBAMR::enum_to_string<IBAMR::TimeSteppingType>(d_time_stepping_type) << "\n"
            << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an updated prediction of the updated positions of the Lagrangian
    // structure.
    if (d_current_num_cycles > 1 && d_current_cycle_num == 0)
    {
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): performing Lagrangian forward-Euler step\n";
        d_ib_method_ops->forwardEulerStep(current_time, new_time);
        d_ibn_method_ops->forwardEulerStep(current_time, new_time);
    }
    else
    {
        switch (d_time_stepping_type)
        {
        case FORWARD_EULER:
            break;
        case BACKWARD_EULER:
            d_ib_method_ops->backwardEulerStep(current_time, new_time);
            d_ibn_method_ops->backwardEulerStep(current_time, new_time);
            break;
        case MIDPOINT_RULE:
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): performing Lagrangian midpoint-rule step\n";
            d_ib_method_ops->midpointStep(current_time, new_time);
            d_ibn_method_ops->midpointStep(current_time, new_time);
            break;
        case TRAPEZOIDAL_RULE:
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): performing Lagrangian trapezoidal-rule step\n";
            d_ib_method_ops->trapezoidalStep(current_time, new_time);
            d_ibn_method_ops->trapezoidalStep(current_time, new_time);
            break;
        default:
            TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                     << "  unsupported time stepping type: "
                                     << IBAMR::enum_to_string<IBAMR::TimeSteppingType>(d_time_stepping_type) << "\n"
                                     << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, "
                                        "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
IBMultiphaseHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                               const double new_time,
                                                               const bool skip_synchronize_new_state_data,
                                                               const int num_cycles)
{
    auto ops = HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(d_u_var, d_hierarchy, true);

    auto velocity_ghost_update = [&](const std::vector<int>& indices)
    {
        using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghostfills;
        for (const int& idx : indices)
        {
            ghostfills.emplace_back(idx,
                                    "CONSERVATIVE_LINEAR_REFINE",
                                    /*use_cf_bdry_interpolation*/ true,
                                    "CONSERVATIVE_COARSEN",
                                    "LINEAR",
                                    false,
                                    d_ins_hier_integrator->getVelocityBoundaryConditions());
        }
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghostfills, d_hierarchy);
        ghost_fill_op.fillData(current_time);
    };

    // The last thing we need to do (before we really postprocess) is update the structure velocity:
    Pointer<MultiphaseStaggeredHierarchyIntegrator> ins_integrator = d_ins_hier_integrator;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int us_new_idx =
        var_db->mapVariableAndContextToIndex(ins_integrator->getSolventVariable(), ins_integrator->getNewContext());
    const int un_new_idx =
        var_db->mapVariableAndContextToIndex(ins_integrator->getNetworkVariable(), ins_integrator->getNewContext());
#ifndef NDEBUG
    ops->setToScalar(d_u_idx, std::numeric_limits<double>::quiet_NaN(), false);
    ops->setToScalar(d_un_idx, std::numeric_limits<double>::quiet_NaN(), false);
#endif
    ops->copyData(d_u_idx, us_new_idx);
    ops->copyData(d_un_idx, un_new_idx);
    if (d_enable_logging)
        plog << d_object_name
             << "::postprocessIntegrateHierarchy(): interpolating Eulerian "
                "velocity to the Lagrangian mesh\n";
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
    d_ib_method_ops->interpolateVelocity(d_u_idx,
                                         getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                         getGhostfillRefineSchedules(d_object_name + "::u"),
                                         new_time);
    d_ibn_method_ops->interpolateVelocity(d_un_idx,
                                          getCoarsenSchedules(d_object_name + "::un::CONSERVATIVE_COARSEN"),
                                          getGhostfillRefineSchedules(d_object_name + "::un"),
                                          new_time);

    // Note IBHierarchyIntegrator calls postprocessIntegrateData on d_ib_method_ops
    d_ibn_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
IBMultiphaseHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                               Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    const IntVector<NDIM> ib_ghosts(d_ib_method_ops->getMinimumGhostCellWidth());

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_un_idx = var_db->registerVariableAndContext(d_un_var, d_ib_context, ib_ghosts);
    d_fn_idx = var_db->registerVariableAndContext(d_fn_var, d_ib_context, ib_ghosts);
    d_cross_links_un_idx = var_db->registerVariableAndContext(d_cross_links_un_var, d_ib_context, ib_ghosts);
    d_cross_links_us_idx = var_db->registerVariableAndContext(d_cross_links_us_var, d_ib_context, ib_ghosts);
    d_ib_data.setFlag(d_un_idx);
    d_ib_data.setFlag(d_fn_idx);
    d_ib_data.setFlag(d_cross_links_un_idx);
    d_ib_data.setFlag(d_cross_links_us_idx);

    if (d_time_stepping_type == FORWARD_EULER || d_time_stepping_type == TRAPEZOIDAL_RULE)
    {
        d_fn_current_idx = var_db->registerClonedPatchDataIndex(d_fn_var, d_fn_idx);
        d_cross_links_current_un_idx = var_db->registerClonedPatchDataIndex(d_cross_links_un_var, d_cross_links_un_idx);
        d_cross_links_current_us_idx = var_db->registerClonedPatchDataIndex(d_cross_links_us_var, d_cross_links_us_idx);
        d_ib_data.setFlag(d_fn_current_idx);
        d_ib_data.setFlag(d_cross_links_current_un_idx);
        d_ib_data.setFlag(d_cross_links_current_us_idx);
    }

    // Setup the fluid solver for IB forces.
    Pointer<MultiphaseStaggeredHierarchyIntegrator> ins_integrator = d_ins_hier_integrator;
#ifndef NDEBUG
    TBOX_ASSERT(ins_integrator);
#endif

    // Note we use setForcingFunctionsScaled() instead of setForcingFunctions() because
    // IBMultiphaseEulerianForceFunction does note scale by volume fraction.
    d_fn_fcn = new IBMultiphaseEulerianForceFunction(this, d_fn_idx);
    d_fs_fcn = new IBMultiphaseEulerianForceFunction(this, d_f_idx);
    ins_integrator->setForcingFunctionsScaled(d_fn_fcn, d_fs_fcn);

    // Note we use setForcingFunctionsScaled because we need to scale by both volume fractions.
    // IBMultiphaseScaledEulerianForceFunction only scales by one volume fraction.
    d_cross_links_un_fcn = new IBMultiphaseEulerianForceFunction(this, d_cross_links_un_idx);
    d_cross_links_us_fcn = new IBMultiphaseEulerianForceFunction(this, d_cross_links_us_idx);
    ins_integrator->setForcingFunctionsScaledByBoth(d_cross_links_un_fcn, d_cross_links_us_fcn);

    // Initialize ibn_method.
    d_ibn_method_ops->registerEulerianVariables();
    d_ibn_method_ops->registerEulerianCommunicationAlgorithms();

    // Create network communications algorithms.
    Pointer<Geometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    d_un_ghostfill_alg = new RefineAlgorithm<NDIM>();
    d_un_ghostfill_op = nullptr;
    d_un_ghostfill_alg->registerRefine(d_un_idx, d_un_idx, d_un_idx, d_un_ghostfill_op);
    std::unique_ptr<RefinePatchStrategy<NDIM>> un_phys_bdry_op_unique(d_u_phys_bdry_op);
    registerGhostfillRefineAlgorithm(d_object_name + "::un", d_un_ghostfill_alg, std::move(un_phys_bdry_op_unique));

    d_un_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_un_coarsen_op = grid_geom->lookupCoarsenOperator(d_un_var, "CONSERVATIVE_COARSEN");
    d_un_coarsen_alg->registerCoarsen(d_un_idx, d_un_idx, d_un_coarsen_op);
    registerCoarsenAlgorithm(d_object_name + "::un::CONSERVATIVE_COARSEN", d_un_coarsen_alg);

    d_fn_prolong_alg = new RefineAlgorithm<NDIM>();
    d_fn_prolong_op = grid_geom->lookupRefineOperator(d_fn_var, "CONSERVATIVE_LINEAR_REFINE");
    d_fn_prolong_alg->registerRefine(d_fn_idx, d_fn_idx, d_fn_idx, d_fn_prolong_op);
    registerProlongRefineAlgorithm(d_object_name + "::fn", d_fn_prolong_alg);

    // Setup tag buffer.
    d_ibn_method_ops->setupTagBuffer(d_tag_buffer, d_gridding_alg);

    d_integrator_is_initialized = true;

    return;
} // initializeHierarchyIntegrator

void
IBMultiphaseHierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                          Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_hierarchy_is_initialized) return;

    IBHierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);

    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        d_ibn_method_ops->beginDataRedistribution(d_hierarchy, d_gridding_alg);
        d_ibn_method_ops->endDataRedistribution(d_hierarchy, d_gridding_alg);
    }

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_un_idx, d_integrator_time);
        level->allocatePatchData(d_scratch_data, d_integrator_time);
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_current_idx = var_db->mapVariableAndContextToIndex(d_un_var, getCurrentContext());
    d_hier_velocity_data_ops->copyData(d_un_idx, un_current_idx);
    const bool initial_time = IBTK::rel_equal_eps(d_integrator_time, d_start_time);
    d_u_phys_bdry_op->setPatchDataIndex(d_un_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
    d_ibn_method_ops->initializePatchHierarchy(hierarchy,
                                               gridding_alg,
                                               d_un_idx,
                                               getCoarsenSchedules(d_object_name + "::un::CONSERVATIVE_COARSEN"),
                                               getGhostfillRefineSchedules(d_object_name + "::un"),
                                               d_integrator_step,
                                               d_integrator_time,
                                               initial_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_un_idx);
        level->deallocatePatchData(d_scratch_data);
    }

    d_hierarchy_is_initialized = true;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBMultiphaseHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
    IBHierarchyIntegrator::regridHierarchyBeginSpecialized();

    d_ibn_method_ops->beginDataRedistribution(d_hierarchy, d_gridding_alg);
} // regridHierarchyBeginSpecialized

void
IBMultiphaseHierarchyIntegrator::regridHierarchyEndSpecialized()
{
    d_ibn_method_ops->endDataRedistribution(d_hierarchy, d_gridding_alg);
    IBHierarchyIntegrator::regridHierarchyEndSpecialized();
} // regridHierarchyEndSpecialized

void
IBMultiphaseHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION", IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION);
    return;
} // putToDatabaseSpecialized

void
IBMultiphaseHierarchyIntegrator::initializeLevelDataSpecialized(Pointer<BasePatchHierarchy<NDIM>> hierarchy,
                                                                const int level_number,
                                                                const double init_data_time,
                                                                const bool can_be_refined,
                                                                const bool initial_time,
                                                                Pointer<BasePatchLevel<NDIM>> old_level,
                                                                const bool allocate_data)
{
    IBHierarchyIntegrator::initializeLevelDataSpecialized(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    d_ibn_method_ops->initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
}

void
IBMultiphaseHierarchyIntegrator::resetHierarchyConfigurationSpecialized(Pointer<BasePatchHierarchy<NDIM>> hierarchy,
                                                                        const int coarsest_level,
                                                                        const int finest_level)
{
    d_ibn_method_ops->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    IBHierarchyIntegrator::resetHierarchyConfigurationSpecialized(hierarchy, coarsest_level, finest_level);
}

/*!
 * Set integer tags to "one" in cells where refinement of the given level
 * should occur according to the magnitude of the fluid vorticity.
 */
void
IBMultiphaseHierarchyIntegrator::applyGradientDetectorSpecialized(Pointer<BasePatchHierarchy<NDIM>> hierarchy,
                                                                  const int level_number,
                                                                  const double error_data_time,
                                                                  const int tag_index,
                                                                  const bool initial_time,
                                                                  const bool uses_richardson_extrapolation_too)
{
    d_ibn_method_ops->applyGradientDetector(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    IBHierarchyIntegrator::applyGradientDetectorSpecialized(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
}

void
IBMultiphaseHierarchyIntegrator::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                     const int workload_data_idx)
{
    d_ibn_method_ops->addWorkloadEstimate(hierarchy, workload_data_idx);
    IBHierarchyIntegrator::addWorkloadEstimate(hierarchy, workload_data_idx);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBMultiphaseHierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
