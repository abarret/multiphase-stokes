/////////////////////////////// INCLUDES /////////////////////////////////////

#include "multiphase/MultiphaseStaggeredStokesOperator.h"
#include "multiphase/fd_operators.h"
#include "multiphase/utility_functions.h"

#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/SideNoCornersFillPattern.h"

#include "CellVariable.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "VariableFillPattern.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace multiphase
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string CC_DATA_REFINE_TYPE = "NONE"; // how to fill in fine cells from coarse cells, how to fill ghost
                                                       // cells on refine patch
static const std::string SC_DATA_REFINE_TYPE = "NONE"; // how to fill in fine cells from coarse cells, how to fill ghost
                                                       // cells on refine patch
static const bool USE_CF_INTERPOLATION = true;         // Refine Patch Strategy: CartSideDoubleQuadraticCFInterpolation.
static const std::string DATA_COARSEN_TYPE =
    "CUBIC_COARSEN"; // going from fine to coarse. fill in coarse cells by whatever is in the fine cells. synchronizing
                     // the hierarchies

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "QUADRATIC"; // these operations are all in IBAMR

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////
MultiphaseStaggeredStokesOperator::MultiphaseStaggeredStokesOperator(
    const std::string& object_name,
    bool homogeneous_bc,
    const MultiphaseParameters& params,
    const std::unique_ptr<VolumeFractionDataManager>& thn_manager)
    : LinearOperator(object_name, homogeneous_bc),
      d_default_un_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_un_bc_coef", Pointer<Database>(nullptr))),
      d_default_us_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_us_bc_coef", Pointer<Database>(nullptr))),
      d_un_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef)),
      d_us_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef)),
      d_default_P_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(nullptr))),
      d_P_bc_coef(d_default_P_bc_coef),
      d_os_var(new OutersideVariable<NDIM, double>(d_object_name + "::outerside_variable")),
      d_params(params),
      d_thn_manager(thn_manager)
{
    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_default_un_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_un_bc_coef);
        auto p_default_us_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_us_bc_coef);
        p_default_un_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_un_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        p_default_us_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_us_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        auto p_default_P_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_P_bc_coef);
        p_default_P_bc_coef->setBoundarySlope(2 * d, 0.0);
        p_default_P_bc_coef->setBoundarySlope(2 * d + 1, 0.0);
    }

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    // Check if we've already made this variable
    if (var_db->checkVariableExists(d_object_name + "::outerside_variable"))
        d_os_var = var_db->getVariable(d_object_name + "::outerside_variable");
    d_os_idx =
        var_db->registerVariableAndContext(d_os_var, var_db->getContext(d_object_name + "::CTX"), IntVector<NDIM>(0));

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef),
                       std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef),
                       d_default_P_bc_coef);

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply =
                      TimerManager::getManager()->getTimer("multiphase::MultiphaseStaggeredStokesOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "multiphase::MultiphaseStaggeredStokesOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "multiphase::MultiphaseStaggeredStokesOperator::deallocateOperatorState()"););
    return;
} // TwoFluidStaggeredStokesOperator

MultiphaseStaggeredStokesOperator::~MultiphaseStaggeredStokesOperator()
{
    deallocateOperatorState();
    delete d_default_un_bc_coef;
    delete d_default_us_bc_coef;
    d_default_un_bc_coef = nullptr;
    d_default_us_bc_coef = nullptr;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = nullptr;
    // Remove internal patch index from variable database.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_os_idx);
    return;
} // ~TwoFluidStaggeredStokesOperator

// Probably need another one for second fluid equation
void
MultiphaseStaggeredStokesOperator::setVelocityPoissonSpecifications(const PoissonSpecifications& coefs)
{
    TBOX_WARNING(d_object_name +
                 "::setVelocityPoissonSpecifications: This function is not used. Use setCandDCoefficients instead.");
    return;
} // setVelocityPoissonSpecifications

void
MultiphaseStaggeredStokesOperator::setCandDCoefficients(const double C,
                                                        const double D_u,
                                                        const double D_p,
                                                        const double D_div)
{
    d_C = C;
    d_D_u = D_u;
    d_D_p = D_p;
    d_D_div = D_div;
}

void
MultiphaseStaggeredStokesOperator::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                                                      const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                                                      RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(un_bc_coefs.size() == NDIM);
    TBOX_ASSERT(us_bc_coefs.size() == NDIM);
#endif
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (un_bc_coefs[d])
        {
            d_un_bc_coefs[d] = un_bc_coefs[d];
        }
        else
        {
            d_un_bc_coefs[d] = d_default_un_bc_coef;
        }
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (us_bc_coefs[d])
        {
            d_us_bc_coefs[d] = us_bc_coefs[d];
        }
        else
        {
            d_us_bc_coefs[d] = d_default_us_bc_coef;
        }
    }

    if (P_bc_coef)
    {
        d_P_bc_coef = P_bc_coef;
    }
    else
    {
        d_P_bc_coef = d_default_P_bc_coef;
    }
    return;
} // setPhysicalBcCoefs

void
MultiphaseStaggeredStokesOperator::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_un_helper,
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_us_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_un_helper);
    TBOX_ASSERT(bc_us_helper);
#endif
    d_bc_un_helper = bc_un_helper;
    d_bc_us_helper = bc_us_helper;
    return;
} // setPhysicalBoundaryHelper

void
MultiphaseStaggeredStokesOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBAMR_TIMER_START(t_apply);

    // Get the vector components. These pull out patch data indices
    const int un_idx = x.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = x.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int P_idx = x.getComponentDescriptorIndex(2);  // pressure
    const int A_un_idx = y.getComponentDescriptorIndex(0);
    const int A_us_idx = y.getComponentDescriptorIndex(1);
    const int A_P_idx = y.getComponentDescriptorIndex(2);

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(3);
    transaction_comps[0] = InterpolationTransactionComponent(un_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_un_bc_coefs,
                                                             nullptr,
                                                             "QUADRATIC");
    transaction_comps[1] = InterpolationTransactionComponent(us_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_us_bc_coefs,
                                                             nullptr,
                                                             "QUADRATIC");
    transaction_comps[2] = InterpolationTransactionComponent(P_idx,
                                                             CC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef,
                                                             nullptr,
                                                             "QUADRATIC");

    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time); // Fills in all of the ghost cells
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    applySpecialized(A_P_idx, A_un_idx, A_us_idx, P_idx, un_idx, us_idx);

    if (d_bc_un_helper)
    {
        d_bc_un_helper->copyDataAtDirichletBoundaries(A_un_idx, un_idx);
    }
    if (d_bc_us_helper)
    {
        d_bc_us_helper->copyDataAtDirichletBoundaries(A_us_idx, us_idx);
    }

    {
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp = { ITC(A_us_idx, "NONE", false, "CONSERVATIVE_COARSEN"),
                                             ITC(A_un_idx, "NONE", false, "CONSERVATIVE_COARSEN") };
        HierarchyGhostCellInterpolation ghost_cell_fill;
        ghost_cell_fill.initializeOperatorState(ghost_cell_comp, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
        ghost_cell_fill.fillData(d_solution_time);
    }

    auto sync_fcn = [&](const int dst_idx) -> void
    {
        for (int ln = d_hierarchy->getFinestLevelNumber(); ln > 0; --ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<SideData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
                Pointer<OutersideData<NDIM, double>> os_data = patch->getPatchData(d_os_idx);
                os_data->copy(*dst_data);
            }
            Pointer<CoarsenAlgorithm<NDIM>> coarsen_alg = new CoarsenAlgorithm<NDIM>();
            coarsen_alg->registerCoarsen(dst_idx, d_os_idx, d_os_coarsen_op);
            coarsen_alg->resetSchedule(d_os_coarsen_scheds[ln - 1]);
            d_os_coarsen_scheds[ln - 1]->coarsenData();
            d_os_coarsen_alg->resetSchedule(d_os_coarsen_scheds[ln - 1]);
        }
    };

    sync_fcn(A_us_idx);
    sync_fcn(A_un_idx);

    IBAMR_TIMER_STOP(t_apply);
    return;
} // apply

void
MultiphaseStaggeredStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup operator state.
    d_hierarchy = in.getPatchHierarchy();

    // Allocate synchronization variable
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    allocate_patch_data({ d_os_idx }, d_hierarchy, d_solution_time, coarsest_ln, finest_ln);

    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    d_os_coarsen_op = grid_geom->lookupCoarsenOperator(d_os_var, "CONSERVATIVE_COARSEN");
    d_os_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_os_coarsen_alg->registerCoarsen(out.getComponentDescriptorIndex(0), d_os_idx, d_os_coarsen_op);
    d_os_coarsen_scheds.resize(finest_ln - coarsest_ln);
    for (int dst_ln = coarsest_ln; dst_ln < finest_ln; ++dst_ln)
    {
        Pointer<PatchLevel<NDIM>> src_level = d_hierarchy->getPatchLevel(dst_ln + 1);
        Pointer<PatchLevel<NDIM>> dst_level = d_hierarchy->getPatchLevel(dst_ln);
        d_os_coarsen_scheds[dst_ln] = d_os_coarsen_alg->createSchedule(dst_level, src_level);
    }

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.resize(3);
    d_transaction_comps[0] = InterpolationTransactionComponent(in.getComponentDescriptorIndex(0),
                                                               SC_DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_un_bc_coefs);
    d_transaction_comps[1] = InterpolationTransactionComponent(in.getComponentDescriptorIndex(1),
                                                               SC_DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_us_bc_coefs);
    d_transaction_comps[2] = InterpolationTransactionComponent(in.getComponentDescriptorIndex(2),
                                                               CC_DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_P_bc_coef); // noFillCorners

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, in.getPatchHierarchy());

    // Initialize hierarchy math ops object.
    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops = new HierarchyMathOps(d_object_name + "::HierarchyMathOps",
                                               in.getPatchHierarchy(),
                                               in.getCoarsestLevelNumber(),
                                               in.getFinestLevelNumber());
    }
#if !defined(NDEBUG)
    else
    {
        TBOX_ASSERT(d_hier_math_ops);
    }
#endif

    if (d_bc_un_helper) d_bc_un_helper->cacheBcCoefData(d_un_bc_coefs, d_solution_time, d_hierarchy);
    if (d_bc_us_helper) d_bc_us_helper->cacheBcCoefData(d_us_bc_coefs, d_solution_time, d_hierarchy);

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
MultiphaseStaggeredStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();
    d_un_fill_pattern.setNull();
    d_us_fill_pattern.setNull();
    d_P_fill_pattern.setNull();

    // Deallocate synchronization variable
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    deallocate_patch_data({ d_os_idx }, d_hierarchy, coarsest_ln, finest_ln);
    d_os_coarsen_scheds.clear();
    d_os_coarsen_alg = nullptr;

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

void
MultiphaseStaggeredStokesOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
{
    if (!d_homogeneous_bc)
    {
        // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
        // inhomogeneous boundary conditions.
        Pointer<SAMRAIVectorReal<NDIM, double>> x = y.cloneVector("");
        Pointer<SAMRAIVectorReal<NDIM, double>> b = y.cloneVector("");
        x->allocateVectorData();
        b->allocateVectorData();
        x->setToScalar(0.0);
        if (d_bc_un_helper)
        {
            const int un_idx = x->getComponentDescriptorIndex(0);
            const int P_idx = x->getComponentDescriptorIndex(2);
            d_bc_un_helper->enforceNormalVelocityBoundaryConditions(
                un_idx, P_idx, d_un_bc_coefs, d_solution_time, d_homogeneous_bc);
        }
        if (d_bc_us_helper)
        {
            const int us_idx = x->getComponentDescriptorIndex(1);
            const int P_idx = x->getComponentDescriptorIndex(2);
            d_bc_us_helper->enforceNormalVelocityBoundaryConditions(
                us_idx, P_idx, d_us_bc_coefs, d_solution_time, d_homogeneous_bc);
        }
        apply(*x, *b);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&y, false), b);
        x->freeVectorComponents();
        b->freeVectorComponents();
    }
    const bool homogeneous_bc = true;
    if (d_bc_un_helper)
    {
        const int un_idx = y.getComponentDescriptorIndex(0);
        const int P_idx = y.getComponentDescriptorIndex(2);
        d_bc_un_helper->enforceNormalVelocityBoundaryConditions(
            un_idx, P_idx, d_un_bc_coefs, d_solution_time, homogeneous_bc);
    }
    if (d_bc_us_helper)
    {
        const int us_idx = y.getComponentDescriptorIndex(1);
        const int P_idx = y.getComponentDescriptorIndex(2);
        d_bc_us_helper->enforceNormalVelocityBoundaryConditions(
            us_idx, P_idx, d_us_bc_coefs, d_solution_time, homogeneous_bc);
    }
    return;
} // modifyRhsForBcs

void
MultiphaseStaggeredStokesOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    if (d_bc_un_helper)
    {
        const int un_idx = u.getComponentDescriptorIndex(0);
        const int P_idx = u.getComponentDescriptorIndex(2);
        d_bc_un_helper->enforceNormalVelocityBoundaryConditions(
            un_idx, P_idx, d_un_bc_coefs, d_solution_time, d_homogeneous_bc);
    }
    if (d_bc_us_helper)
    {
        const int us_idx = u.getComponentDescriptorIndex(1);
        const int P_idx = u.getComponentDescriptorIndex(2);
        d_bc_us_helper->enforceNormalVelocityBoundaryConditions(
            us_idx, P_idx, d_us_bc_coefs, d_solution_time, d_homogeneous_bc);
    }
    return;
} // imposeSolBcs

void
MultiphaseStaggeredStokesOperator::applySpecialized(const int A_P_idx,
                                                    const int A_un_idx,
                                                    const int A_us_idx,
                                                    const int p_idx,
                                                    const int un_idx,
                                                    const int us_idx)
{
    const int thn_idx = d_thn_manager->getCellIndex();
    const int thn_sc_idx = d_thn_manager->getSideIndex();
    const int thn_nc_idx = d_thn_manager->getNodeIndex();
    // Compute volume average velocity and compute divergence.
    pre_div_interp(A_un_idx, thn_idx, un_idx, us_idx, d_hierarchy);
    d_hier_math_ops->div(A_P_idx,
                         Pointer<CellVariable<NDIM, double>>(nullptr),
                         d_D_div,
                         A_un_idx,
                         Pointer<SideVariable<NDIM, double>>(nullptr),
                         nullptr,
                         d_solution_time,
                         true);

    if (d_params.isVariableDrag())
        accumulateMomentumForcesVariableDrag(*d_hierarchy,
                                             A_un_idx,
                                             A_us_idx,
                                             p_idx,
                                             un_idx,
                                             us_idx,
                                             thn_idx,
                                             thn_nc_idx,
                                             thn_sc_idx,
                                             d_params,
                                             d_C,
                                             d_D_u,
                                             d_D_p);
    else
        accumulateMomentumForcesConstantCoefficient(*d_hierarchy,
                                                    A_un_idx,
                                                    A_us_idx,
                                                    p_idx,
                                                    un_idx,
                                                    us_idx,
                                                    thn_idx,
                                                    thn_nc_idx,
                                                    thn_sc_idx,
                                                    d_params,
                                                    d_C,
                                                    d_D_u,
                                                    d_D_p);
}

} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
