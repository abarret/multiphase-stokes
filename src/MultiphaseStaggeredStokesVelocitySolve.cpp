/////////////////////////////// INCLUDES /////////////////////////////////////

#include "multiphase/MultiphaseStaggeredStokesVelocitySolve.h"
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
MultiphaseStaggeredStokesVelocitySolve::MultiphaseStaggeredStokesVelocitySolve(const std::string& object_name,
                                                                               bool homogeneous_bc,
                                                                               const MultiphaseParameters& params)
    : LinearOperator(object_name, homogeneous_bc),
      d_default_un_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_un_bc_coef", Pointer<Database>(nullptr))),
      d_default_us_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_us_bc_coef", Pointer<Database>(nullptr))),
      d_un_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef)),
      d_us_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef)),
      d_default_thn_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_thn_bc_coef", Pointer<Database>(nullptr))),
      d_thn_bc_coef(d_default_thn_bc_coef),
      d_os_var(new OutersideVariable<NDIM, double>(d_object_name + "::outerside_variable")),
      d_nc_scr_var(new NodeVariable<NDIM, double>(d_object_name + "::ThnNode", 1, false)),
      d_cc_ndim_var(new CellVariable<NDIM, double>(d_object_name + "::ThnCell", NDIM)),
      d_sc_scr_var(new SideVariable<NDIM, double>(d_object_name + "::VelAvg", 1, false)),
      d_params(params)
{
    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the volume fraction(?)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_default_un_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_un_bc_coef);
        auto p_default_us_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_us_bc_coef);
        p_default_un_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_un_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        p_default_us_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_us_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        auto p_default_thn_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_thn_bc_coef);
        p_default_thn_bc_coef->setBoundarySlope(2 * d, 0.0);
        p_default_thn_bc_coef->setBoundarySlope(2 * d + 1, 0.0);
    }

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    // Check if we've already made this variable
    if (var_db->checkVariableExists(d_object_name + "::outerside_variable"))
        d_os_var = var_db->getVariable(d_object_name + "::outerside_variable");
    d_os_idx =
        var_db->registerVariableAndContext(d_os_var, var_db->getContext(d_object_name + "::CTX"), IntVector<NDIM>(0));

    if (var_db->checkVariableExists(d_object_name + "::ThnNode"))
        d_nc_scr_var = var_db->getVariable(d_object_name + "::ThnNode");
    d_nc_scr_idx = var_db->registerVariableAndContext(d_nc_scr_var, var_db->getContext(d_object_name + "::CTX"));

    if (var_db->checkVariableExists(d_object_name + "::ThnCell"))
        d_cc_ndim_var = var_db->getVariable(d_object_name + "::ThnCell");
    d_cc_ndim_idx = var_db->registerVariableAndContext(
        d_cc_ndim_var, var_db->getContext(d_object_name + "::CTX"), IntVector<NDIM>(1));

    if (var_db->checkVariableExists(d_object_name + "::VelAvg"))
        d_sc_scr_var = var_db->getVariable(d_object_name + "::VelAvg");
    d_sc_scr_idx = var_db->registerVariableAndContext(d_sc_scr_var, var_db->getContext(d_object_name + "::CTX"));

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef),
                       std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef),
                       d_default_thn_bc_coef);

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::TwoFluidStaggeredStokesOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::TwoFluidStaggeredStokesOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::TwoFluidStaggeredStokesOperator::deallocateOperatorState()"););
    return;
} // TwoFluidStaggeredStokesOperator

MultiphaseStaggeredStokesVelocitySolve::~MultiphaseStaggeredStokesVelocitySolve()
{
    deallocateOperatorState();
    delete d_default_un_bc_coef;
    delete d_default_us_bc_coef;
    d_default_un_bc_coef = nullptr;
    d_default_us_bc_coef = nullptr;
    // Remove internal patch index from variable database.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_os_idx);
    return;
} // ~TwoFluidStaggeredStokesOperator

// Probably need another one for second fluid equation
void
MultiphaseStaggeredStokesVelocitySolve::setVelocityPoissonSpecifications(const PoissonSpecifications& coefs)
{
    TBOX_WARNING(d_object_name +
                 "::setVelocityPoissonSpecifications: This function is not used. Use setCandDCoefficients instead.");
    return;
} // setVelocityPoissonSpecifications

void
MultiphaseStaggeredStokesVelocitySolve::setCandDCoefficients(const double C, const double D_u)
{
    d_C = C;
    d_D_u = D_u;
}

void
MultiphaseStaggeredStokesVelocitySolve::setThnIdx(const int thn_idx)
{
    d_thn_idx = thn_idx;
}

void
MultiphaseStaggeredStokesVelocitySolve::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                                                           const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                                                           RobinBcCoefStrategy<NDIM>* thn_bc_coef)
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

    if (thn_bc_coef)
        d_thn_bc_coef = thn_bc_coef;
    else
        d_thn_bc_coef = d_default_thn_bc_coef;
    return;
} // setPhysicalBcCoefs

void
MultiphaseStaggeredStokesVelocitySolve::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_bc_helper = bc_helper;
    return;
} // setPhysicalBoundaryHelper

void
MultiphaseStaggeredStokesVelocitySolve::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBAMR_TIMER_START(t_apply);

    // Get the vector components. These pull out patch data indices
    const int un_idx = x.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = x.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int A_un_idx = y.getComponentDescriptorIndex(0);
    const int A_us_idx = y.getComponentDescriptorIndex(1);
    const int un_scratch_idx = d_x->getComponentDescriptorIndex(0);
    const int us_scratch_idx = d_x->getComponentDescriptorIndex(1);
    const int thn_idx = d_thn_idx;

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = InterpolationTransactionComponent(un_scratch_idx,
                                                             un_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_un_bc_coefs);
    transaction_comps[1] = InterpolationTransactionComponent(us_scratch_idx,
                                                             us_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_us_bc_coefs);

    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time); // Fills in all of the ghost cells
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    {
        // Note that thn ghost cells are always filled under inhomogeneous conditions.
        std::vector<InterpolationTransactionComponent> thn_ghost_comps = { InterpolationTransactionComponent(
            thn_idx,
            "CONSERVATIVE_LINEAR_REFINE",
            USE_CF_INTERPOLATION,
            DATA_COARSEN_TYPE,
            BDRY_EXTRAP_TYPE,
            CONSISTENT_TYPE_2_BDRY,
            d_thn_bc_coef) };
        HierarchyGhostCellInterpolation hier_bdry_fill;
        hier_bdry_fill.initializeOperatorState(thn_ghost_comps, d_hierarchy);
        hier_bdry_fill.setHomogeneousBc(false);
        hier_bdry_fill.fillData(d_solution_time);
    }

    applySpecialized(A_un_idx, A_us_idx, un_scratch_idx, us_scratch_idx, thn_idx);

    if (d_bc_helper)
    {
        d_bc_helper->copyDataAtDirichletBoundaries(A_un_idx, un_scratch_idx);
        d_bc_helper->copyDataAtDirichletBoundaries(A_us_idx, us_scratch_idx);
    }

    {
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp = { ITC(A_us_idx, "NONE", false, "CONSERVATIVE_COARSEN"),
                                             ITC(A_un_idx, "NONE", false, "CONSERVATIVE_COARSEN") };
        HierarchyGhostCellInterpolation ghost_cell_fill;
        ghost_cell_fill.initializeOperatorState(ghost_cell_comp, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
        ghost_cell_fill.fillData(d_new_time);
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
MultiphaseStaggeredStokesVelocitySolve::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                                const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in.cloneVector(in.getName());
    d_b = out.cloneVector(out.getName());

    // Setup operator state.
    d_hierarchy = in.getPatchHierarchy();

    // Allocate scratch data.
    d_x->allocateVectorData();
    const int thn_idx = d_thn_idx;

    // Allocate synchronization variable
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_os_idx)) level->allocatePatchData(d_os_idx);
        if (!level->checkAllocated(d_sc_scr_idx)) level->allocatePatchData(d_sc_scr_idx);
        if (!level->checkAllocated(d_nc_scr_idx)) level->allocatePatchData(d_nc_scr_idx);
        if (!level->checkAllocated(d_cc_ndim_idx)) level->allocatePatchData(d_cc_ndim_idx);
    }

    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    d_os_coarsen_op = grid_geom->lookupCoarsenOperator(d_os_var, "CONSERVATIVE_COARSEN");
    d_os_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_os_coarsen_alg->registerCoarsen(d_b->getComponentDescriptorIndex(0), d_os_idx, d_os_coarsen_op);
    d_os_coarsen_scheds.resize(finest_ln - coarsest_ln);
    for (int dst_ln = coarsest_ln; dst_ln < finest_ln; ++dst_ln)
    {
        Pointer<PatchLevel<NDIM>> src_level = d_hierarchy->getPatchLevel(dst_ln + 1);
        Pointer<PatchLevel<NDIM>> dst_level = d_hierarchy->getPatchLevel(dst_ln);
        d_os_coarsen_scheds[dst_ln] = d_os_coarsen_alg->createSchedule(dst_level, src_level);
    }

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.resize(2);
    d_transaction_comps[0] = InterpolationTransactionComponent(d_x->getComponentDescriptorIndex(0),
                                                               in.getComponentDescriptorIndex(0),
                                                               SC_DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_un_bc_coefs);
    d_transaction_comps[1] = InterpolationTransactionComponent(d_x->getComponentDescriptorIndex(1),
                                                               in.getComponentDescriptorIndex(1),
                                                               SC_DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_us_bc_coefs); // noFillCorners

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_x->getPatchHierarchy());

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

    if (d_bc_helper) d_bc_helper->cacheBcCoefData(d_us_bc_coefs, d_solution_time, d_hierarchy);

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
MultiphaseStaggeredStokesVelocitySolve::deallocateOperatorState()
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

    // Deallocate scratch data.
    d_x->deallocateVectorData();
    d_b->deallocateVectorData();

    // Delete the solution and rhs vectors.
    d_x->resetLevels(d_x->getCoarsestLevelNumber(),
                     std::min(d_x->getFinestLevelNumber(), d_x->getPatchHierarchy()->getFinestLevelNumber()));
    d_x->freeVectorComponents();

    d_b->resetLevels(d_b->getCoarsestLevelNumber(),
                     std::min(d_b->getFinestLevelNumber(), d_b->getPatchHierarchy()->getFinestLevelNumber()));
    d_b->freeVectorComponents();

    d_x.setNull();
    d_b.setNull();

    // Deallocate synchronization variable
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_os_idx)) level->deallocatePatchData(d_os_idx);
        if (level->checkAllocated(d_sc_scr_idx)) level->deallocatePatchData(d_sc_scr_idx);
        if (level->checkAllocated(d_nc_scr_idx)) level->deallocatePatchData(d_nc_scr_idx);
        if (level->checkAllocated(d_cc_ndim_idx)) level->deallocatePatchData(d_cc_ndim_idx);
    }
    d_os_coarsen_scheds.clear();
    d_os_coarsen_alg = nullptr;

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

void
MultiphaseStaggeredStokesVelocitySolve::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
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
        if (d_bc_helper)
        {
            const int un_idx = x->getComponentDescriptorIndex(0);
            const int us_idx = x->getComponentDescriptorIndex(1);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                un_idx, IBTK::invalid_index, d_un_bc_coefs, d_new_time, d_homogeneous_bc);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                us_idx, IBTK::invalid_index, d_us_bc_coefs, d_new_time, d_homogeneous_bc);
        }
        apply(*x, *b);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&y, false), b);
        x->freeVectorComponents();
        b->freeVectorComponents();
    }
    const bool homogeneous_bc = true;
    if (d_bc_helper)
    {
        const int un_idx = y.getComponentDescriptorIndex(0);
        const int us_idx = y.getComponentDescriptorIndex(1);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(
            un_idx, IBTK::invalid_index, d_un_bc_coefs, d_new_time, homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(
            us_idx, IBTK::invalid_index, d_us_bc_coefs, d_new_time, homogeneous_bc);
    }
    return;
} // modifyRhsForBcs

void
MultiphaseStaggeredStokesVelocitySolve::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    if (d_bc_helper)
    {
        const int un_idx = u.getComponentDescriptorIndex(0);
        const int us_idx = u.getComponentDescriptorIndex(1);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(
            un_idx, IBTK::invalid_index, d_un_bc_coefs, d_new_time, d_homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(
            us_idx, IBTK::invalid_index, d_us_bc_coefs, d_new_time, d_homogeneous_bc);
    }
    return;
} // imposeSolBcs

void
MultiphaseStaggeredStokesVelocitySolve::applySpecialized(const int A_un_idx,
                                                         const int A_us_idx,
                                                         const int un_idx,
                                                         const int us_idx,
                                                         const int thn_idx)
{
    // Interpolate and synchronize volume fraction
    d_hier_math_ops->interp(
        d_nc_scr_idx, d_nc_scr_var, true, thn_idx, Pointer<CellVariable<NDIM, double>>(nullptr), nullptr, d_new_time);
    // Interpolate to cell sides
    convert_to_ndim_cc(d_cc_ndim_idx, thn_idx, *d_hierarchy);
    d_hier_math_ops->interp(d_sc_scr_idx, d_sc_scr_var, true, d_cc_ndim_idx, d_cc_ndim_var, nullptr, d_new_time);

    const double xi = d_params.xi;
    const double eta_n = d_params.eta_n;
    const double eta_s = d_params.eta_s;
    const double C = d_C;
    const double D = d_D_u;

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            Pointer<SideData<NDIM, double>> A_un_data = patch->getPatchData(A_un_idx);
            Pointer<SideData<NDIM, double>> A_us_data = patch->getPatchData(A_us_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);

            Pointer<NodeData<NDIM, double>> thn_nc_data = patch->getPatchData(d_nc_scr_idx);
            Pointer<SideData<NDIM, double>> thn_sc_data = patch->getPatchData(d_sc_scr_idx);
            Pointer<CellData<NDIM, double>> thn_cc_data = patch->getPatchData(thn_idx);

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            IntVector<NDIM> xp(1, 0), yp(0, 1);

            for (SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
            {
                const SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)

                CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
                CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
                SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
                SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)
                SideIndex<NDIM> l_y_idx(idx_c_low, 1, 0);    // (i-1,j-1/2)
                SideIndex<NDIM> u_y_idx(idx_c_low, 1, 1);    // (i-1,j+1/2)

                NodeIndex<NDIM> idx_n_l(idx.toCell(1), NodeIndex<NDIM>::LowerLeft);
                NodeIndex<NDIM> idx_n_u(idx.toCell(1), NodeIndex<NDIM>::UpperLeft);
                double thn_lower = (*thn_sc_data)(idx);
                double thn_imhalf_jphalf = (*thn_nc_data)(idx_n_u);
                double thn_imhalf_jmhalf = (*thn_nc_data)(idx_n_l);

                // components of first row (x-component of network vel) of network equation
                double ddx_Thn_dx_un = eta_n / (dx[0] * dx[0]) *
                                       ((*thn_cc_data)(idx_c_up) * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                        (*thn_cc_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - xp)));
                double ddy_Thn_dy_un = eta_n / (dx[1] * dx[1]) *
                                       (thn_imhalf_jphalf * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - yp)));
                double ddy_Thn_dx_vn = eta_n / (dx[1] * dx[0]) *
                                       (thn_imhalf_jphalf * ((*un_data)(upper_y_idx) - (*un_data)(u_y_idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(lower_y_idx) - (*un_data)(l_y_idx)));
                double ddx_Thn_dy_vn = -eta_n / (dx[0] * dx[1]) *
                                       ((*thn_cc_data)(idx_c_up) * ((*un_data)(upper_y_idx) - (*un_data)(lower_y_idx)) -
                                        (*thn_cc_data)(idx_c_low) * ((*un_data)(u_y_idx) - (*un_data)(l_y_idx)));
                double drag_n = -xi * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
                const double A_un = C * thn_lower * (*un_data)(idx) +
                                    D * (ddx_Thn_dx_un + ddy_Thn_dy_un + ddy_Thn_dx_vn + ddx_Thn_dy_vn);
                const double xi_un = D * drag_n;
                (*A_un_data)(idx) = A_un + xi_un;

                // solvent equation
                double ddx_Ths_dx_us =
                    eta_s / (dx[0] * dx[0]) *
                    (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                     convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - xp)));
                double ddy_Ths_dy_us = eta_s / (dx[1] * dx[1]) *
                                       (convertToThs(thn_imhalf_jphalf) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                        convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - yp)));
                double ddy_Ths_dx_vs =
                    eta_s / (dx[1] * dx[0]) *
                    (convertToThs(thn_imhalf_jphalf) * ((*us_data)(upper_y_idx) - (*us_data)(u_y_idx)) -
                     convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_y_idx) - (*us_data)(l_y_idx)));
                double ddx_Ths_dy_vs =
                    -eta_s / (dx[0] * dx[1]) *
                    (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(upper_y_idx) - (*us_data)(lower_y_idx)) -
                     convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(u_y_idx) - (*us_data)(l_y_idx)));

                const double A_us = D * (ddx_Ths_dx_us + ddy_Ths_dy_us + ddy_Ths_dx_vs + ddx_Ths_dy_vs) +
                                    C * convertToThs(thn_lower) * (*us_data)(idx);
                (*A_us_data)(idx) = A_un + A_us;
            }

            for (SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
            {
                const SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)

                CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
                CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
                SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
                SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
                SideIndex<NDIM> l_x_idx(idx_c_low, 0, 0);    // (i-1/2,j-1)
                SideIndex<NDIM> u_x_idx(idx_c_low, 0, 1);    // (i+1/2,j-1)

                NodeIndex<NDIM> idx_n_l(idx.toCell(1), NodeIndex<NDIM>::LowerLeft);
                NodeIndex<NDIM> idx_n_u(idx.toCell(1), NodeIndex<NDIM>::LowerRight);
                double thn_lower = (*thn_sc_data)(idx);
                double thn_iphalf_jmhalf = (*thn_nc_data)(idx_n_u);
                double thn_imhalf_jmhalf = (*thn_nc_data)(idx_n_l);

                // components of second row (y-component of network vel) of network equation
                double ddy_Thn_dy_un = eta_n / (dx[1] * dx[1]) *
                                       ((*thn_cc_data)(idx_c_up) * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                        (*thn_cc_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - yp)));
                double ddx_Thn_dx_un = eta_n / (dx[0] * dx[0]) *
                                       (thn_iphalf_jmhalf * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - xp)));
                double ddx_Thn_dy_vn = eta_n / (dx[1] * dx[0]) *
                                       (thn_iphalf_jmhalf * ((*un_data)(upper_x_idx) - (*un_data)(u_x_idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(lower_x_idx) - (*un_data)(l_x_idx)));
                double ddy_Thn_dx_vn = -eta_n / (dx[0] * dx[1]) *
                                       ((*thn_cc_data)(idx_c_up) * ((*un_data)(upper_x_idx) - (*un_data)(lower_x_idx)) -
                                        (*thn_cc_data)(idx_c_low) * ((*un_data)(u_x_idx) - (*un_data)(l_x_idx)));
                double drag_n = -xi * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));

                const double A_un = D * (ddy_Thn_dy_un + ddx_Thn_dx_un + ddx_Thn_dy_vn + ddy_Thn_dx_vn) +
                                    C * thn_lower * (*un_data)(idx);
                const double xi_un = D * drag_n;
                (*A_un_data)(idx) = A_un + xi_un;

                // Solvent equation
                double ddy_Ths_dy_us =
                    eta_s / (dx[1] * dx[1]) *
                    (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                     convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - yp)));
                double ddx_Ths_dx_us = eta_s / (dx[0] * dx[0]) *
                                       (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                        convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - xp)));
                double ddx_Ths_dy_vs =
                    eta_s / (dx[1] * dx[0]) *
                    (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(upper_x_idx) - (*us_data)(u_x_idx)) -
                     convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_x_idx) - (*us_data)(l_x_idx)));
                double ddy_Ths_dx_vs =
                    -eta_s / (dx[0] * dx[1]) *
                    (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(upper_x_idx) - (*us_data)(lower_x_idx)) -
                     convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(u_x_idx) - (*us_data)(l_x_idx)));
                const double A_us = D * (ddy_Ths_dy_us + ddx_Ths_dx_us + ddx_Ths_dy_vs + ddy_Ths_dx_vs) +
                                    C * convertToThs(thn_lower) * (*us_data)(idx);

                (*A_us_data)(idx) = A_un + A_us;
            }
        }
    }
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
