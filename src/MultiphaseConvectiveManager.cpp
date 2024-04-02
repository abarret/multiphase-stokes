
// Local includes
#include "multiphase/MultiphaseConvectiveManager.h"
#include "multiphase/utility_functions.h"

#include <ibamr/INSStaggeredConvectiveOperatorManager.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>

#include <memory>

// Fortran routines
extern "C"
{
    void vc_navier_stokes_upwind_quantity2d_(const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const double*,
                                             const double*,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const double*,
                                             const double*,
                                             const int&,
                                             const int&,
                                             double*,
                                             double*,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const double*,
                                             const double*,
                                             const int&,
                                             const int&,
                                             double*,
                                             double*);

    void convect_derivative2d_(const double*,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const double*,
                               const double*,
                               const double*,
                               const double*,
                               const int&,
                               const int&,
                               double*);

    void navier_stokes_interp_comps2d_(const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const double*,
                                       const double*,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       double*,
                                       double*,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       const int&,
                                       double*,
                                       double*);

    void vc_navier_stokes_cui_quantity2d_(const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const double*,
                                          const double*,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const double*,
                                          const double*,
                                          const int&,
                                          const int&,
                                          double*,
                                          double*,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const int&,
                                          const double*,
                                          const double*,
                                          const int&,
                                          const int&,
                                          double*,
                                          double*);

    void vc_navier_stokes_mgamma_quantity2d_(const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const double*,
                                             const double*,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const double*,
                                             const double*,
                                             const int&,
                                             const int&,
                                             double*,
                                             double*,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const int&,
                                             const double*,
                                             const double*,
                                             const int&,
                                             const int&,
                                             double*,
                                             double*);

    void vc_navier_stokes_fbics_quantity2d_(const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const double*,
                                            const double*,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const double*,
                                            const double*,
                                            const int&,
                                            const int&,
                                            double*,
                                            double*,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const int&,
                                            const double*,
                                            const double*,
                                            const int&,
                                            const int&,
                                            double*,
                                            double*);
}
namespace multiphase
{
// Function to convert limiter to required ghost cell width
int
get_limiter_gcw(IBAMR::LimiterType limiter)
{
    switch (limiter)
    {
    case LimiterType::UPWIND:
        return 2;
        break;
    case LimiterType::CUI:
    case LimiterType::FBICS:
    case LimiterType::MGAMMA:
        return 3;
        break;
    case LimiterType::PPM:
        return 4;
        break;
    default:
        TBOX_ERROR("Limiter type " << IBAMR::enum_to_string(limiter) << " not supported!\n");
        break;
    }
    // Should not reach this statement.
    return -1;
}

template <typename VarType>
Pointer<VarType>
get_var(const std::string& var_name)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VarType> var;
    if (var_db->checkVariableExists(var_name))
        var = var_db->getVariable(var_name);
    else
        var = new VarType(var_name);
    return var;
}

MultiphaseConvectiveManager::MultiphaseConvectiveManager(std::string object_name,
                                                         Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                         Pointer<Database> input_db,
                                                         const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                                                         const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                                                         RobinBcCoefStrategy<NDIM>* thn_bc_coef)
    : d_object_name(std::move(object_name)),
      d_hierarchy(hierarchy),
      d_un_bc_coefs(un_bc_coefs),
      d_us_bc_coefs(us_bc_coefs),
      d_thn_bc_coef(thn_bc_coef)
{
    d_limiter = IBAMR::string_to_enum<IBAMR::LimiterType>(input_db->getString("limiter_type"));
    d_bdry_interp_order = input_db->getStringWithDefault("bdry_interp_order", d_bdry_interp_order);
    commonConstructor();
}

MultiphaseConvectiveManager::MultiphaseConvectiveManager(std::string object_name,
                                                         Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                         LimiterType limiter_type,
                                                         const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                                                         const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                                                         RobinBcCoefStrategy<NDIM>* thn_bc_coef)
    : d_object_name(std::move(object_name)),
      d_hierarchy(hierarchy),
      d_limiter(limiter_type),
      d_un_bc_coefs(un_bc_coefs),
      d_us_bc_coefs(us_bc_coefs),
      d_thn_bc_coef(thn_bc_coef)
{
    commonConstructor();
}

void
MultiphaseConvectiveManager::commonConstructor()
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> network_ctx = var_db->getContext(d_object_name + "::NetworkCTX");
    Pointer<VariableContext> solvent_ctx = var_db->getContext(d_object_name + "::SolventCTX");
    const int gcw = get_limiter_gcw(d_limiter);

    // Grab the variables.
    d_mom_var = get_var<SideVariable<NDIM, double>>(d_object_name + "::Mom_var");
    d_N_var = get_var<SideVariable<NDIM, double>>(d_object_name + "::N_var");
    d_N0_var = get_var<SideVariable<NDIM, double>>(d_object_name + "::N0_var");
    d_U_var = get_var<SideVariable<NDIM, double>>(d_object_name + "::U_var");
    d_thn_var = get_var<CellVariable<NDIM, double>>(d_object_name + "::thn_var");

    // Create variable indices
    d_mom_un_idx = var_db->registerVariableAndContext(d_mom_var, network_ctx, gcw);
    d_mom_us_idx = var_db->registerVariableAndContext(d_mom_var, solvent_ctx, gcw);
    d_N_un_idx = var_db->registerVariableAndContext(d_N_var, network_ctx);
    d_N_us_idx = var_db->registerVariableAndContext(d_N_var, solvent_ctx);
    d_N0_un_idx = var_db->registerVariableAndContext(d_N0_var, network_ctx);
    d_N0_us_idx = var_db->registerVariableAndContext(d_N0_var, solvent_ctx);
    d_un_scr_idx = var_db->registerVariableAndContext(d_U_var, network_ctx, gcw);
    d_us_scr_idx = var_db->registerVariableAndContext(d_U_var, solvent_ctx, gcw);
    d_thn_scr_idx = var_db->registerVariableAndContext(d_thn_var, network_ctx, gcw + 1);
}

MultiphaseConvectiveManager::~MultiphaseConvectiveManager()
{
    if (getIsAllocated()) deallocateData();
    // Note that we can not remove indices from the variable database here. The GridGeometry object in the
    // PatchHierarchy requires that the maximum ghost cell width of all the patch indices be constant, and removing a
    // patch index could change the maximum ghost cell width (especially because the ghost cell widths used here are so
    // large).
}

void
MultiphaseConvectiveManager::setBoundaryConditions(const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                                                   const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                                                   RobinBcCoefStrategy<NDIM>* thn_bc_coef)
{
    d_un_bc_coefs = un_bc_coefs;
    d_us_bc_coefs = us_bc_coefs;
    d_thn_bc_coef = thn_bc_coef;
}

void
MultiphaseConvectiveManager::deallocateData()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    deallocate_patch_data({ d_mom_un_idx,
                            d_mom_us_idx,
                            d_N_un_idx,
                            d_N_us_idx,
                            d_N0_un_idx,
                            d_N0_us_idx,
                            d_un_scr_idx,
                            d_us_scr_idx,
                            d_thn_scr_idx },
                          d_hierarchy,
                          coarsest_ln,
                          finest_ln);

    d_thn_ghost_fill.deallocateOperatorState();
    d_u_ghost_fill.deallocateOperatorState();

    d_is_allocated = false;
}

void
MultiphaseConvectiveManager::resetData()
{
    d_is_initial_approximation_filled = false;
}

void
MultiphaseConvectiveManager::allocateData(const double time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    if (!getIsAllocated())
    {
        allocate_patch_data({ d_mom_un_idx,
                              d_mom_us_idx,
                              d_N_un_idx,
                              d_N_us_idx,
                              d_N0_un_idx,
                              d_N0_us_idx,
                              d_un_scr_idx,
                              d_us_scr_idx,
                              d_thn_scr_idx },
                            d_hierarchy,
                            time,
                            coarsest_ln,
                            finest_ln);

        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> thn_ghost_comp = { ITC(d_thn_scr_idx,
                                                "CONSERVATIVE_LINEAR_REFINE",
                                                true,
                                                "NONE",
                                                "LINEAR",
                                                false,
                                                d_thn_bc_coef,
                                                nullptr,
                                                d_bdry_interp_order) };
        d_thn_ghost_fill.initializeOperatorState(thn_ghost_comp, d_hierarchy, coarsest_ln, finest_ln);

        std::vector<ITC> u_ghost_fill_itc = { ITC(d_un_scr_idx,
                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                  true,
                                                  "NONE",
                                                  "LINEAR",
                                                  false,
                                                  d_un_bc_coefs,
                                                  nullptr,
                                                  d_bdry_interp_order),
                                              ITC(d_us_scr_idx,
                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                  true,
                                                  "NONE",
                                                  "LINEAR",
                                                  false,
                                                  d_us_bc_coefs,
                                                  nullptr,
                                                  d_bdry_interp_order) };
        d_u_ghost_fill.initializeOperatorState(u_ghost_fill_itc, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());

        // Allocate d_hier_sc_data_ops
        d_hier_sc_data_ops = std::make_unique<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double>>(d_hierarchy);

        // Allocate d_hier_cc_data_ops
        d_hier_cc_data_ops = std::make_unique<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double>>(d_hierarchy);

        d_is_allocated = true;
    }
}

void
MultiphaseConvectiveManager::approximateConvectiveOperator(const int dst_un_idx,
                                                           const int dst_us_idx,
                                                           IBAMR::TimeSteppingType ts_type,
                                                           const double current_time,
                                                           const double new_time,
                                                           const int un_cur_idx,
                                                           const int us_cur_idx,
                                                           const int thn_cur_idx,
                                                           const int un_new_idx,
                                                           const int us_new_idx,
                                                           const int thn_new_idx)
{
    approximateConvectiveOperator(
        ts_type, current_time, new_time, un_cur_idx, us_cur_idx, thn_cur_idx, un_new_idx, us_new_idx, thn_new_idx);
    fillWithConvectiveOperator(dst_un_idx, dst_us_idx);
}

void
MultiphaseConvectiveManager::approximateConvectiveOperator(IBAMR::TimeSteppingType ts_type,
                                                           const double current_time,
                                                           const double new_time,
                                                           const int un_cur_idx,
                                                           const int us_cur_idx,
                                                           const int thn_cur_idx,
                                                           const int un_new_idx,
                                                           const int us_new_idx,
                                                           const int thn_new_idx)
{
    allocateData(current_time);
    switch (ts_type)
    {
    case FORWARD_EULER:
        approximateForwardEuler(current_time, new_time, un_cur_idx, us_cur_idx, thn_cur_idx);
        break;
    case TRAPEZOIDAL_RULE:
        approximateTrapezoidalRule(
            current_time, new_time, un_cur_idx, us_cur_idx, thn_cur_idx, un_new_idx, us_new_idx, thn_new_idx);
        break;
    case MIDPOINT_RULE:
        approximateMidpointRule(
            current_time, new_time, un_cur_idx, us_cur_idx, thn_cur_idx, un_new_idx, us_new_idx, thn_new_idx);
        break;
    default:
        TBOX_ERROR("Unknown time stepping type "
                   << IBAMR::enum_to_string(ts_type) << "\n"
                   << "Valid options are FORWARD_EULER, TRAPEZOIDAL_RULE, or MIDPOINT_RULE.\n");
    }
}

void
MultiphaseConvectiveManager::fillWithConvectiveOperator(const int dst_un_idx, const int dst_us_idx)
{
    if (dst_un_idx != IBTK::invalid_index) d_hier_sc_data_ops->copyData(dst_un_idx, d_N_un_idx);
    if (dst_us_idx != IBTK::invalid_index) d_hier_sc_data_ops->copyData(dst_us_idx, d_N_us_idx);
}

std::pair<int, int>
MultiphaseConvectiveManager::getConvectiveIndices() const
{
    return std::make_pair(d_N_un_idx, d_N_us_idx);
}

bool
MultiphaseConvectiveManager::getIsAllocated() const
{
    return d_is_allocated;
}

void
MultiphaseConvectiveManager::approximateOperator(const int dst_un_idx,
                                                 const int dst_us_idx,
                                                 const double eval_time,
                                                 const int un_idx,
                                                 const int us_idx,
                                                 const int thn_idx)
{
    // Fill in ghost cells for thn. Because thn may have changed indices, we need to reinitialize the ghost filling
    // routines.
    d_hier_cc_data_ops->copyData(d_thn_scr_idx, thn_idx);
    d_thn_ghost_fill.fillData(eval_time);

    // Fill in velocity ghost cells. Needed to compute staggered control volume velocities
    d_hier_sc_data_ops->copyData(d_un_scr_idx, un_idx);
    d_hier_sc_data_ops->copyData(d_us_scr_idx, us_idx);
    d_u_ghost_fill.fillData(eval_time);

    // Fill in N0 approximations and N approximations.
    // First find the respective momentums. Note that this should also fill in ghost cells of the momentum operator.
    findNetworkMomentum(d_mom_un_idx, d_thn_scr_idx, d_un_scr_idx);
    findSolventMomentum(d_mom_us_idx, d_thn_scr_idx, d_us_scr_idx);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(d_un_scr_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(d_us_scr_idx);
            Pointer<SideData<NDIM, double>> mom_un_data = patch->getPatchData(d_mom_un_idx);
            Pointer<SideData<NDIM, double>> mom_us_data = patch->getPatchData(d_mom_us_idx);
            Pointer<SideData<NDIM, double>> N_un_data = patch->getPatchData(dst_un_idx);
            Pointer<SideData<NDIM, double>> N_us_data = patch->getPatchData(dst_us_idx);

            // Allocate the interpolation and advection velocity data.
            std::array<std::unique_ptr<FaceData<NDIM, double>>, NDIM> interp_mom_data;
            std::array<std::unique_ptr<FaceData<NDIM, double>>, NDIM> u_adv_data;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                static const IntVector<NDIM> ghosts = 1;
                interp_mom_data[axis] = std::make_unique<FaceData<NDIM, double>>(side_box, 1, ghosts);
                u_adv_data[axis] = std::make_unique<FaceData<NDIM, double>>(side_box, 1, ghosts);
            }

            // Compute advection velocity
            computeAdvectionVelocity(patch, u_adv_data, *un_data);

            // Interpolate momentum to sides using appropriate limiter
            interpolateToSides(patch, interp_mom_data, *mom_un_data, u_adv_data, d_limiter);

            // Finally do flux differencing
            fluxDifference(patch, *N_un_data, interp_mom_data, u_adv_data);

            // Now do the same for the solvent
            computeAdvectionVelocity(patch, u_adv_data, *us_data);
            interpolateToSides(patch, interp_mom_data, *mom_us_data, u_adv_data, d_limiter);
            fluxDifference(patch, *N_us_data, interp_mom_data, u_adv_data);
        }
    }
}

void
MultiphaseConvectiveManager::approximateForwardEuler(const double current_time,
                                                     const double new_time,
                                                     const int un_cur_idx,
                                                     const int us_cur_idx,
                                                     const int thn_cur_idx)
{
    // We only recompute ForwardEuler if we haven't performed it yet this time step.
    if (d_is_initial_approximation_filled) return;
    approximateOperator(d_N_un_idx, d_N_us_idx, current_time, un_cur_idx, us_cur_idx, thn_cur_idx);

    // Now copy the data to N0.
    d_hier_sc_data_ops->copyData(d_N0_un_idx, d_N_un_idx);
    d_hier_sc_data_ops->copyData(d_N0_us_idx, d_N_us_idx);
}

void
MultiphaseConvectiveManager::approximateTrapezoidalRule(const double current_time,
                                                        const double new_time,
                                                        const int un_cur_idx,
                                                        const int us_cur_idx,
                                                        const int thn_cur_idx,
                                                        const int un_new_idx,
                                                        const int us_new_idx,
                                                        const int thn_new_idx)
{
    // First approximate using forward euler
    approximateForwardEuler(current_time, new_time, un_cur_idx, us_cur_idx, thn_cur_idx);

    // Now compute approximate at end time
    approximateOperator(d_N_un_idx, d_N_us_idx, new_time, un_new_idx, us_new_idx, thn_new_idx);

    // Now average N and N0
    d_hier_sc_data_ops->linearSum(d_N_un_idx, 0.5, d_N_un_idx, 0.5, d_N0_un_idx);
    d_hier_sc_data_ops->linearSum(d_N_us_idx, 0.5, d_N_us_idx, 0.5, d_N0_us_idx);
}

void
MultiphaseConvectiveManager::approximateMidpointRule(const double current_time,
                                                     const double new_time,
                                                     const int un_cur_idx,
                                                     const int us_cur_idx,
                                                     const int thn_cur_idx,
                                                     const int un_new_idx,
                                                     const int us_new_idx,
                                                     const int thn_new_idx)
{
    double half_time = 0.5 * (current_time + new_time);
    // First compute "half" values
    d_hier_sc_data_ops->linearSum(d_un_scr_idx, 0.5, un_cur_idx, 0.5, un_new_idx);
    d_hier_sc_data_ops->linearSum(d_us_scr_idx, 0.5, us_cur_idx, 0.5, us_new_idx);
    d_hier_cc_data_ops->linearSum(d_thn_scr_idx, 0.5, thn_cur_idx, 0.5, thn_new_idx);

    // Now approximate operator
    approximateOperator(d_N_un_idx, d_N_us_idx, half_time, d_un_scr_idx, d_us_scr_idx, d_thn_scr_idx);

    // Note that N0 does not have a value!
}

void
MultiphaseConvectiveManager::computeAdvectionVelocity(
    const Pointer<Patch<NDIM>>& patch,
    std::array<std::unique_ptr<FaceData<NDIM, double>>, NDIM>& u_adv_data,
    const SideData<NDIM, double>& U_data)
{
    const Box<NDIM>& box = patch->getBox();
    const Box<NDIM>& u0_box = u_adv_data[0]->getBox();
    const Box<NDIM>& u1_box = u_adv_data[1]->getBox();
    navier_stokes_interp_comps2d_(box.lower(0),
                                  box.upper(0),
                                  box.lower(1),
                                  box.upper(1),
                                  U_data.getGhostCellWidth()(0),
                                  U_data.getGhostCellWidth()(1),
                                  U_data.getPointer(0),
                                  U_data.getPointer(1),
                                  u0_box.lower(0),
                                  u0_box.upper(0),
                                  u0_box.lower(1),
                                  u0_box.upper(1),
                                  u_adv_data[0]->getGhostCellWidth()(0),
                                  u_adv_data[0]->getGhostCellWidth()(1),
                                  u_adv_data[0]->getPointer(0),
                                  u_adv_data[0]->getPointer(1),
                                  u1_box.lower(0),
                                  u1_box.upper(0),
                                  u1_box.lower(1),
                                  u1_box.upper(1),
                                  u_adv_data[1]->getGhostCellWidth()(0),
                                  u_adv_data[1]->getGhostCellWidth()(1),
                                  u_adv_data[1]->getPointer(0),
                                  u_adv_data[1]->getPointer(1));
}

void
MultiphaseConvectiveManager::findNetworkMomentum(const int dst_idx, const int thn_idx, const int u_idx)
{
    // Interpolate thn to sides.
    multiply_sc_and_thn(dst_idx, u_idx, thn_idx, d_hierarchy, true);
}

void
MultiphaseConvectiveManager::findSolventMomentum(const int dst_idx, const int thn_idx, const int u_idx)
{
    multiply_sc_and_ths(dst_idx, u_idx, thn_idx, d_hierarchy, true);
}

void
MultiphaseConvectiveManager::interpolateToSides(Pointer<Patch<NDIM>>& patch,
                                                std::array<std::unique_ptr<FaceData<NDIM, double>>, NDIM>& interp_data,
                                                const SideData<NDIM, double>& mom_data,
                                                const std::array<std::unique_ptr<FaceData<NDIM, double>>, NDIM>& u_data,
                                                const LimiterType& limiter)
{
    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM>& u0_box = u_data[0]->getBox();
    const Box<NDIM>& u1_box = u_data[1]->getBox();
    switch (limiter)
    {
    case UPWIND:
        vc_navier_stokes_upwind_quantity2d_(patch_box.lower(0),
                                            patch_box.upper(0),
                                            patch_box.lower(1),
                                            patch_box.upper(1),
                                            mom_data.getGhostCellWidth()(0),
                                            mom_data.getGhostCellWidth()(1),
                                            mom_data.getPointer(0),
                                            mom_data.getPointer(1),
                                            u0_box.lower(0),
                                            u0_box.upper(0),
                                            u0_box.lower(1),
                                            u0_box.upper(1),
                                            u_data[0]->getGhostCellWidth()(0),
                                            u_data[0]->getGhostCellWidth()(1),
                                            u_data[0]->getPointer(0),
                                            u_data[0]->getPointer(1),
                                            interp_data[0]->getGhostCellWidth()(0),
                                            interp_data[0]->getGhostCellWidth()(1),
                                            interp_data[0]->getPointer(0),
                                            interp_data[0]->getPointer(1),
                                            u1_box.lower(0),
                                            u1_box.upper(0),
                                            u1_box.lower(1),
                                            u1_box.upper(1),
                                            u_data[1]->getGhostCellWidth()(0),
                                            u_data[1]->getGhostCellWidth()(1),
                                            u_data[1]->getPointer(0),
                                            u_data[1]->getPointer(1),
                                            interp_data[1]->getGhostCellWidth()(0),
                                            interp_data[1]->getGhostCellWidth()(1),
                                            interp_data[1]->getPointer(0),
                                            interp_data[1]->getPointer(1));
        break;
    case CUI:
        vc_navier_stokes_cui_quantity2d_(patch_box.lower(0),
                                         patch_box.upper(0),
                                         patch_box.lower(1),
                                         patch_box.upper(1),
                                         mom_data.getGhostCellWidth()(0),
                                         mom_data.getGhostCellWidth()(1),
                                         mom_data.getPointer(0),
                                         mom_data.getPointer(1),
                                         u0_box.lower(0),
                                         u0_box.upper(0),
                                         u0_box.lower(1),
                                         u0_box.upper(1),
                                         u_data[0]->getGhostCellWidth()(0),
                                         u_data[0]->getGhostCellWidth()(1),
                                         u_data[0]->getPointer(0),
                                         u_data[0]->getPointer(1),
                                         interp_data[0]->getGhostCellWidth()(0),
                                         interp_data[0]->getGhostCellWidth()(1),
                                         interp_data[0]->getPointer(0),
                                         interp_data[0]->getPointer(1),
                                         u1_box.lower(0),
                                         u1_box.upper(0),
                                         u1_box.lower(1),
                                         u1_box.upper(1),
                                         u_data[1]->getGhostCellWidth()(0),
                                         u_data[1]->getGhostCellWidth()(1),
                                         u_data[1]->getPointer(0),
                                         u_data[1]->getPointer(1),
                                         interp_data[1]->getGhostCellWidth()(0),
                                         interp_data[1]->getGhostCellWidth()(1),
                                         interp_data[1]->getPointer(0),
                                         interp_data[1]->getPointer(1));
        break;
    case FBICS:
        vc_navier_stokes_fbics_quantity2d_(patch_box.lower(0),
                                           patch_box.upper(0),
                                           patch_box.lower(1),
                                           patch_box.upper(1),
                                           mom_data.getGhostCellWidth()(0),
                                           mom_data.getGhostCellWidth()(1),
                                           mom_data.getPointer(0),
                                           mom_data.getPointer(1),
                                           u0_box.lower(0),
                                           u0_box.upper(0),
                                           u0_box.lower(1),
                                           u0_box.upper(1),
                                           u_data[0]->getGhostCellWidth()(0),
                                           u_data[0]->getGhostCellWidth()(1),
                                           u_data[0]->getPointer(0),
                                           u_data[0]->getPointer(1),
                                           interp_data[0]->getGhostCellWidth()(0),
                                           interp_data[0]->getGhostCellWidth()(1),
                                           interp_data[0]->getPointer(0),
                                           interp_data[0]->getPointer(1),
                                           u1_box.lower(0),
                                           u1_box.upper(0),
                                           u1_box.lower(1),
                                           u1_box.upper(1),
                                           u_data[1]->getGhostCellWidth()(0),
                                           u_data[1]->getGhostCellWidth()(1),
                                           u_data[1]->getPointer(0),
                                           u_data[1]->getPointer(1),
                                           interp_data[1]->getGhostCellWidth()(0),
                                           interp_data[1]->getGhostCellWidth()(1),
                                           interp_data[1]->getPointer(0),
                                           interp_data[1]->getPointer(1));
        break;
    case MGAMMA:
        vc_navier_stokes_mgamma_quantity2d_(patch_box.lower(0),
                                            patch_box.upper(0),
                                            patch_box.lower(1),
                                            patch_box.upper(1),
                                            mom_data.getGhostCellWidth()(0),
                                            mom_data.getGhostCellWidth()(1),
                                            mom_data.getPointer(0),
                                            mom_data.getPointer(1),
                                            u0_box.lower(0),
                                            u0_box.upper(0),
                                            u0_box.lower(1),
                                            u0_box.upper(1),
                                            u_data[0]->getGhostCellWidth()(0),
                                            u_data[0]->getGhostCellWidth()(1),
                                            u_data[0]->getPointer(0),
                                            u_data[0]->getPointer(1),
                                            interp_data[0]->getGhostCellWidth()(0),
                                            interp_data[0]->getGhostCellWidth()(1),
                                            interp_data[0]->getPointer(0),
                                            interp_data[0]->getPointer(1),
                                            u1_box.lower(0),
                                            u1_box.upper(0),
                                            u1_box.lower(1),
                                            u1_box.upper(1),
                                            u_data[1]->getGhostCellWidth()(0),
                                            u_data[1]->getGhostCellWidth()(1),
                                            u_data[1]->getPointer(0),
                                            u_data[1]->getPointer(1),
                                            interp_data[1]->getGhostCellWidth()(0),
                                            interp_data[1]->getGhostCellWidth()(1),
                                            interp_data[1]->getPointer(0),
                                            interp_data[1]->getPointer(1));
        break;
    default:
        TBOX_ERROR("Unknown limiter " << IBAMR::enum_to_string(limiter) << "\n");
    }
}

void
MultiphaseConvectiveManager::fluxDifference(Pointer<Patch<NDIM>>& patch,
                                            SideData<NDIM, double>& N_data,
                                            const std::array<std::unique_ptr<FaceData<NDIM, double>>, NDIM>& mom_data,
                                            const std::array<std::unique_ptr<FaceData<NDIM, double>>, NDIM>& U_adv_data)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM>& shft_box = U_adv_data[axis]->getBox();
        convect_derivative2d_(dx,
                              shft_box.lower(0),
                              shft_box.upper(0),
                              shft_box.lower(1),
                              shft_box.upper(1),
                              U_adv_data[axis]->getGhostCellWidth()(0),
                              U_adv_data[axis]->getGhostCellWidth()(1),
                              mom_data[axis]->getGhostCellWidth()(0),
                              mom_data[axis]->getGhostCellWidth()(1),
                              U_adv_data[axis]->getPointer(0),
                              U_adv_data[axis]->getPointer(1),
                              mom_data[axis]->getPointer(0),
                              mom_data[axis]->getPointer(1),
                              N_data.getGhostCellWidth()(0),
                              N_data.getGhostCellWidth()(1),
                              N_data.getPointer(axis));
    }
}
} // namespace multiphase
