
// Local includes
#include "multiphase/INSVCTwoFluidConvectiveManager.h"
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
INSVCTwoFluidConvectiveManager::INSVCTwoFluidConvectiveManager(std::string object_name,
                                                               Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                               Pointer<Database> input_db)
    : d_object_name(std::move(object_name)),
      d_hierarchy(hierarchy),
      d_mom_var(new SideVariable<NDIM, double>(d_object_name + "::Mom_var")),
      d_N_var(new SideVariable<NDIM, double>(d_object_name + "::N_var")),
      d_N0_var(new SideVariable<NDIM, double>(d_object_name + "::N0_var")),
      d_hier_sc_data_ops(hierarchy)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> network_ctx = var_db->getContext(d_object_name + "::NetworkCTX");
    Pointer<VariableContext> solvent_ctx = var_db->getContext(d_object_name + "::SolventCTX");
    static const int one_gcw = 3;
    d_mom_un_idx = var_db->registerVariableAndContext(d_mom_var, network_ctx, one_gcw);
    d_mom_us_idx = var_db->registerVariableAndContext(d_mom_var, solvent_ctx, one_gcw);
    d_N_un_idx = var_db->registerVariableAndContext(d_N_var, network_ctx);
    d_N_us_idx = var_db->registerVariableAndContext(d_N_var, solvent_ctx);
    d_N0_un_idx = var_db->registerVariableAndContext(d_N0_var, network_ctx);
    d_N0_us_idx = var_db->registerVariableAndContext(d_N0_var, solvent_ctx);

    d_limiter = IBAMR::string_to_enum<IBAMR::LimiterType>(input_db->getString("limiter_type"));
}

INSVCTwoFluidConvectiveManager::~INSVCTwoFluidConvectiveManager()
{
    if (getIsAllocated()) deallocateData();
    // Remove patch indices from variable database
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    std::array<int, 6> idxs{ d_mom_un_idx, d_mom_us_idx, d_N_un_idx, d_N_us_idx, d_N0_un_idx, d_N0_us_idx };
    for (const auto& idx : idxs) var_db->removePatchDataIndex(idx);
}

void
INSVCTwoFluidConvectiveManager::deallocateData()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    deallocate_patch_data({ d_mom_un_idx, d_mom_us_idx, d_N_un_idx, d_N_us_idx, d_N0_un_idx, d_N0_us_idx },
                          d_hierarchy,
                          coarsest_ln,
                          finest_ln);

    d_thn_ghost_fill.deallocateOperatorState();
    d_mom_ghost_fill.deallocateOperatorState();

    d_is_allocated = false;
}

void
INSVCTwoFluidConvectiveManager::allocateData(const double time, const int thn_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    if (!getIsAllocated())
    {
        allocate_patch_data({ d_mom_un_idx, d_mom_us_idx, d_N_un_idx, d_N_us_idx, d_N0_un_idx, d_N0_us_idx },
                            d_hierarchy,
                            time,
                            coarsest_ln,
                            finest_ln);

        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> thn_ghost_comp = { ITC(thn_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE") };
        d_thn_ghost_fill.initializeOperatorState(thn_ghost_comp, d_hierarchy, coarsest_ln, finest_ln);

        std::vector<ITC> mom_ghost_comp = {
            ITC(d_mom_un_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN"),
            ITC(d_mom_us_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN")
        };
        d_mom_ghost_fill.initializeOperatorState(mom_ghost_comp, d_hierarchy, coarsest_ln, finest_ln);
    }
}

void
INSVCTwoFluidConvectiveManager::approximateConvectiveOperator(const int dst_un_idx,
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
INSVCTwoFluidConvectiveManager::approximateConvectiveOperator(IBAMR::TimeSteppingType ts_type,
                                                              const double current_time,
                                                              const double new_time,
                                                              const int un_cur_idx,
                                                              const int us_cur_idx,
                                                              const int thn_cur_idx,
                                                              const int un_new_idx,
                                                              const int us_new_idx,
                                                              const int thn_new_idx)
{
    allocateData(current_time, thn_cur_idx);
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
        TBOX_ERROR("Unknown time stepping type " << IBAMR::enum_to_string(ts_type) << "\n");
    }
}

void
INSVCTwoFluidConvectiveManager::fillWithConvectiveOperator(const int dst_un_idx, const int dst_us_idx)
{
    if (dst_un_idx != IBTK::invalid_index) d_hier_sc_data_ops.copyData(dst_un_idx, d_N_un_idx);
    if (dst_us_idx != IBTK::invalid_index) d_hier_sc_data_ops.copyData(dst_us_idx, d_N_us_idx);
}

std::pair<int, int>
INSVCTwoFluidConvectiveManager::getConvectiveIndices() const
{
    return std::make_pair(d_N_un_idx, d_N_us_idx);
    ;
}

bool
INSVCTwoFluidConvectiveManager::getIsAllocated() const
{
    return d_is_allocated;
}

void
INSVCTwoFluidConvectiveManager::approximateOperator(const int dst_un_idx,
                                                    const int dst_us_idx,
                                                    const double eval_time,
                                                    const int un_idx,
                                                    const int us_idx,
                                                    const int thn_idx)
{
    // Fill in ghost cells for thn. Because thn may have changed types, we need to reinitialize the ghost filling
    // routines.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> thn_ghost_fill_itc = { ITC(thn_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE") };
    d_thn_ghost_fill.resetTransactionComponents(thn_ghost_fill_itc);
    d_thn_ghost_fill.fillData(eval_time);

    std::vector<ITC> u_ghost_fill_itc = { ITC(us_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE"),
                                          ITC(un_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE") };
    HierarchyGhostCellInterpolation ghost_fill;
    ghost_fill.initializeOperatorState(u_ghost_fill_itc, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    ghost_fill.fillData(eval_time);

    // Fill in N0 approximations and N approximations.
    // First find the respective momentums.
    findNetworkMomentum(d_mom_un_idx, thn_idx, un_idx);
    findSolventMomentum(d_mom_us_idx, thn_idx, us_idx);

    // Fill in ghost cells for momentum. Because the indices filled in here are owned by this class, we do not need to
    // reinitialize the ghost filling routines.
    d_mom_ghost_fill.fillData(eval_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
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
INSVCTwoFluidConvectiveManager::approximateForwardEuler(const double current_time,
                                                        const double new_time,
                                                        const int un_cur_idx,
                                                        const int us_cur_idx,
                                                        const int thn_cur_idx)
{
    approximateOperator(d_N_un_idx, d_N_us_idx, current_time, un_cur_idx, us_cur_idx, thn_cur_idx);

    // Now copy the data to N0.
    d_hier_sc_data_ops.copyData(d_N0_un_idx, d_N_un_idx);
    d_hier_sc_data_ops.copyData(d_N0_us_idx, d_N_us_idx);
}

void
INSVCTwoFluidConvectiveManager::computeAdvectionVelocity(
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
INSVCTwoFluidConvectiveManager::findNetworkMomentum(const int dst_idx, const int thn_idx, const int u_idx)
{
    // Interpolate thn to sides.
    multiply_sc_and_thn(dst_idx, u_idx, thn_idx, d_hierarchy);
}

void
INSVCTwoFluidConvectiveManager::findSolventMomentum(const int dst_idx, const int thn_idx, const int u_idx)
{
    multiply_sc_and_ths(dst_idx, u_idx, thn_idx, d_hierarchy);
}

void
INSVCTwoFluidConvectiveManager::interpolateToSides(
    Pointer<Patch<NDIM>>& patch,
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
        //    case PPM:
        //        for (unsigned int axis = 0; axis < NDIM; ++axis)
        //        {
        //            Pointer<SideData<NDIM, double> > dQ_data =
        //                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
        //            Pointer<SideData<NDIM, double> > Q_L_data =
        //                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
        //            Pointer<SideData<NDIM, double> > Q_R_data =
        //                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
        //            Pointer<SideData<NDIM, double> > Q_scratch1_data =
        //                new SideData<NDIM, double>(Q_data->getBox(), Q_data->getDepth(), Q_data->getGhostCellWidth());
        //            GODUNOV_EXTRAPOLATE_FC(side_boxes[axis].lower(0),
        //                                   side_boxes[axis].upper(0),
        //                                   side_boxes[axis].lower(1),
        //                                   side_boxes[axis].upper(1),
        //                                   Q_data->getGhostCellWidth()(0),
        //                                   Q_data->getGhostCellWidth()(1),
        //                                   Q_data->getPointer(axis),
        //                                   Q_scratch1_data->getPointer(axis),
        //                                   dQ_data->getPointer(axis),
        //                                   Q_L_data->getPointer(axis),
        //                                   Q_R_data->getPointer(axis),
        //                                   U_adv_data[axis]->getGhostCellWidth()(0),
        //                                   U_adv_data[axis]->getGhostCellWidth()(1),
        //                                   Q_half_data[axis]->getGhostCellWidth()(0),
        //                                   Q_half_data[axis]->getGhostCellWidth()(1),
        //                                   U_adv_data[axis]->getPointer(0),
        //                                   U_adv_data[axis]->getPointer(1),
        //                                   Q_half_data[axis]->getPointer(0),
        //                                   Q_half_data[axis]->getPointer(1));
        //        }
        //        break;
    default:
        TBOX_ERROR("Unknown limiter " << IBAMR::enum_to_string(limiter) << "\n");
    }
}

void
INSVCTwoFluidConvectiveManager::fluxDifference(
    Pointer<Patch<NDIM>>& patch,
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
