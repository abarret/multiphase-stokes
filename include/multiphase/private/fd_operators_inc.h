#ifndef included_multiphase_fd_operators_inc
#define included_multiphase_fd_operators_inc

#include "multiphase/fd_operators.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "Box.h"
#include "Patch.h"

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <SideData.h>

#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

// Define Fortran routines to accumulate momentum without pressure
#define m_W_P_P_C_C IBTK_FC_FUNC_(m_w_p_p_c_c, FDOPS)
#define m_W_P_P_Var_Drag IBTK_FC_FUNC_(m_w_p_p_var_drag, FDOPS)
#define co_div IBTK_FC_FUNC_(co_div, FDOPS)
#define multiphase_grad_fcn IBTK_FC_FUNC_(multiphase_grad, FDOPS)
#define multiphase_grad_accum IBTK_FC_FUNC_(multiphase_grad_accum, FDOPS)

extern "C"
{
    void m_W_P_P_C_C(const double*, // dx
                     const int&,    // ilower0
                     const int&,    // iupper0
                     const int&,    // ilower1
                     const int&,    // iupper1
                     double* const, // un_data_0
                     double* const, // un_data_1
                     const int&,    // un_gcw
                     double* const, // us_data_0
                     double* const, // us_data_0
                     const int&,    // us_gcw
                     double* const, // f_un_data_0
                     double* const, // f_un_data_1
                     const int&,    // f_un_gcw
                     double* const, // f_us_data_0
                     double* const, // f_us_data_1
                     const int&,    // f_us_gcw
                     double* const, // thn_data
                     const int&,    // thn_gcw
                     double* const, // thn_nc
                     const int&,
                     double* const, // thn_sc_0
                     double* const, // thn_sc_1,
                     const int&,
                     const double&,  // eta_n
                     const double&,  // eta_s
                     const double&,  // lambda_n
                     const double&,  // lambda_s
                     const double&,  // nu_n
                     const double&,  // nu_s
                     const double&,  // xi
                     const double&,  // C in C*u term
                     const double&); // D

    void m_W_P_P_Var_Drag(const double*, // dx
                          const int&,    // ilower0
                          const int&,    // iupper0
                          const int&,    // ilower1
                          const int&,    // iupper1
                          double* const, // un_data_0
                          double* const, // un_data_1
                          const int&,    // un_gcw
                          double* const, // us_data_0
                          double* const, // us_data_0
                          const int&,    // us_gcw
                          double* const, // f_un_data_0
                          double* const, // f_un_data_1
                          const int&,    // f_un_gcw
                          double* const, // f_us_data_0
                          double* const, // f_us_data_1
                          const int&,    // f_us_gcw
                          double* const, // thn_data
                          const int&,    // thn_gcw
                          double* const, // thn_nc
                          const int&,
                          double* const, // thn_sc_0
                          double* const, // thn_sc_1,
                          const int&,
                          const double&,  // eta_n
                          const double&,  // eta_s
                          const double&,  // lambda_n
                          const double&,  // lambda_s
                          double* const,  // xi_data_0
                          double* const,  // xi_data_1
                          const int&,     // xi_gcw
                          const double&,  // C in C*u term
                          const double&); // D
    void co_div(double* const,
                const int&,
                double* const,
                double* const,
                const int&,
                double* const,
                double* const,
                const int&,
                double* const,
                double* const,
                const int&,
                const int&,
                const int&,
                const int&,
                const int&,
                const double* const,
                const double&);
    void multiphase_grad_fcn(double* const,
                             double* const,
                             const int&,
                             double* const,
                             double* const,
                             const int&,
                             double* const,
                             double* const,
                             const int&,
                             double* const,
                             const int&,
                             const int&,
                             const int&,
                             const int&,
                             const int&,
                             const double* const,
                             const double&);
    void multiphase_grad_accum(double* const,
                               double* const,
                               const int&,
                               double* const,
                               double* const,
                               const int&,
                               double* const,
                               double* const,
                               const int&,
                               double* const,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const int&,
                               const double* const,
                               const double&);
}
namespace multiphase
{
inline void
accumulateMomentumForcesConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                            const int A_un_idx,
                                            const int A_us_idx,
                                            const int p_idx,
                                            const int un_idx,
                                            const int us_idx,
                                            const int thn_idx,
                                            const int thn_nc_idx,
                                            const int thn_sc_idx,
                                            const MultiphaseParameters& params,
                                            const double C,
                                            const double D_u,
                                            const double D_p,
                                            int coarsest_ln,
                                            int finest_ln)
{
    set_valid_level_numbers(hierarchy, coarsest_ln, finest_ln);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            accumulateMomentumForcesOnPatchConstantCoefficient(
                patch, A_un_idx, A_us_idx, p_idx, un_idx, us_idx, thn_idx, thn_nc_idx, thn_sc_idx, params, C, D_u, D_p);
        }
    }
}

inline void
accumulateMomentumForcesConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                            const int A_un_idx,
                                            const int A_us_idx,
                                            const int p_idx,
                                            const int un_idx,
                                            const int us_idx,
                                            const VolumeFractionDataManager& thn_manager,
                                            const MultiphaseParameters& params,
                                            const double C,
                                            const double D_u,
                                            const double D_p,
                                            int coarsest_ln,
                                            int finest_ln)
{
    accumulateMomentumForcesConstantCoefficient(hierarchy,
                                                A_un_idx,
                                                A_us_idx,
                                                p_idx,
                                                un_idx,
                                                us_idx,
                                                thn_manager.getCellIndex(),
                                                thn_manager.getNodeIndex(),
                                                thn_manager.getSideIndex(),
                                                params,
                                                C,
                                                D_u,
                                                D_p,
                                                coarsest_ln,
                                                finest_ln);
}

inline void
accumulateMomentumForcesOnPatchConstantCoefficient(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                   const int A_un_idx,
                                                   const int A_us_idx,
                                                   const int p_idx,
                                                   const int un_idx,
                                                   const int us_idx,
                                                   const int thn_idx,
                                                   const int thn_nc_idx,
                                                   const int thn_sc_idx,
                                                   const MultiphaseParameters& params,
                                                   const double C,
                                                   const double D_u,
                                                   const double D_p)
{
#ifndef NDEBUG
    TBOX_ASSERT(!params.isVariableDrag());
#endif

    accumulateMomentumWithoutPressureOnPatchConstantCoefficient(
        patch, A_un_idx, A_us_idx, un_idx, us_idx, thn_idx, thn_nc_idx, thn_sc_idx, params, C, D_u);

    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_sc_data = patch->getPatchData(thn_sc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data =
        patch->getPatchData(A_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data =
        patch->getPatchData(A_us_idx); // Forces on solvent
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y

    // Now add the pressure force to the momentum forces computed from above.
    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si();         // axis = 0, (i-1/2,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0); // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);  // (i,j)

        double thn_lower = (*thn_sc_data)(idx);
        double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) += D_p * (pressure_n);
        (*A_us_data)(idx) += D_p * (pressure_s);
    }

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si();         // axis = 1, (i,j-1/2)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0); // (i,j-1)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);  // (i,j)

        double thn_lower = (*thn_sc_data)(idx);
        double pressure_n = -thn_lower / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        double pressure_s = -convertToThs(thn_lower) / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) += D_p * (pressure_n);
        (*A_us_data)(idx) += D_p * (pressure_s);
    }
}

inline void
accumulateMomentumForcesVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                     const int A_un_idx,
                                     const int A_us_idx,
                                     const int p_idx,
                                     const int un_idx,
                                     const int us_idx,
                                     const int thn_idx,
                                     const int thn_nc_idx,
                                     const int thn_sc_idx,
                                     const MultiphaseParameters& params,
                                     const double C,
                                     const double D_u,
                                     const double D_p,
                                     int coarsest_ln,
                                     int finest_ln)
{
    set_valid_level_numbers(hierarchy, coarsest_ln, finest_ln);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            accumulateMomentumForcesOnPatchVariableDrag(
                patch, A_un_idx, A_us_idx, p_idx, un_idx, us_idx, thn_idx, thn_nc_idx, thn_sc_idx, params, C, D_u, D_p);
        }
    }
}

inline void
accumulateMomentumForcesVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                     const int A_un_idx,
                                     const int A_us_idx,
                                     const int p_idx,
                                     const int un_idx,
                                     const int us_idx,
                                     const VolumeFractionDataManager& thn_manager,
                                     const MultiphaseParameters& params,
                                     const double C,
                                     const double D_u,
                                     const double D_p,
                                     int coarsest_ln,
                                     int finest_ln)
{
    accumulateMomentumForcesVariableDrag(hierarchy,
                                         A_un_idx,
                                         A_us_idx,
                                         p_idx,
                                         un_idx,
                                         us_idx,
                                         thn_manager.getCellIndex(),
                                         thn_manager.getNodeIndex(),
                                         thn_manager.getSideIndex(),
                                         params,
                                         C,
                                         D_u,
                                         D_p,
                                         coarsest_ln,
                                         finest_ln);
}

inline void
accumulateMomentumForcesOnPatchVariableDrag(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                            const int A_un_idx,
                                            const int A_us_idx,
                                            const int p_idx,
                                            const int un_idx,
                                            const int us_idx,
                                            const int thn_idx,
                                            const int thn_nc_idx,
                                            const int thn_sc_idx,
                                            const MultiphaseParameters& params,
                                            const double C,
                                            const double D_u,
                                            const double D_p)
{
#ifndef NDEBUG
    TBOX_ASSERT(params.isVariableDrag());
#endif
    accumulateMomentumWithoutPressureOnPatchVariableDrag(
        patch, A_un_idx, A_us_idx, un_idx, us_idx, thn_idx, thn_nc_idx, thn_sc_idx, params, C, D_u);

    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_data = patch->getPatchData(thn_sc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data =
        patch->getPatchData(A_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data =
        patch->getPatchData(A_us_idx); // Forces on solvent
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);

    // Now add the pressure force to the momentum forces computed from above.
    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si();         // axis = 0, (i-1/2,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0); // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);  // (i,j)

        double thn_lower = (*thn_data)(idx);
        double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) += D_p * (pressure_n);
        (*A_us_data)(idx) += D_p * (pressure_s);
    }

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si();         // axis = 1, (i,j-1/2)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0); // (i,j-1)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);  // (i,j)

        double thn_lower = (*thn_data)(idx);
        double pressure_n = -thn_lower / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        double pressure_s = -convertToThs(thn_lower) / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) += D_p * (pressure_n);
        (*A_us_data)(idx) += D_p * (pressure_s);
    }
}

inline void
accumulateMomentumWithoutPressureVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                              const int F_un_idx,
                                              const int F_us_idx,
                                              const int un_idx,
                                              const int us_idx,
                                              const VolumeFractionDataManager& thn_manager,
                                              const MultiphaseParameters& params,
                                              const double C,
                                              const double D_u,
                                              int coarsest_ln,
                                              int finest_ln)
{
    accumulateMomentumWithoutPressureVariableDrag(hierarchy,
                                                  F_un_idx,
                                                  F_us_idx,
                                                  un_idx,
                                                  us_idx,
                                                  thn_manager.getCellIndex(),
                                                  thn_manager.getNodeIndex(),
                                                  thn_manager.getSideIndex(),
                                                  params,
                                                  C,
                                                  D_u,
                                                  coarsest_ln,
                                                  finest_ln);
}

inline void
accumulateMomentumWithoutPressureVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                              const int F_un_idx,
                                              const int F_us_idx,
                                              const int un_idx,
                                              const int us_idx,
                                              const int thn_idx,
                                              const int thn_nc_idx,
                                              const int thn_sc_idx,
                                              const MultiphaseParameters& params,
                                              const double C,
                                              const double D_u,
                                              int coarsest_ln,
                                              int finest_ln)
{
#ifndef NDEBUG
    TBOX_ASSERT(params.isVariableDrag());
#endif
    set_valid_level_numbers(hierarchy, coarsest_ln, finest_ln);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            accumulateMomentumWithoutPressureOnPatchVariableDrag(
                patch, F_un_idx, F_us_idx, un_idx, us_idx, thn_idx, thn_nc_idx, thn_sc_idx, params, C, D_u);
        }
    }
}

inline void
accumulateMomentumWithoutPressureOnPatchVariableDrag(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                     const int F_un_idx,
                                                     const int F_us_idx,
                                                     const int un_idx,
                                                     const int us_idx,
                                                     const int thn_idx,
                                                     const int thn_nc_idx,
                                                     const int thn_sc_idx,
                                                     const MultiphaseParameters& params,
                                                     const double C,
                                                     const double D_u)
{
#ifndef NDEBUG
    TBOX_ASSERT(params.isVariableDrag());
#endif
    const double eta_n = params.eta_n;
    const double eta_s = params.eta_s;
    const double lambda_n = params.lambda_n;
    const double lambda_s = params.lambda_s;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double>> thn_nc_data = patch->getPatchData(thn_nc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_sc_data = patch->getPatchData(thn_sc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> F_un_data =
        patch->getPatchData(F_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> F_us_data =
        patch->getPatchData(F_us_idx); // Forces on solvent
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> xi_data = patch->getPatchData(params.xi_idx);

    const Box<NDIM>& patch_box = patch->getBox();
    // patch_lower(0), patch_lower(1) are min indices in x and y-dir
    const IntVector<NDIM>& patch_lower = patch_box.lower();
    // patch_upper(0), patch_upper(1) are max indices in x and y-dir
    const IntVector<NDIM>& patch_upper = patch_box.upper();

    m_W_P_P_Var_Drag(dx,
                     patch_lower(0), // ilower0
                     patch_upper(0), // iupper0
                     patch_lower(1), // ilower1
                     patch_upper(1), // iupper1
                     un_data->getPointer(0),
                     un_data->getPointer(1),
                     un_data->getGhostCellWidth().min(),
                     us_data->getPointer(0),
                     us_data->getPointer(1),
                     us_data->getGhostCellWidth().min(),
                     F_un_data->getPointer(0),
                     F_un_data->getPointer(1),
                     F_un_data->getGhostCellWidth().min(),
                     F_us_data->getPointer(0),
                     F_us_data->getPointer(1),
                     F_us_data->getGhostCellWidth().min(),
                     thn_data->getPointer(),
                     thn_data->getGhostCellWidth().min(),
                     thn_nc_data->getPointer(),
                     thn_nc_data->getGhostCellWidth().min(),
                     thn_sc_data->getPointer(0),
                     thn_sc_data->getPointer(1),
                     thn_sc_data->getGhostCellWidth().min(),
                     eta_n,
                     eta_s,
                     lambda_n,
                     lambda_s,
                     xi_data->getPointer(0),
                     xi_data->getPointer(1),
                     xi_data->getGhostCellWidth().min(),
                     C,
                     D_u);
}

inline void
accumulateMomentumWithoutPressureConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                                     const int F_un_idx,
                                                     const int F_us_idx,
                                                     const int un_idx,
                                                     const int us_idx,
                                                     const int thn_idx,
                                                     const int thn_nc_idx,
                                                     const int thn_sc_idx,
                                                     const MultiphaseParameters& params,
                                                     const double C,
                                                     const double D_u,
                                                     int coarsest_ln,
                                                     int finest_ln)
{
    set_valid_level_numbers(hierarchy, coarsest_ln, finest_ln);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            accumulateMomentumWithoutPressureOnPatchConstantCoefficient(
                patch, F_un_idx, F_us_idx, un_idx, us_idx, thn_idx, thn_nc_idx, thn_sc_idx, params, C, D_u);
        }
    }
}

inline void
accumulateMomentumWithoutPressureConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                                     const int F_un_idx,
                                                     const int F_us_idx,
                                                     const int un_idx,
                                                     const int us_idx,
                                                     const VolumeFractionDataManager& thn_manager,
                                                     const MultiphaseParameters& params,
                                                     const double C,
                                                     const double D_u,
                                                     int coarsest_ln,
                                                     int finest_ln)
{
    accumulateMomentumWithoutPressureConstantCoefficient(hierarchy,
                                                         F_un_idx,
                                                         F_us_idx,
                                                         un_idx,
                                                         us_idx,
                                                         thn_manager.getCellIndex(),
                                                         thn_manager.getNodeIndex(),
                                                         thn_manager.getSideIndex(),
                                                         params,
                                                         C,
                                                         D_u,
                                                         coarsest_ln,
                                                         finest_ln);
}

// Accumulate the momentum forces with thn interpolated to cell nodes and cell sides.
inline void
accumulateMomentumWithoutPressureOnPatchConstantCoefficient(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                            const int A_un_idx,
                                                            const int A_us_idx,
                                                            const int un_idx,
                                                            const int us_idx,
                                                            const int thn_idx,
                                                            const int thn_nc_idx,
                                                            const int thn_sc_idx,
                                                            const MultiphaseParameters& params,
                                                            const double C,
                                                            const double D_u)
{
#ifndef NDEBUG
    TBOX_ASSERT(!params.isVariableDrag());
#endif
    const double xi = params.xi;
    const double nu_n = params.nu_n;
    const double nu_s = params.nu_s;

    const double eta_n = params.eta_n;
    const double eta_s = params.eta_s;
    const double lambda_n = params.lambda_n;
    const double lambda_s = params.lambda_s;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double>> thn_nc_data = patch->getPatchData(thn_nc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_sc_data = patch->getPatchData(thn_sc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data =
        patch->getPatchData(A_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data =
        patch->getPatchData(A_us_idx); // Forces on solvent
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y

    const Box<NDIM>& patch_box = patch->getBox();
    // patch_lower(0), patch_lower(1) are min indices in x and y-dir
    const IntVector<NDIM>& patch_lower = patch_box.lower();
    // patch_upper(0), patch_upper(1) are max indices in x and y-dir
    const IntVector<NDIM>& patch_upper = patch_box.upper();

    m_W_P_P_C_C(dx,
                patch_lower(0), // ilower0
                patch_upper(0), // iupper0
                patch_lower(1), // ilower1
                patch_upper(1), // iupper1
                un_data->getPointer(0),
                un_data->getPointer(1),
                un_data->getGhostCellWidth().min(),
                us_data->getPointer(0),
                us_data->getPointer(1),
                us_data->getGhostCellWidth().min(),
                A_un_data->getPointer(0),
                A_un_data->getPointer(1),
                A_un_data->getGhostCellWidth().min(),
                A_us_data->getPointer(0),
                A_us_data->getPointer(1),
                A_us_data->getGhostCellWidth().min(),
                thn_data->getPointer(0),
                thn_data->getGhostCellWidth().min(),
                thn_nc_data->getPointer(),
                thn_nc_data->getGhostCellWidth().min(),
                thn_sc_data->getPointer(0),
                thn_sc_data->getPointer(1),
                thn_sc_data->getGhostCellWidth().min(),
                eta_n,
                eta_s,
                lambda_n,
                lambda_s,
                nu_n,
                nu_s,
                xi,
                C,
                D_u);
}

inline void
applyCoincompressibility(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                         const int A_idx,
                         const int un_idx,
                         const int us_idx,
                         const int thn_sc_idx,
                         const double D)
{
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_data = patch->getPatchData(thn_sc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> A_data = patch->getPatchData(A_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

    co_div(A_data->getPointer(),
           A_data->getGhostCellWidth().min(),
           un_data->getPointer(0),
           un_data->getPointer(1),
           un_data->getGhostCellWidth().min(),
           us_data->getPointer(0),
           us_data->getPointer(1),
           us_data->getGhostCellWidth().min(),
           thn_data->getPointer(0),
           thn_data->getPointer(1),
           thn_data->getGhostCellWidth().min(),
           patch_box.lower(0),
           patch_box.upper(0),
           patch_box.lower(1),
           patch_box.upper(1),
           dx,
           D);
}

inline void
applyCoincompressibility(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                         int A_idx,
                         int un_idx,
                         int us_idx,
                         int thn_sc_idx,
                         double D,
                         int coarsest_ln,
                         int finest_ln)
{
    set_valid_level_numbers(hierarchy, coarsest_ln, finest_ln);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            applyCoincompressibility(patch, A_idx, un_idx, us_idx, thn_sc_idx, D);
        }
    }
}

inline void
multiphase_grad_on_hierarchy(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                             const int Gun_idx,
                             const int Gus_idx,
                             const int thn_idx,
                             const int p_idx,
                             const double C,
                             const bool do_accumulate,
                             int coarsest_ln,
                             int finest_ln)
{
    set_valid_level_numbers(hierarchy, coarsest_ln, finest_ln);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            multiphase_grad(*patch, Gun_idx, Gus_idx, thn_idx, p_idx, C, do_accumulate);
        }
    }
}

inline void
multiphase_grad(const SAMRAI::hier::Patch<NDIM>& patch,
                const int Gun_idx,
                const int Gus_idx,
                const int thn_sc_idx,
                const int p_idx,
                const double C,
                const bool do_accumulate)
{
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> Gun_data = patch.getPatchData(Gun_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> Gus_data = patch.getPatchData(Gus_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch.getPatchData(p_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_data = patch.getPatchData(thn_sc_idx);

    const Box<NDIM>& box = patch.getBox();
    if (do_accumulate)
    {
        multiphase_grad_accum(Gun_data->getPointer(0),
                              Gun_data->getPointer(1),
                              Gun_data->getGhostCellWidth().max(),
                              Gus_data->getPointer(0),
                              Gus_data->getPointer(1),
                              Gus_data->getGhostCellWidth().max(),
                              thn_data->getPointer(0),
                              thn_data->getPointer(1),
                              thn_data->getGhostCellWidth().max(),
                              p_data->getPointer(),
                              p_data->getGhostCellWidth().max(),
                              box.lower(0),
                              box.lower(1),
                              box.upper(0),
                              box.upper(1),
                              dx,
                              C);
    }
    else
    {
        multiphase_grad_fcn(Gun_data->getPointer(0),
                            Gun_data->getPointer(1),
                            Gun_data->getGhostCellWidth().max(),
                            Gus_data->getPointer(0),
                            Gus_data->getPointer(1),
                            Gus_data->getGhostCellWidth().max(),
                            thn_data->getPointer(0),
                            thn_data->getPointer(1),
                            thn_data->getGhostCellWidth().max(),
                            p_data->getPointer(),
                            p_data->getGhostCellWidth().max(),
                            box.lower(0),
                            box.lower(1),
                            box.upper(0),
                            box.upper(1),
                            dx,
                            C);
    }
}

inline void
computeVelocitySubBlockOnPatch(const SAMRAI::hier::Patch<NDIM>& patch,
                               const int A_un_idx,
                               const int A_us_idx,
                               const int un_idx,
                               const int us_idx,
                               const int thn_idx,
                               const int thn_nc_idx,
                               const int thn_sc_idx,
                               const MultiphaseParameters& params,
                               const double C,
                               const double D)
{
    const double eta_s = params.eta_s;
    const double lambda_s = params.lambda_s;
    const double eta_n = params.eta_n;
    const double lambda_n = params.lambda_n;

    const double xi = params.xi;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data = patch.getPatchData(A_un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data = patch.getPatchData(A_us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch.getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch.getPatchData(us_idx);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double>> thn_nc_data = patch.getPatchData(thn_nc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_sc_data = patch.getPatchData(thn_sc_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_cc_data = patch.getPatchData(thn_idx);

    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch.getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)

        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)
        SAMRAI::pdat::SideIndex<NDIM> l_y_idx(idx_c_low, 1, 0);    // (i-1,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> u_y_idx(idx_c_low, 1, 1);    // (i-1,j+1/2)

        SAMRAI::pdat::NodeIndex<NDIM> idx_n_l(idx.toCell(1), NodeIndex<NDIM>::LowerLeft);
        SAMRAI::pdat::NodeIndex<NDIM> idx_n_u(idx.toCell(1), NodeIndex<NDIM>::UpperLeft);
        double thn_lower = (*thn_sc_data)(idx);
        double thn_imhalf_jphalf = (*thn_nc_data)(idx_n_u);
        double thn_imhalf_jmhalf = (*thn_nc_data)(idx_n_l);

        // components of first row (x-component of network vel) of network equation
        double ddx_Thn_dx_un = (2.0 * eta_n - lambda_n) / (dx[0] * dx[0]) *
                               ((*thn_cc_data)(idx_c_up) * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                (*thn_cc_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddy_Thn_dy_un = eta_n / (dx[1] * dx[1]) *
                               (thn_imhalf_jphalf * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddy_Thn_dx_vn = eta_n / (dx[1] * dx[0]) *
                               (thn_imhalf_jphalf * ((*un_data)(upper_y_idx) - (*un_data)(u_y_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_y_idx) - (*un_data)(l_y_idx)));
        double ddx_Thn_dy_vn = -lambda_n / (dx[0] * dx[1]) *
                               ((*thn_cc_data)(idx_c_up) * ((*un_data)(upper_y_idx) - (*un_data)(lower_y_idx)) -
                                (*thn_cc_data)(idx_c_low) * ((*un_data)(u_y_idx) - (*un_data)(l_y_idx)));
        double drag_n = -xi * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
        const double A_un =
            C * thn_lower * (*un_data)(idx) + D * (ddx_Thn_dx_un + ddy_Thn_dy_un + ddy_Thn_dx_vn + ddx_Thn_dy_vn);
        const double xi_un = D * drag_n;
        (*A_un_data)(idx) = A_un + xi_un;

        // solvent equation
        double ddx_Ths_dx_us = (2.0 * eta_s - lambda_s) / (dx[0] * dx[0]) *
                               (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddy_Ths_dy_us = eta_s / (dx[1] * dx[1]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddy_Ths_dx_vs = eta_s / (dx[1] * dx[0]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(upper_y_idx) - (*us_data)(u_y_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_y_idx) - (*us_data)(l_y_idx)));
        double ddx_Ths_dy_vs =
            -lambda_s / (dx[0] * dx[1]) *
            (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(upper_y_idx) - (*us_data)(lower_y_idx)) -
             convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(u_y_idx) - (*us_data)(l_y_idx)));

        const double A_us = D * (ddx_Ths_dx_us + ddy_Ths_dy_us + ddy_Ths_dx_vs + ddx_Ths_dy_vs) +
                            C * convertToThs(thn_lower) * (*us_data)(idx);
        (*A_us_data)(idx) = A_un + A_us;
    }

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch.getBox(), 1); si; si++) // side-centers in y-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)

        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> l_x_idx(idx_c_low, 0, 0);    // (i-1/2,j-1)
        SAMRAI::pdat::SideIndex<NDIM> u_x_idx(idx_c_low, 0, 1);    // (i+1/2,j-1)

        SAMRAI::pdat::NodeIndex<NDIM> idx_n_l(idx.toCell(1), SAMRAI::pdat::NodeIndex<NDIM>::LowerLeft);
        SAMRAI::pdat::NodeIndex<NDIM> idx_n_u(idx.toCell(1), SAMRAI::pdat::NodeIndex<NDIM>::LowerRight);
        double thn_lower = (*thn_sc_data)(idx);
        double thn_iphalf_jmhalf = (*thn_nc_data)(idx_n_u);
        double thn_imhalf_jmhalf = (*thn_nc_data)(idx_n_l);

        // components of second row (y-component of network vel) of network equation
        double ddy_Thn_dy_un = (2.0 * eta_n - lambda_n) / (dx[1] * dx[1]) *
                               ((*thn_cc_data)(idx_c_up) * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                (*thn_cc_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddx_Thn_dx_un = eta_n / (dx[0] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddx_Thn_dy_vn = eta_n / (dx[1] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(upper_x_idx) - (*un_data)(u_x_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_x_idx) - (*un_data)(l_x_idx)));
        double ddy_Thn_dx_vn = -lambda_n / (dx[0] * dx[1]) *
                               ((*thn_cc_data)(idx_c_up) * ((*un_data)(upper_x_idx) - (*un_data)(lower_x_idx)) -
                                (*thn_cc_data)(idx_c_low) * ((*un_data)(u_x_idx) - (*un_data)(l_x_idx)));
        double drag_n = -xi * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));

        const double A_un =
            D * (ddy_Thn_dy_un + ddx_Thn_dx_un + ddx_Thn_dy_vn + ddy_Thn_dx_vn) + C * thn_lower * (*un_data)(idx);
        const double xi_un = D * drag_n;
        (*A_un_data)(idx) = A_un + xi_un;

        // Solvent equation
        double ddy_Ths_dy_us = (2.0 * eta_s - lambda_s) / (dx[1] * dx[1]) *
                               (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddx_Ths_dx_us = eta_s / (dx[0] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddx_Ths_dy_vs = eta_s / (dx[1] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(upper_x_idx) - (*us_data)(u_x_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_x_idx) - (*us_data)(l_x_idx)));
        double ddy_Ths_dx_vs =
            -lambda_s / (dx[0] * dx[1]) *
            (convertToThs((*thn_cc_data)(idx_c_up)) * ((*us_data)(upper_x_idx) - (*us_data)(lower_x_idx)) -
             convertToThs((*thn_cc_data)(idx_c_low)) * ((*us_data)(u_x_idx) - (*us_data)(l_x_idx)));
        const double A_us = D * (ddy_Ths_dy_us + ddx_Ths_dx_us + ddx_Ths_dy_vs + ddy_Ths_dx_vs) +
                            C * convertToThs(thn_lower) * (*us_data)(idx);

        (*A_us_data)(idx) = A_un + A_us;
    }
}

} // namespace multiphase
#endif
