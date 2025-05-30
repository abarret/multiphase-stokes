#ifndef included_multiphase_fd_operators_inc
#define included_multiphase_fd_operators_inc

#include "multiphase/fd_operators.h"
#include "multiphase/utility_functions.h"

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <SideData.h>

namespace multiphase
{
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
    const double xi = params.xi;
    const double nu_n = params.nu_n;
    const double nu_s = params.nu_s;

    const double eta_n = params.eta_n;
    const double eta_s = params.eta_s;
    const double lambda_n = params.lambda_n;
    const double lambda_s = params.lambda_s;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
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

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)

        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)
        SAMRAI::pdat::SideIndex<NDIM> l_y_idx(idx_c_low, 1, 0);    // (i-1,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> u_y_idx(idx_c_low, 1, 1);    // (i-1,j+1/2)

        SAMRAI::pdat::NodeIndex<NDIM> idx_n_l(idx.toCell(1), SAMRAI::pdat::NodeIndex<NDIM>::LowerLeft);
        SAMRAI::pdat::NodeIndex<NDIM> idx_n_u(idx.toCell(1), SAMRAI::pdat::NodeIndex<NDIM>::UpperLeft);
        double thn_lower = (*thn_sc_data)(idx);
        double thn_imhalf_jphalf = (*thn_nc_data)(idx_n_u);
        double thn_imhalf_jmhalf = (*thn_nc_data)(idx_n_l);

        // components of first row (x-component of network vel) of network equation
        double ddx_Thn_dx_un = 1.0 / (dx[0] * dx[0]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddy_Thn_dy_un = 1.0 / (dx[1] * dx[1]) *
                               (thn_imhalf_jphalf * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddy_Thn_dx_vn = 1.0 / (dx[1] * dx[0]) *
                               (thn_imhalf_jphalf * ((*un_data)(upper_y_idx) - (*un_data)(u_y_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_y_idx) - (*un_data)(l_y_idx)));
        double ddx_Thn_dy_vn = 1.0 / (dx[0] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(upper_y_idx) - (*un_data)(lower_y_idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(u_y_idx) - (*un_data)(l_y_idx)));
        double visc_stress_n = (2.0 * eta_n + lambda_n) * ddx_Thn_dx_un + eta_n * ddy_Thn_dy_un +
                               eta_n * ddy_Thn_dx_vn + lambda_n * ddx_Thn_dy_vn;
        double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
        double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) = D_u * (visc_stress_n + drag_n) + D_p * (pressure_n) + C * thn_lower * (*un_data)(idx);

        // solvent equation
        double ddx_Ths_dx_us = 1.0 / (dx[0] * dx[0]) *
                               (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddy_Ths_dy_us = 1.0 / (dx[1] * dx[1]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddy_Ths_dx_vs = 1.0 / (dx[1] * dx[0]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(upper_y_idx) - (*us_data)(u_y_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_y_idx) - (*us_data)(l_y_idx)));
        double ddx_Ths_dy_vs =
            1.0 / (dx[0] * dx[1]) *
            (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_y_idx) - (*us_data)(lower_y_idx)) -
             convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_y_idx) - (*us_data)(l_y_idx)));
        double visc_stress_s = (2.0 * eta_s + lambda_s) * ddx_Ths_dx_us + eta_s * ddy_Ths_dy_us +
                               eta_s * ddy_Ths_dx_vs + lambda_s * ddx_Ths_dy_vs;
        double drag_s = -xi / nu_s * thn_lower * convertToThs(thn_lower) * ((*us_data)(idx) - (*un_data)(idx));
        double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_us_data)(idx) =
            D_u * (visc_stress_s + drag_s) + D_p * (pressure_s) + C * convertToThs(thn_lower) * (*us_data)(idx);
    }

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
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
        double ddy_Thn_dy_vn = 1.0 / (dx[1] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddx_Thn_dx_vn = 1.0 / (dx[0] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddx_Thn_dy_un = 1.0 / (dx[1] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(upper_x_idx) - (*un_data)(u_x_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_x_idx) - (*un_data)(l_x_idx)));
        double ddy_Thn_dx_un = 1.0 / (dx[0] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(upper_x_idx) - (*un_data)(lower_x_idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(u_x_idx) - (*un_data)(l_x_idx)));
        double visc_stress_n = (2.0 * eta_n + lambda_n) * ddy_Thn_dy_vn + eta_n * ddx_Thn_dy_un +
                               eta_n * ddx_Thn_dx_vn + lambda_n * ddy_Thn_dx_un;
        double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
        double pressure_n = -thn_lower / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) = D_u * (visc_stress_n + drag_n) + D_p * (pressure_n) + C * thn_lower * (*un_data)(idx);

        // Solvent equation
        double ddy_Ths_dy_vs = 1.0 / (dx[1] * dx[1]) *
                               (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddx_Ths_dx_vs = 1.0 / (dx[0] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddx_Ths_dy_us = 1.0 / (dx[1] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(upper_x_idx) - (*us_data)(u_x_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_x_idx) - (*us_data)(l_x_idx)));
        double ddy_Ths_dx_us =
            1.0 / (dx[0] * dx[1]) *
            (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_x_idx) - (*us_data)(lower_x_idx)) -
             convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_x_idx) - (*us_data)(l_x_idx)));
        double visc_stress_s = (2.0 * eta_s + lambda_s) * ddy_Ths_dy_vs + eta_s * ddx_Ths_dy_us +
                               eta_s * ddx_Ths_dx_vs + lambda_s * ddy_Ths_dx_us;
        double drag_s = -xi / nu_s * thn_lower * convertToThs(thn_lower) * ((*us_data)(idx) - (*un_data)(idx));
        double pressure_s = -convertToThs(thn_lower) / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_us_data)(idx) =
            D_u * (visc_stress_s + drag_s) + D_p * (pressure_s) + C * convertToThs(thn_lower) * (*us_data)(idx);
    }
}

inline void
accumulateMomentumForcesOnPatchConstantCoefficient(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                   const int A_un_idx,
                                                   const int A_us_idx,
                                                   const int p_idx,
                                                   const int un_idx,
                                                   const int us_idx,
                                                   const int thn_idx,
                                                   const MultiphaseParameters& params,
                                                   const double C,
                                                   const double D_u,
                                                   const double D_p)
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
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data =
        patch->getPatchData(A_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data =
        patch->getPatchData(A_us_idx); // Forces on solvent
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)

        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)
        SAMRAI::pdat::SideIndex<NDIM> l_y_idx(idx_c_low, 1, 0);    // (i-1,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> u_y_idx(idx_c_low, 1, 1);    // (i-1,j+1/2)

        // thn at sides
        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i-1/2,j)
        // thn at corners
        double thn_imhalf_jphalf = 0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + yp) +
                                           (*thn_data)(idx_c_low + yp)); // thn(i-1/2,j+1/2)
        double thn_imhalf_jmhalf = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up - yp) +
                                           (*thn_data)(idx_c_low - yp)); // thn(i-1/2,j-1/2)

        // components of first row (x-component of network vel) of network equation
        double ddx_Thn_dx_un = 1.0 / (dx[0] * dx[0]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddy_Thn_dy_un = 1.0 / (dx[1] * dx[1]) *
                               (thn_imhalf_jphalf * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddy_Thn_dx_vn = 1.0 / (dx[1] * dx[0]) *
                               (thn_imhalf_jphalf * ((*un_data)(upper_y_idx) - (*un_data)(u_y_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_y_idx) - (*un_data)(l_y_idx)));
        double ddx_Thn_dy_vn = 1.0 / (dx[0] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(upper_y_idx) - (*un_data)(lower_y_idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(u_y_idx) - (*un_data)(l_y_idx)));
        double visc_stress_n = (2.0 * eta_n + lambda_n) * ddx_Thn_dx_un + eta_n * ddy_Thn_dy_un +
                               eta_n * ddy_Thn_dx_vn + lambda_n * ddx_Thn_dy_vn;
        double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
        double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) = D_u * (visc_stress_n + drag_n) + D_p * (pressure_n) + C * thn_lower * (*un_data)(idx);

        // solvent equation
        double ddx_Ths_dx_us = 1.0 / (dx[0] * dx[0]) *
                               (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddy_Ths_dy_us = 1.0 / (dx[1] * dx[1]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddy_Ths_dx_vs = 1.0 / (dx[1] * dx[0]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(upper_y_idx) - (*us_data)(u_y_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_y_idx) - (*us_data)(l_y_idx)));
        double ddx_Ths_dy_vs =
            1.0 / (dx[0] * dx[1]) *
            (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_y_idx) - (*us_data)(lower_y_idx)) -
             convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_y_idx) - (*us_data)(l_y_idx)));
        double visc_stress_s = (2.0 * eta_s + lambda_s) * ddx_Ths_dx_us + eta_s * ddy_Ths_dy_us +
                               eta_s * ddy_Ths_dx_vs + lambda_s * ddx_Ths_dy_vs;
        double drag_s = -xi / nu_s * thn_lower * convertToThs(thn_lower) * ((*us_data)(idx) - (*un_data)(idx));
        double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_us_data)(idx) =
            D_u * (visc_stress_s + drag_s) + D_p * (pressure_s) + C * convertToThs(thn_lower) * (*us_data)(idx);
    }

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)

        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> l_x_idx(idx_c_low, 0, 0);    // (i-1/2,j-1)
        SAMRAI::pdat::SideIndex<NDIM> u_x_idx(idx_c_low, 0, 1);    // (i+1/2,j-1)

        // thn at sides
        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i,j-1/2)

        // thn at corners
        double thn_imhalf_jmhalf = 0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up - xp) +
                                           (*thn_data)(idx_c_low - xp)); // thn(i-1/2,j-1/2)
        double thn_iphalf_jmhalf = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up + xp) +
                                           (*thn_data)(idx_c_low + xp)); // thn(i+1/2,j-1/2)

        // components of second row (y-component of network vel) of network equation
        double ddy_Thn_dy_vn = 1.0 / (dx[1] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddx_Thn_dx_vn = 1.0 / (dx[0] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddx_Thn_dy_un = 1.0 / (dx[1] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(upper_x_idx) - (*un_data)(u_x_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_x_idx) - (*un_data)(l_x_idx)));
        double ddy_Thn_dx_un = 1.0 / (dx[0] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(upper_x_idx) - (*un_data)(lower_x_idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(u_x_idx) - (*un_data)(l_x_idx)));
        double visc_stress_n = (2.0 * eta_n + lambda_n) * ddy_Thn_dy_vn + eta_n * ddx_Thn_dy_un +
                               eta_n * ddx_Thn_dx_vn + lambda_n * ddy_Thn_dx_un;
        double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
        double pressure_n = -thn_lower / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) = D_u * (visc_stress_n + drag_n) + D_p * (pressure_n) + C * thn_lower * (*un_data)(idx);

        // Solvent equation
        double ddy_Ths_dy_vs = 1.0 / (dx[1] * dx[1]) *
                               (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddx_Ths_dx_vs = 1.0 / (dx[0] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddx_Ths_dy_us = 1.0 / (dx[1] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(upper_x_idx) - (*us_data)(u_x_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_x_idx) - (*us_data)(l_x_idx)));
        double ddy_Ths_dx_us =
            1.0 / (dx[0] * dx[1]) *
            (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_x_idx) - (*us_data)(lower_x_idx)) -
             convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_x_idx) - (*us_data)(l_x_idx)));
        double visc_stress_s = (2.0 * eta_s + lambda_s) * ddy_Ths_dy_vs + eta_s * ddx_Ths_dy_us +
                               eta_s * ddx_Ths_dx_vs + lambda_s * ddy_Ths_dx_us;
        double drag_s = -xi / nu_s * thn_lower * convertToThs(thn_lower) * ((*us_data)(idx) - (*un_data)(idx));
        double pressure_s = -convertToThs(thn_lower) / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_us_data)(idx) =
            D_u * (visc_stress_s + drag_s) + D_p * (pressure_s) + C * convertToThs(thn_lower) * (*us_data)(idx);
    }
}

inline void
accumulateMomentumForcesOnPatchVariableDrag(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                            const int A_un_idx,
                                            const int A_us_idx,
                                            const int p_idx,
                                            const int un_idx,
                                            const int us_idx,
                                            const int thn_idx,
                                            const MultiphaseParameters& params,
                                            const double C,
                                            const double D_u,
                                            const double D_p)
{
#ifndef NDEBUG
    TBOX_ASSERT(params.isVariableDrag());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> xi_data = patch->getPatchData(params.xi_idx);

    const double eta_n = params.eta_n;
    const double eta_s = params.eta_s;
    const double lambda_n = params.lambda_n;
    const double lambda_s = params.lambda_s;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data =
        patch->getPatchData(A_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data =
        patch->getPatchData(A_us_idx); // Forces on solvent
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)

        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)
        SAMRAI::pdat::SideIndex<NDIM> l_y_idx(idx_c_low, 1, 0);    // (i-1,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> u_y_idx(idx_c_low, 1, 1);    // (i-1,j+1/2)

        // thn at sides
        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i-1/2,j)
        // thn at corners
        double thn_imhalf_jphalf = 0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + yp) +
                                           (*thn_data)(idx_c_low + yp)); // thn(i-1/2,j+1/2)
        double thn_imhalf_jmhalf = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up - yp) +
                                           (*thn_data)(idx_c_low - yp)); // thn(i-1/2,j-1/2)

        // components of first row (x-component of network vel) of network equation
        double ddx_Thn_dx_un = 1.0 / (dx[0] * dx[0]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddy_Thn_dy_un = 1.0 / (dx[1] * dx[1]) *
                               (thn_imhalf_jphalf * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddy_Thn_dx_vn = 1.0 / (dx[1] * dx[0]) *
                               (thn_imhalf_jphalf * ((*un_data)(upper_y_idx) - (*un_data)(u_y_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_y_idx) - (*un_data)(l_y_idx)));
        double ddx_Thn_dy_vn = 1.0 / (dx[0] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(upper_y_idx) - (*un_data)(lower_y_idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(u_y_idx) - (*un_data)(l_y_idx)));

        double visc_stress_n = (2.0 * eta_n + lambda_n) * ddx_Thn_dx_un + eta_n * ddy_Thn_dy_un +
                               eta_n * ddy_Thn_dx_vn + lambda_n * ddx_Thn_dy_vn;
        double drag_n = -(*xi_data)(idx) * ((*un_data)(idx) - (*us_data)(idx));
        double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) = D_u * (visc_stress_n + drag_n) + D_p * (pressure_n) + C * thn_lower * (*un_data)(idx);

        // solvent equation
        double ddx_Ths_dx_us = 1.0 / (dx[0] * dx[0]) *
                               (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddy_Ths_dy_us = 1.0 / (dx[1] * dx[1]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddy_Ths_dx_vs = 1.0 / (dx[1] * dx[0]) *
                               (convertToThs(thn_imhalf_jphalf) * ((*us_data)(upper_y_idx) - (*us_data)(u_y_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_y_idx) - (*us_data)(l_y_idx)));
        double ddx_Ths_dy_vs =
            1.0 / (dx[0] * dx[1]) *
            (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_y_idx) - (*us_data)(lower_y_idx)) -
             convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_y_idx) - (*us_data)(l_y_idx)));
        double visc_stress_s = (2.0 * eta_s + lambda_s) * ddx_Ths_dx_us + eta_s * ddy_Ths_dy_us +
                               eta_s * ddy_Ths_dx_vs + lambda_s * ddx_Ths_dy_vs;
        double drag_s = -(*xi_data)(idx) * ((*us_data)(idx) - (*un_data)(idx));
        double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_us_data)(idx) =
            D_u * (visc_stress_s + drag_s) + D_p * (pressure_s) + C * convertToThs(thn_lower) * (*us_data)(idx);
    }

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)

        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> l_x_idx(idx_c_low, 0, 0);    // (i-1/2,j-1)
        SAMRAI::pdat::SideIndex<NDIM> u_x_idx(idx_c_low, 0, 1);    // (i+1/2,j-1)

        // thn at sides
        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i,j-1/2)

        // thn at corners
        double thn_imhalf_jmhalf = 0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up - xp) +
                                           (*thn_data)(idx_c_low - xp)); // thn(i-1/2,j-1/2)
        double thn_iphalf_jmhalf = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up + xp) +
                                           (*thn_data)(idx_c_low + xp)); // thn(i+1/2,j-1/2)

        // components of second row (y-component of network vel) of network equation
        double ddy_Thn_dy_vn = 1.0 / (dx[1] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - yp)));
        double ddx_Thn_dx_vn = 1.0 / (dx[0] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - xp)));
        double ddx_Thn_dy_un = 1.0 / (dx[1] * dx[0]) *
                               (thn_iphalf_jmhalf * ((*un_data)(upper_x_idx) - (*un_data)(u_x_idx)) -
                                thn_imhalf_jmhalf * ((*un_data)(lower_x_idx) - (*un_data)(l_x_idx)));
        double ddy_Thn_dx_un = 1.0 / (dx[0] * dx[1]) *
                               ((*thn_data)(idx_c_up) * ((*un_data)(upper_x_idx) - (*un_data)(lower_x_idx)) -
                                (*thn_data)(idx_c_low) * ((*un_data)(u_x_idx) - (*un_data)(l_x_idx)));
        double visc_stress_n = (2.0 * eta_n + lambda_n) * ddy_Thn_dy_vn + eta_n * ddx_Thn_dy_un +
                               eta_n * ddx_Thn_dx_vn + lambda_n * ddy_Thn_dx_un;
        double drag_n = -(*xi_data)(idx) * ((*un_data)(idx) - (*us_data)(idx));
        double pressure_n = -thn_lower / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) = D_u * (visc_stress_n + drag_n) + D_p * (pressure_n) + C * thn_lower * (*un_data)(idx);

        // Solvent equation
        double ddy_Ths_dy_vs = 1.0 / (dx[1] * dx[1]) *
                               (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - yp)));
        double ddx_Ths_dx_vs = 1.0 / (dx[0] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - xp)));
        double ddx_Ths_dy_us = 1.0 / (dx[1] * dx[0]) *
                               (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(upper_x_idx) - (*us_data)(u_x_idx)) -
                                convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_x_idx) - (*us_data)(l_x_idx)));
        double ddy_Ths_dx_us =
            1.0 / (dx[0] * dx[1]) *
            (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_x_idx) - (*us_data)(lower_x_idx)) -
             convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_x_idx) - (*us_data)(l_x_idx)));
        double visc_stress_s = (2.0 * eta_s + lambda_s) * ddy_Ths_dy_vs + eta_s * ddx_Ths_dy_us +
                               eta_s * ddx_Ths_dx_vs + lambda_s * ddy_Ths_dx_us;
        double drag_s = -(*xi_data)(idx) * ((*us_data)(idx) - (*un_data)(idx));
        double pressure_s = -convertToThs(thn_lower) / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_us_data)(idx) =
            D_u * (visc_stress_s + drag_s) + D_p * (pressure_s) + C * convertToThs(thn_lower) * (*us_data)(idx);
    }
}

inline void
applyCoincompressibility(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                         const int A_idx,
                         const int un_idx,
                         const int us_idx,
                         const int thn_idx,
                         const double D)
{
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> A_data = patch->getPatchData(A_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);
    for (SAMRAI::pdat::CellIterator<NDIM> ci(patch->getBox()); ci; ci++) // cell-centers
    {
        const SAMRAI::pdat::CellIndex<NDIM>& idx = ci();

        SAMRAI::pdat::SideIndex<NDIM> lower_x_idx(idx, 0, 0); // (i-1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> upper_x_idx(idx, 0, 1); // (i+1/2,j)
        SAMRAI::pdat::SideIndex<NDIM> lower_y_idx(idx, 1, 0); // (i,j-1/2)
        SAMRAI::pdat::SideIndex<NDIM> upper_y_idx(idx, 1, 1); // (i,j+1/2)

        // thn at sides
        double thn_lower_x = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx - xp)); // thn(i-1/2,j)
        double thn_upper_x = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx + xp)); // thn(i+1/2,j)
        double thn_lower_y = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx - yp)); // thn(i,j-1/2)
        double thn_upper_y = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx + yp)); // thn(i,j+1/2)

        double div_un_dot_thn_dx =
            ((thn_upper_x * (*un_data)(upper_x_idx)) - (thn_lower_x * (*un_data)(lower_x_idx))) / dx[0];
        double div_un_dot_thn_dy =
            ((thn_upper_y * (*un_data)(upper_y_idx)) - (thn_lower_y * (*un_data)(lower_y_idx))) / dx[1];
        double div_un_thn = div_un_dot_thn_dx + div_un_dot_thn_dy;

        double div_us_dot_ths_dx = ((convertToThs(thn_upper_x) * (*us_data)(upper_x_idx)) -
                                    (convertToThs(thn_lower_x) * (*us_data)(lower_x_idx))) /
                                   dx[0];
        double div_us_dot_ths_dy = ((convertToThs(thn_upper_y) * (*us_data)(upper_y_idx)) -
                                    (convertToThs(thn_lower_y) * (*us_data)(lower_y_idx))) /
                                   dx[1];
        double div_us_ths = div_us_dot_ths_dx + div_us_dot_ths_dy;
        (*A_data)(idx) = D * (div_un_thn + div_us_ths);
    }
}
} // namespace multiphase

#endif
