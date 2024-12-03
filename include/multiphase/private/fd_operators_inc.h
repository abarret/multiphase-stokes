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
#define m_W_P_P_C_C_Thn_Nodes IBTK_FC_FUNC_(m_w_p_p_c_c_thn_nodes, FDOPS)
#define m_W_P_P_Var_Drag IBTK_FC_FUNC_(m_w_p_p_var_drag, FDOPS)

extern "C"
{
    void m_W_P_P_C_C(const double*,  // dx
                     const int&,     // ilower0
                     const int&,     // iupper0
                     const int&,     // ilower1
                     const int&,     // iupper1
                     double* const,  // un_data_0
                     double* const,  // un_data_1
                     const int&,     // un_gcw
                     double* const,  // us_data_0
                     double* const,  // us_data_0
                     const int&,     // us_gcw
                     double* const,  // f_un_data_0
                     double* const,  // f_un_data_1
                     const int&,     // f_un_gcw
                     double* const,  // f_us_data_0
                     double* const,  // f_us_data_1
                     const int&,     // f_us_gcw
                     double* const,  // thn_data
                     const int&,     // thn_gcw
                     const double&,  // eta_n
                     const double&,  // eta_s
                     const double&,  // nu_n
                     const double&,  // nu_s
                     const double&,  // xi
                     const double&,  // C in C*u term
                     const double&); // D

    void m_W_P_P_C_C_Thn_Nodes(const double*,  // dx
                               const int&,     // ilower0
                               const int&,     // iupper0
                               const int&,     // ilower1
                               const int&,     // iupper1
                               double* const,  // un_data_0
                               double* const,  // un_data_1
                               const int&,     // un_gcw
                               double* const,  // us_data_0
                               double* const,  // us_data_0
                               const int&,     // us_gcw
                               double* const,  // f_un_data_0
                               double* const,  // f_un_data_1
                               const int&,     // f_un_gcw
                               double* const,  // f_us_data_0
                               double* const,  // f_us_data_1
                               const int&,     // f_us_gcw
                               double* const,  // thn_data
                               double* const,  // thn_nc_data
                               double* const,  // thn_sc_data_0
                               double* const,  // thn_sc_data_1
                               const int&,     // thn_gcw
                               const int&,     // thn_nc_gcw
                               const int&,     // thn_sc_gcw
                               const double&,  // eta_n
                               const double&,  // eta_s
                               const double&,  // nu_n
                               const double&,  // nu_s
                               const double&,  // xi
                               const double&,  // C in C*u term
                               const double&); // D

    void m_W_P_P_Var_Drag(const double*,  // dx
                          const int&,     // ilower0
                          const int&,     // iupper0
                          const int&,     // ilower1
                          const int&,     // iupper1
                          double* const,  // un_data_0
                          double* const,  // un_data_1
                          const int&,     // un_gcw
                          double* const,  // us_data_0
                          double* const,  // us_data_0
                          const int&,     // us_gcw
                          double* const,  // f_un_data_0
                          double* const,  // f_un_data_1
                          const int&,     // f_un_gcw
                          double* const,  // f_us_data_0
                          double* const,  // f_us_data_1
                          const int&,     // f_us_gcw
                          double* const,  // thn_data
                          const int&,     // thn_gcw
                          const double&,  // eta_n
                          const double&,  // eta_s
                          const double&,  // nu_n
                          const double&,  // nu_s
                          double* const,  // xi_data_0
                          double* const,  // xi_data_1
                          const int&,     // xi_gcw
                          const double&,  // C in C*u term
                          const double&); // D
}
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

    accumulateMomentumWithoutPressureOnPatchConstantCoefficient(patch, 
                                                                A_un_idx, A_us_idx, un_idx, us_idx, 
                                                                thn_idx, thn_nc_idx, thn_sc_idx,
                                                                params, C, D_u);

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
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        
        double thn_lower = (*thn_sc_data)(idx);
        double pressure_n = -thn_lower / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        double pressure_s = -convertToThs(thn_lower) / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) += D_p * (pressure_n);
        (*A_us_data)(idx) += D_p * (pressure_s);
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

    accumulateMomentumWithoutPressureOnPatchConstantCoefficient(patch, A_un_idx, A_us_idx, un_idx, us_idx, thn_idx, params, C, D_u);
    
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data =
        patch->getPatchData(A_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data =
        patch->getPatchData(A_us_idx); // Forces on solvent
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);
    
    // Now add the pressure force to the momentum forces computed from above.
    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        
        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i-1/2,j)
        double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) += D_p * (pressure_n); 
        (*A_us_data)(idx) += D_p * (pressure_s);
    }

    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
        
        // components of second row (y-component of network vel) of network equation
        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i,j-1/2)
        double pressure_n = -thn_lower / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        double pressure_s = -convertToThs(thn_lower) / dx[1] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
        (*A_un_data)(idx) += D_p * (pressure_n);
        (*A_us_data)(idx) += D_p * (pressure_s);
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
    accumulateMomentumWithoutPressureOnPatchVariableDrag(patch, A_un_idx, A_us_idx, un_idx, us_idx, thn_idx, params, C, D_u);


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

    // Now add the pressure force to the momentum forces computed from above.
    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
    {
        const SAMRAI::pdat::SideIndex<NDIM>& idx = si();         // axis = 0, (i-1/2,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0); // (i-1,j)
        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);  // (i,j)

        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i-1/2,j)
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

        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i,j-1/2)
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
                                              const int thn_idx,
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
                patch, F_un_idx, F_us_idx, un_idx, us_idx, thn_idx, params, C, D_u);
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
                                                     const MultiphaseParameters& params,
                                                     const double C,
                                                     const double D_u)
{
#ifndef NDEBUG
    TBOX_ASSERT(params.isVariableDrag());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> xi_data = patch->getPatchData(params.xi_idx);

    const double eta_n = params.eta_n;
    const double eta_s = params.eta_s;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> F_un_data =
        patch->getPatchData(F_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> F_us_data =
        patch->getPatchData(F_us_idx); // Forces on solvent
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);

    double* const un_data_0 = un_data->getPointer(0);
    double* const un_data_1 = un_data->getPointer(1);
    double* const us_data_0 = us_data->getPointer(0);
    double* const us_data_1 = us_data->getPointer(1);
    double* const thn_ptr_data = thn_data->getPointer(0);
    double* const f_un_data_0 = F_un_data->getPointer(0);
    double* const f_un_data_1 = F_un_data->getPointer(1);
    double* const f_us_data_0 = F_us_data->getPointer(0);
    double* const f_us_data_1 = F_us_data->getPointer(1);

    double* const xi_data_0 = xi_data->getPointer(0);
    double* const xi_data_1 = xi_data->getPointer(1);

    const Box<NDIM>& patch_box = patch->getBox();
    const IntVector<NDIM>& patch_lower =
        patch_box.lower(); // patch_lower(0), patch_lower(1) are min indices in x and y-dir
    const IntVector<NDIM>& patch_upper =
        patch_box.upper(); // patch_upper(0), patch_upper(1) are max indices in x and y-dir

    const IntVector<NDIM>& thn_gcw = thn_data->getGhostCellWidth();
    const IntVector<NDIM>& un_gcw = un_data->getGhostCellWidth();
    const IntVector<NDIM>& us_gcw = us_data->getGhostCellWidth();
    const IntVector<NDIM>& f_un_gcw = F_un_data->getGhostCellWidth();
    const IntVector<NDIM>& f_us_gcw = F_us_data->getGhostCellWidth();
    const IntVector<NDIM>& xi_gcw = xi_data->getGhostCellWidth();

    m_W_P_P_Var_Drag(dx,
                     patch_lower(0), // ilower0
                     patch_upper(0), // iupper0
                     patch_lower(1), // ilower1
                     patch_upper(1), // iupper1
                     un_data_0,
                     un_data_1,
                     un_gcw.min(),
                     us_data_0,
                     us_data_1,
                     us_gcw.min(),
                     f_un_data_0,
                     f_un_data_1,
                     f_un_gcw.min(),
                     f_us_data_0,
                     f_us_data_1,
                     f_us_gcw.min(),
                     thn_ptr_data,
                     thn_gcw.min(),
                     params.eta_n,
                     params.eta_s,
                     params.nu_n,
                     params.nu_s,
                     xi_data_0,
                     xi_data_1,
                     xi_gcw.min(),
                     C,
                     D_u);
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

    double* const un_data_0 = un_data->getPointer(0);
    double* const un_data_1 = un_data->getPointer(1);
    double* const us_data_0 = us_data->getPointer(0);
    double* const us_data_1 = us_data->getPointer(1);
    double* const thn_ptr_data = thn_data->getPointer(0);
    double* const thn_nc_ptr = thn_nc_data->getPointer(0);
    double* const thn_sc_ptr_0 = thn_sc_data->getPointer(0);
    double* const thn_sc_ptr_1 = thn_sc_data->getPointer(1);
    double* const f_un_data_0 = A_un_data->getPointer(0);
    double* const f_un_data_1 = A_un_data->getPointer(1);
    double* const f_us_data_0 = A_us_data->getPointer(0);
    double* const f_us_data_1 = A_us_data->getPointer(1);

    const Box<NDIM>& patch_box = patch->getBox();
    const IntVector<NDIM>& patch_lower =
        patch_box.lower(); // patch_lower(0), patch_lower(1) are min indices in x and y-dir
    const IntVector<NDIM>& patch_upper =
        patch_box.upper(); // patch_upper(0), patch_upper(1) are max indices in x and y-dir

    const IntVector<NDIM>& thn_gcw = thn_data->getGhostCellWidth();
    const IntVector<NDIM>& thn_nc_gcw = thn_nc_data->getGhostCellWidth();
    const IntVector<NDIM>& thn_sc_gcw = thn_sc_data->getGhostCellWidth();
    const IntVector<NDIM>& un_gcw = un_data->getGhostCellWidth();
    const IntVector<NDIM>& us_gcw = us_data->getGhostCellWidth();
    const IntVector<NDIM>& f_un_gcw = A_un_data->getGhostCellWidth();
    const IntVector<NDIM>& f_us_gcw = A_us_data->getGhostCellWidth();

    m_W_P_P_C_C_Thn_Nodes(dx,
                          patch_lower(0), // ilower0
                          patch_upper(0), // iupper0
                          patch_lower(1), // ilower1
                          patch_upper(1), // iupper1
                          un_data_0,
                          un_data_1,
                          un_gcw.min(),
                          us_data_0,
                          us_data_1,
                          us_gcw.min(),
                          f_un_data_0,
                          f_un_data_1,
                          f_un_gcw.min(),
                          f_us_data_0,
                          f_us_data_1,
                          f_us_gcw.min(),
                          thn_ptr_data,
                          thn_nc_ptr,
                          thn_sc_ptr_0,
                          thn_sc_ptr_1,
                          thn_gcw.min(), // TODO: assert that they're same for both dir
                          thn_nc_gcw.min(),
                          thn_sc_gcw.min(),
                          params.eta_n,
                          params.eta_s,
                          params.nu_n,
                          params.nu_s,
                          params.xi,
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
                patch, F_un_idx, F_us_idx, un_idx, us_idx, thn_idx, params, C, D_u);
        }
    }
}

// Accumulate the momentum forces with thn only provided at cell centers.
inline void
accumulateMomentumWithoutPressureOnPatchConstantCoefficient(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                            const int A_un_idx,
                                                            const int A_us_idx,
                                                            const int un_idx,
                                                            const int us_idx,
                                                            const int thn_idx,
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
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_un_data =
        patch->getPatchData(A_un_idx); // Forces on network
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> A_us_data =
        patch->getPatchData(A_us_idx); // Forces on solvent
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y

    double* const un_data_0 = un_data->getPointer(0);
    double* const un_data_1 = un_data->getPointer(1);
    double* const us_data_0 = us_data->getPointer(0);
    double* const us_data_1 = us_data->getPointer(1);
    double* const thn_ptr_data = thn_data->getPointer(0);
    double* const f_un_data_0 = A_un_data->getPointer(0);
    double* const f_un_data_1 = A_un_data->getPointer(1);
    double* const f_us_data_0 = A_us_data->getPointer(0);
    double* const f_us_data_1 = A_us_data->getPointer(1);

    const Box<NDIM>& patch_box = patch->getBox();
    const IntVector<NDIM>& patch_lower =
        patch_box.lower(); // patch_lower(0), patch_lower(1) are min indices in x and y-dir
    const IntVector<NDIM>& patch_upper =
        patch_box.upper(); // patch_upper(0), patch_upper(1) are max indices in x and y-dir

    const IntVector<NDIM>& thn_gcw = thn_data->getGhostCellWidth();
    const IntVector<NDIM>& un_gcw = un_data->getGhostCellWidth();
    const IntVector<NDIM>& us_gcw = us_data->getGhostCellWidth();
    const IntVector<NDIM>& f_un_gcw = A_un_data->getGhostCellWidth();
    const IntVector<NDIM>& f_us_gcw = A_us_data->getGhostCellWidth();

    m_W_P_P_C_C(dx,
                patch_lower(0), // ilower0
                patch_upper(0), // iupper0
                patch_lower(1), // ilower1
                patch_upper(1), // iupper1
                un_data_0,
                un_data_1,
                un_gcw.min(),
                us_data_0,
                us_data_1,
                us_gcw.min(),
                f_un_data_0,
                f_un_data_1,
                f_un_gcw.min(),
                f_us_data_0,
                f_us_data_1,
                f_us_gcw.min(),
                thn_ptr_data,
                thn_gcw.min(),
                params.eta_n,
                params.eta_s,
                params.nu_n,
                params.nu_s,
                params.xi,
                C,
                D_u);

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

inline void
applyCoincompressibility(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                         int A_idx,
                         int un_idx,
                         int us_idx,
                         int thn_idx,
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
            applyCoincompressibility(patch, A_idx, un_idx, us_idx, thn_idx, D);
        }
    }
}

// TODO: This routine will be removed as we'll use IABMR's CCLaplaceOperator class to create this operator.
inline void
preconditonerBlockGTGOperator(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                              const int GtG_idx,
                              const int u_idx, // size of pressure vector
                              const int thn_idx,
                              const MultiphaseParameters& params,
                              const double C,
                              const double D_div,
                              const double D_p)
{
    // This sets up the G^T*G operator which acts on a vector u.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> GtG_data = patch->getPatchData(GtG_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> u_data = patch->getPatchData(u_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
    const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
    SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);

    for (SAMRAI::pdat::CellIterator<NDIM> ci(patch->getBox()); ci; ci++) // cell-centers
    {
        const SAMRAI::pdat::CellIndex<NDIM>& idx = ci();

        double uxx = ((*u_data)(idx-xp) - 2*(*u_data)(idx) + (*u_data)(idx+xp))/(dx[0]*dx[0]);
        double uyy = ((*u_data)(idx-yp) - 2*(*u_data)(idx) + (*u_data)(idx+yp))/(dx[0]*dx[1]);

        double ux = ((*u_data)(idx+xp) -(*u_data)(idx-xp))/dx[0];
        double uy = ((*u_data)(idx+yp) -(*u_data)(idx-yp))/dx[1];

        double thn_sq_ths_sq = (*thn_data)(idx)*(*thn_data)(idx) + convertToThs((*thn_data)(idx))*convertToThs((*thn_data)(idx));
        double th_sq_uxx = thn_sq_ths_sq * uxx;
        double th_sq_uyy = thn_sq_ths_sq * uyy;

        double ddx_th_sq = (((*thn_data)(idx+xp)*(*thn_data)(idx+xp) + convertToThs((*thn_data)(idx+xp))*convertToThs((*thn_data)(idx+xp)))
                            - ((*thn_data)(idx-xp)*(*thn_data)(idx-xp) + convertToThs((*thn_data)(idx-xp))*convertToThs((*thn_data)(idx-xp))))/dx[0];

        double ddy_th_sq = (((*thn_data)(idx+yp)*(*thn_data)(idx+yp) + convertToThs((*thn_data)(idx+yp))*convertToThs((*thn_data)(idx+yp)))
                            - ((*thn_data)(idx-yp)*(*thn_data)(idx-yp) + convertToThs((*thn_data)(idx-yp))*convertToThs((*thn_data)(idx-yp))))/dx[1];

        (*GtG_data)(idx) = D_div*D_p*(th_sq_uxx + th_sq_uyy + ddx_th_sq*ux + ddy_th_sq*uy);
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
                const int thn_idx,
                const int p_idx,
                const double C,
                const bool do_accumulate)
{
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> Gun_data = patch.getPatchData(Gun_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> Gus_data = patch.getPatchData(Gus_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> p_data = patch.getPatchData(p_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch.getPatchData(thn_idx);

    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SAMRAI::pdat::SideIterator<NDIM> si(patch.getBox(), axis); si; si++)
        {
            const SAMRAI::pdat::SideIndex<NDIM>& idx = si();
            // Note: This does not account for any synchronization that needs to occur at coarse-fine
            // interfaces. Consider using HierarchyMathOps::grad() instead.
            const double dp = ((*p_data)(idx.toCell(1)) - (*p_data)(idx.toCell(0))) / dx[axis];
            const double thn = 0.5 * ((*thn_data)(idx.toCell(1)) + (*thn_data)(idx.toCell(0)));
            if (do_accumulate)
            {
                (*Gun_data)(idx) += C * thn * dp;
                (*Gus_data)(idx) += C * convertToThs(thn) * dp;
            }
            else
            {
                (*Gun_data)(idx) = C * thn * dp;
                (*Gus_data)(idx) = C * convertToThs(thn) * dp;
            }
        }
    }
}
} // namespace multiphase
#endif
