// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/VCTwoFluidTwoFluidStaggeredStokesOperator.h"
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
#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
} // namespace

double convertToThs(double Thn);
/////////////////////////////// PUBLIC ///////////////////////////////////////
VCTwoFluidStaggeredStokesOperator::VCTwoFluidStaggeredStokesOperator(const std::string& object_name,
                                                                     bool homogeneous_bc)
    : LinearOperator(object_name, homogeneous_bc),
      d_un_problem_coefs(d_object_name + "::un_problem_coefs"),
      d_us_problem_coefs(d_object_name + "::us_problem_coefs"),
      d_default_un_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_un_bc_coef", Pointer<Database>(nullptr))),
      d_default_us_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_us_bc_coef", Pointer<Database>(nullptr))),
      d_un_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef)),
      d_us_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef)),
      d_default_P_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(nullptr))),
      d_P_bc_coef(d_default_P_bc_coef)
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

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef),
                       std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef) d_default_P_bc_coef);

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::TwoFluidStaggeredStokesOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::TwoFluidStaggeredStokesOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::TwoFluidStaggeredStokesOperator::deallocateOperatorState()"););
    return;
} // TwoFluidStaggeredStokesOperator

VCTwoFluidStaggeredStokesOperator::~VCTwoFluidStaggeredStokesOperator()
{
    deallocateOperatorState();
    delete d_default_un_bc_coef;
    delete d_default_us_bc_coef;
    d_default_un_bc_coef = nullptr;
    d_default_us_bc_coef = nullptr;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = nullptr;
    return;
} // ~TwoFluidStaggeredStokesOperator

// Probably need another one for second fluid equation
void
VCTwoFluidStaggeredStokesOperator::setVelocityPoissonSpecifications(const PoissonSpecifications& un_problem_coefs,
                                                                    const PoissonSpecifications& us_problem_coefs)
{
    d_un_problem_coefs = un_problem_coefs;
    d_us_problem_coefs = us_problem_coefs;
    return;
} // setVelocityPoissonSpecifications

const PoissonSpecifications&
VCTwoFluidStaggeredStokesOperator::getNetworkVelocityPoissonSpecifications() const
{
    return d_un_problem_coefs;
} // getNetworkVelocityPoissonSpecifications

const PoissonSpecifications&
VCTwoFluidStaggeredStokesOperator::getSolventVelocityPoissonSpecifications() const
{
    return d_us_problem_coefs;
} // getSolventVelocityPoissonSpecifications

void
VCTwoFluidStaggeredStokesOperator::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
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
VCTwoFluidStaggeredStokesOperator::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_bc_helper = bc_helper;
    return;
} // setPhysicalBoundaryHelper

// create another member function to set-up Thn
// Thn is defined in the input file and read in using muParserCartGridFunction
void
VCTwoFluidStaggeredStokesOperator::setThnIdx(int thn_idx)
{
    d_thn_idx = thn_idx;
}

void
VCTwoFluidStaggeredStokesOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBAMR_TIMER_START(t_apply);

    // Get the vector components. These pull out patch data indices
    const int un_idx = x.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = x.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int P_idx = x.getComponentDescriptorIndex(2);  // pressure
    const int A_un_idx = y.getComponentDescriptorIndex(0);
    const int A_us_idx = y.getComponentDescriptorIndex(1);
    const int A_P_idx = y.getComponentDescriptorIndex(2);
    const int un_scratch_idx = d_x->getComponentDescriptorIndex(0);
    const int us_scratch_idx = d_x->getComponentDescriptorIndex(1);

    Pointer<SideVariable<NDIM, double>> un_sc_var = x.getComponentVariable(0);
    Pointer<SideVariable<NDIM, double>> us_sc_var = x.getComponentVariable(1);
    Pointer<CellVariable<NDIM, double>> P_cc_var = x.getComponentVariable(2);
    Pointer<SideVariable<NDIM, double>> A_un_sc_var = y.getComponentVariable(0);
    Pointer<SideVariable<NDIM, double>> A_us_sc_var = y.getComponentVariable(1);
    Pointer<CellVariable<NDIM, double>> A_P_cc_var = y.getComponentVariable(2);

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(4);
    transaction_comps[0] = InterpolationTransactionComponent(un_scratch_idx,
                                                             un_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_un_bc_coefs, // modifiy?
                                                             d_un_fill_pattern);
    transaction_comps[1] = InterpolationTransactionComponent(us_scratch_idx,
                                                             us_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_us_bc_coefs, // modify?
                                                             d_us_fill_pattern);
    transaction_comps[2] = InterpolationTransactionComponent(P_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef,
                                                             d_P_fill_pattern);
    transaction_comps[3] = InterpolationTransactionComponent(thn_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef,
                                                             d_P_fill_pattern);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_un_bc_coefs, d_P_bc_coef, un_scratch_idx, P_idx, d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time); // Fills in all of the ghost cells
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_un_bc_coefs, d_P_bc_coef);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);
    const double eta_n = 1.0;
    const double eta_s = 1.0;
    const double xi = 1.0;
    const double nu_n = 1.0;

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
            const double* const xlow =
                pgeom->getXLower(); // {xlow[0], xlow[1]} -> physical location of bottom left of box.
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            Pointer<CellData<NDIM, double>> A_P_data =
                patch->getPatchData(A_P_idx); // result of applying operator (eqn 3)
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> A_un_data =
                patch->getPatchData(A_un_idx); // result of applying operator (eqn 1)
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<SideData<NDIM, double>> A_us_data =
                patch->getPatchData(A_us_idx); // result of applying operator (eqn 2)
            IntVector<NDIM> xp(1, 0), yp(0, 1);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++) // cell-centers
            {
                const CellIndex<NDIM>& idx = ci();
                VectorNd x; // <- Eigen3 vector
                for (int d = 0; d < NDIM; ++d)
                    x[d] = xlow[d] + dx[d] * (idx(d) - idx_low(d) + 0.5); // Get's physical location of idx.

                SideIndex<NDIM> lower_x_idx(idx, 0, 0); // (i-1/2,j)
                SideIndex<NDIM> lower_y_idx(idx, 1, 0); // (i,j-1/2)
                SideIndex<NDIM> upper_x_idx(idx, 0, 1); // (i+1/2,j)
                SideIndex<NDIM> upper_y_idx(idx, 1, 1); // (i,j+1/2)

                // thn at sidess
                double thn_lower_x = 0.5 * (*thn_data)(idx) + (*thn_data)(idx - xp); // thn(i-1/2,j)
                double thn_upper_x = 0.5 * (*thn_data)(idx) + (*thn_data)(idx + xp); // thn(i+1/2,j)
                double thn_lower_y = 0.5 * (*thn_data)(idx) + (*thn_data)(idx - yp); // thn(i,j-1/2)
                double thn_upper_y = 0.5 * (*thn_data)(idx) + (*thn_data)(idx + yp); // thn(i,j+1/2)

                // conservation of mass
                double div_un_dot_thn_dx = (thn_upper_x * (*un_data)(upper_x_idx) - thn_lower_x * ( *un_data)(lower_x_idx))) / dx[0];
                double div_un_dot_thn_dy = (thn_upper_y * (*un_data)(upper_y_idx) - thn_lower_y * ( *un_data)(lower_y_idx))) / dx[1];
                double div_un_thn = div_un_dot_thn_dx + div_un_dot_thn_dy;

                double div_us_dot_ths_dx = (convertToThs(thn_upper_x) * (*us_data)(upper_x_idx)-convertToThs(thn_lower_x) * (*us_data)(lower_x_idx))) /
                                           dx[0];
                double div_us_dot_ths_dy = (convertToThs(thn_upper_y) * (*us_data)(upper_y_idx)-convertToThs(thn_lower_y) * (*us_data)(lower_y_idx))) /
                                           dx[1];
                double div_us_ths = div_us_dot_ths_dx + div_us_dot_ths_dy;
                double div_u_theta = div_un_thn + div_us_ths;
                (*A_P_data)(idx) = div_u_theta;
            }

            for (SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
            {
                const SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)

                CellIndex<NDIM> idx_c_low = idx.toLower();   // (i-1,j)
                CellIndex<NDIM> idx_c_up = idx.toUpper();    // (i,j)
                SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
                SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)
                SideIndex<NDIM> l_y_idx(idx_c_low, 1, 0);    // (i-1,j-1/2)
                SideIndex<NDIM> u_y_idx(idx_c_low, 1, 1);    // (i-1,j+1/2)

                // thn at sides
                double thn_lower = 0.5 * (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up);     // thn(i-1/2,j)
                double thn_upper = 0.5 * (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + xp); // thn(i+1/2,j)

                // thn at corners
                double thn_imhalf_jphalf =
                    0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + yp) +
                            (*thn_data)(idx_c_low + yp)); // thn(i-1/2,j+1/2)
                double thn_imhalf_jmhalf =
                    0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up - yp) +
                            (*thn_data)(idx_c_low - yp)); // thn(i-1/2,j-1/2)

                // components of first row (x-component of network vel) of network equation
                double ddx_Thn_dx_un = eta_n / (dx[0] * dx[0]) *
                                       ((*thn_data)(idx_c_up) * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                        (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - xp)));
                double ddy_Thn_dy_un = eta_n / (dx[1] * dx[1]) *
                                       (thn_imhalf_jphalf * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - yp)));
                double ddy_Thn_dx_vn = eta_n / (dx[1] * dx[0]) *
                                       (thn_imhalf_jphalf * ((*un_data)(upper_y_idx) - (*un_data)(u_y_idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(lower_y_idx) - (*un_data)(l_y_idx)));
                double ddx_Thn_dy_vn = eta_n / (dx[0] * dx[1]) *
                                       ((*thn_data)(idx_c_up) * ((*un_data)(upper_y_idx) - (*un_data)(lower_y_idx)) -
                                        (*thn_data)(idx_c_low) * ((*un_data)(u_y_idx) - (*un_data)(l_y_idx)));

                double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
                double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*A_un_data)(idx) = ddx_Thn_dx_un + ddy_Thn_dy_un + ddy_Thn_dx_vn + ddx_Thn_dy_vn + drag_n + pressure_n;

                // solvent equation
                double ddx_Ths_dx_us =
                    eta_s / (dx[0] * dx[0]) *
                    (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                     convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - xp)));
                double ddy_Ths_dy_us =
                    eta_s / (dx[1] * dx[1]) *
                    convertToThs((thn_imhalf_jphalf) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                 convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - yp)));
                double ddy_Ths_dx_vs =
                    eta_s / (dx[1] * dx[0]) *
                    convertToThs((thn_imhalf_jphalf) * ((*us_data)(upper_y_idx) - (*us_data)(u_y_idx)) -
                                 convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_y_idx) - (*us_data)(l_y_idx)));
                double ddx_Ths_dy_vs =
                    eta_s / (dx[0] * dx[1]) *
                    (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_y_idx) - (*us_data)(lower_y_idx)) -
                     convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_y_idx) - (*us_data)(l_y_idx)));

                double drag_s = -xi / nu_n * thn_lower * convertToThs(ths_lower) * ((*us_data)(idx) - (*un_data)(idx));
                double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*A_us_data)(idx) = ddx_Ths_dx_us + ddy_Ths_dy_us + ddy_Ths_dx_vs + ddx_Ths_dy_vs + drag_s + pressure_s;
            }

            for (SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
            {
                const SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)

                CellIndex<NDIM> idx_c_low = idx.toLower();   // (i,j-1)
                CellIndex<NDIM> idx_c_up = idx.toUpper();    // (i,j)
                SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
                SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
                SideIndex<NDIM> l_x_idx(idx_c_low, 0, 0);    // (i-1/2,j-1)
                SideIndex<NDIM> u_x_idx(idx_c_low, 0, 1);    // (i+1/2,j-1)

                // thn at sides
                double thn_lower = 0.5 * (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up);     // thn(i,j-1/2)
                double thn_upper = 0.5 * (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + yp); // thn(i,j+1/2)

                // thn at corners
                double thn_iphalf_jmhalf =
                    0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + xp) +
                            (*thn_data)(idx_c_low + xp)); // thn(i+1/2,j-1/2)
                double thn_imhalf_jmhalf =
                    0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up - xp) +
                            (*thn_data)(idx_c_low - xp)); // thn(i-1/2,j-1/2)

                // components of second row (y-component of network vel) of network equation
                double ddy_Thn_dy_un = eta_n / (dx[0] * dx[0]) *
                                       ((*thn_data)(idx_c_up) * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                        (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - yp)));
                double ddx_Thn_dx_un = eta_n / (dx[1] * dx[1]) *
                                       (thn_iphalf_jmhalf * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - xp)));
                double ddx_Thn_dy_vn = eta_n / (dx[1] * dx[0]) *
                                       (thn_iphalf_jphalf * ((*un_data)(upper_x_idx) - (*un_data)(u_x_idx)) -
                                        thn_imhalf_jphalf * ((*un_data)(lower_x_idx) - (*un_data)(l_x_idx)));
                double ddy_Thn_dx_vn = eta_n / (dx[0] * dx[1]) *
                                       ((*thn_data)(idx_c_up) * ((*un_data)(upper_x_idx) - (*un_data)(lower_x_idx)) -
                                        (*thn_data)(idx_c_low) * ((*un_data)(u_x_idx) - (*un_data)(l_x_idx)));

                double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
                double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*A_un_data)(idx) = ddy_Thn_dy_un + ddx_Thn_dx_un + ddx_Thn_dy_vn + ddy_Thn_dx_vn + drag_n + pressure_n;

                // Solvent equation
                double ddy_Ths_dy_us =
                    eta_s / (dx[0] * dx[0]) *
                    (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                     convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - yp)));
                double ddx_Ths_dx_us = eta_s / (dx[1] * dx[1]) *
                                       (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                                        convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - xp)));
                double ddx_Ths_dy_vs =
                    eta_s / (dx[1] * dx[0]) *
                    (convertToThs(thn_iphalf_jphalf) * ((*us_data)(upper_x_idx) - (*us_data)(u_x_idx)) -
                     convertToThs(thn_imhalf_jphalf) * ((*us_data)(lower_x_idx) - (*us_data)(l_x_idx)));
                double ddy_Ths_dx_vs =
                    eta_s / (dx[0] * dx[1]) *
                    (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_x_idx) - (*us_data)(lower_x_idx)) -
                     convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_x_idx) - (*us_data)(l_x_idx)));

                double drag_s = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*us_data)(idx) - (*un_data)(idx));
                double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*A_us_data)(idx) = ddy_Ths_dy_us + ddx_Ths_dx_us + ddx_Ths_dy_vs + ddy_Ths_dx_vs + drag_s + pressure_s;
            }
        }
    }
    if (d_bc_helper) d_bc_helper->copyDataAtDirichletBoundaries(A_un_idx, un_scratch_idx);

    IBAMR_TIMER_STOP(t_apply);
    return;
} // apply

double
convertToThs(double Thn)
{
    return 1 - Thn; // Thn+Ths = 1
}

void
VCTwoFluidStaggeredStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in.cloneVector(in.getName());
    d_b = out.cloneVector(out.getName());

    // Allocate scratch data.
    d_x->allocateVectorData();

    // Setup the interpolation transaction information.
    d_un_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    d_us_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    d_P_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.resize(4);
    d_transaction_comps[0] = InterpolationTransactionComponent(d_x->getComponentDescriptorIndex(0),
                                                               in.getComponentDescriptorIndex(0),
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_un_bc_coefs,
                                                               d_un_fill_pattern);
    d_transaction_comps[1] = InterpolationTransactionComponent(d_x->getComponentDescriptorIndex(1),
                                                               in.getComponentDescriptorIndex(1),
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_us_bc_coefs,
                                                               d_us_fill_pattern);
    d_transaction_comps[2] = InterpolationTransactionComponent(in.getComponentDescriptorIndex(2),
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_P_bc_coef,
                                                               d_P_fill_pattern);
    d_transaction_comps[3] = InterpolationTransactionComponent(thn_idx,
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               nullptr,
                                                               d_P_fill_pattern);

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

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
VCTwoFluidStaggeredStokesOperator::deallocateOperatorState()
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

    // Delete the solution and rhs vectors.
    d_x->resetLevels(d_x->getCoarsestLevelNumber(),
                     std::min(d_x->getFinestLevelNumber(), d_x->getPatchHierarchy()->getFinestLevelNumber()));
    d_x->freeVectorComponents();

    d_b->resetLevels(d_b->getCoarsestLevelNumber(),
                     std::min(d_b->getFinestLevelNumber(), d_b->getPatchHierarchy()->getFinestLevelNumber()));
    d_b->freeVectorComponents();

    d_x.setNull();
    d_b.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

void
VCTwoFluidStaggeredStokesOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
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
            const int P_idx = x->getComponentDescriptorIndex(2);
            StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
                d_un_bc_coefs, d_P_bc_coef, un_idx, P_idx, d_homogeneous_bc);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                un_idx, P_idx, d_un_bc_coefs, d_new_time, d_homogeneous_bc);
        }
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_un_bc_coefs, d_P_bc_coef);
        apply(*x, *b);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&y, false), b);
        x->freeVectorComponents();
        b->freeVectorComponents();
    }
    const bool homogeneous_bc = true;
    if (d_bc_helper)
    {
        const int un_idx = y.getComponentDescriptorIndex(0);
        const int P_idx = y.getComponentDescriptorIndex(2);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_un_bc_coefs, d_P_bc_coef, un_idx, P_idx, homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(un_idx, P_idx, d_un_bc_coefs, d_new_time, homogeneous_bc);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_un_bc_coefs, d_P_bc_coef);
    }
    return;
} // modifyRhsForBcs

void
VCTwoFluidStaggeredStokesOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    if (d_bc_helper)
    {
        const int un_idx = u.getComponentDescriptorIndex(0);
        const int P_idx = u.getComponentDescriptorIndex(2);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_un_bc_coefs, d_P_bc_coef, un_idx, P_idx, d_homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(
            un_idx, P_idx, d_un_bc_coefs, d_new_time, d_homogeneous_bc);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_un_bc_coefs, d_P_bc_coef);
    }
    return;
} // imposeSolBcs

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
