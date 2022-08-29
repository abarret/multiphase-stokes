// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2020 by the IBAMR developers
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

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/CartCellDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "Box.h"
#include "BoxList.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "ProcessorMapping.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"

#include "petscksp.h"

#include <Eigen/LU>

#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

// Local includes
#include "VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////
namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
// Timers.
static Timer* t_smooth_error;
static Timer* t_solve_coarsest_level;
static Timer* t_compute_residual;

// Default data depth.
static const int DEFAULT_DATA_DEPTH = 1;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string CC_DATA_REFINE_TYPE =
    "CONSERVATIVE_LINEAR_REFINE"; // how to fill in fine cells from coarse cells, how to fill ghost cells on refine
                                  // patch
static const std::string SC_DATA_REFINE_TYPE =
    "CONSERVATIVE_LINEAR_REFINE"; // how to fill in fine cells from coarse cells, how to fill ghost cells on refine
                                  // patch
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;
} // namespace
double convertToThs(double Thn);
/////////////////////////////// PUBLIC ///////////////////////////////////////

VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : FACPreconditionerStrategy(object_name)
{
    // Create variables and register them with the variable database.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> d_ctx = var_db->getContext("context");

    // State scratch variables: Velocity and pressure.
    Pointer<SideVariable<NDIM, double>> un_scr_var = new SideVariable<NDIM, double>("un_scr");
    Pointer<SideVariable<NDIM, double>> us_scr_var = new SideVariable<NDIM, double>("us_scr");
    Pointer<CellVariable<NDIM, double>> p_scr_var = new CellVariable<NDIM, double>("p_scr");

    // Register patch data indices
    const int un_scr_idx = var_db->registerVariableAndContext(un_scr_var, d_ctx, IntVector<NDIM>(1));
    const int us_scr_idx = var_db->registerVariableAndContext(us_scr_var, d_ctx, IntVector<NDIM>(1));
    const int p_scr_idx = var_db->registerVariableAndContext(p_scr_var, d_ctx, IntVector<NDIM>(1));

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error = TimerManager::getManager()->getTimer(
                     "IBTK::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::computeResidual()"););
    return;
}

VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::~VCTwoFluidStaggeredStokesBoxRelaxationFACOperator()
{
    return;
}

// create another member function to set-up Thn
// Thn is defined in the input file and read in using muParserCartGridFunction
void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::setThnIdx(int thn_idx)
{
    d_thn_idx = thn_idx;
    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::restrictResidual(const SAMRAIVectorReal<NDIM, double>& source,
                                                                    SAMRAIVectorReal<NDIM, double>& dest,
                                                                    int dest_level_num)
{
    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::prolongError(const SAMRAIVectorReal<NDIM, double>& source,
                                                                SAMRAIVectorReal<NDIM, double>& dest,
                                                                int dest_level_num)
{
    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::prolongErrorAndCorrect(const SAMRAIVectorReal<NDIM, double>& source,
                                                                          SAMRAIVectorReal<NDIM, double>& dest,
                                                                          int dest_level_num)
{
    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::smoothError(
    SAMRAIVectorReal<NDIM, double>& error,          // Solution
    const SAMRAIVectorReal<NDIM, double>& residual, // RHS - A*error on the entire patch level
    int level_num,
    int num_sweeps,
    bool /*performing_pre_sweeps*/,
    bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

    // Get the vector components. These pull out patch data indices
    const int un_idx = error.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = error.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int P_idx = error.getComponentDescriptorIndex(2);  // pressure
    const int thn_idx = d_thn_idx;
    const int f_un_idx = residual.getComponentDescriptorIndex(0); // RHS_Un
    const int f_us_idx = residual.getComponentDescriptorIndex(1); // RHS_Us
    const int f_P_idx = residual.getComponentDescriptorIndex(2);  // RHS_pressure

    // Create scratch indices
    // Copy values from the residual at step 0 into scratch indices
    // We want to solve for Un, Us and P at step k+1 using data from step k
    // Allocate data on all patches. (better to do this in InitializeSolverState)
    // d_hierarchy = error.getPatchHierarchy();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);
    level->allocatePatchData(un_scr_idx, 0.0);
    level->allocatePatchData(us_scr_idx, 0.0);
    level->allocatePatchData(p_scr_idx, 0.0);

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(4);
    transaction_comps[0] = InterpolationTransactionComponent(un_idx,
                                                             un_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_un_bc_coefs, // modifiy?
                                                             d_un_fill_pattern);
    transaction_comps[1] = InterpolationTransactionComponent(us_idx,
                                                             us_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_us_bc_coefs, // modify?
                                                             d_us_fill_pattern);
    transaction_comps[2] = InterpolationTransactionComponent(P_idx,
                                                             CC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef,
                                                             d_P_fill_pattern);
    transaction_comps[3] = InterpolationTransactionComponent(thn_idx,
                                                             CC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef); // defaults to fill corner

    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time); // Fills in all of the ghost cells

    const double eta_n = 1.0;
    const double eta_s = 1.0;
    const double xi = 1.0;
    const double nu_n = 1.0;
    const double nu_s = 1.0;
    // outer for loop for number of sweeps
    for (int sweep = 0; sweep <= num_sweeps; sweep++)
    {
        // loop through all patches on this level
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
            const double* const xlow =
                pgeom->getXLower(); // {xlow[0], xlow[1]} -> physical location of bottom left of box.
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(P_idx);
            Pointer<SideData<NDIM, double>> un_scr_data = patch->getPatchData(un_scr_idx);
            Pointer<SideData<NDIM, double>> us_scr_data = patch->getPatchData(us_scr_idx);
            Pointer<CellData<NDIM, double>> p_scr_data = patch->getPatchData(p_scr_idx);
            Pointer<SideData<NDIM, double>> f_un_data = patch->getPatchData(f_un_idx);
            Pointer<SideData<NDIM, double>> f_us_data = patch->getPatchData(f_us_idx);
            Pointer<CellData<NDIM, double>> f_p_data = patch->getPatchData(f_P_idx);
            IntVector<NDIM> xp(1, 0), yp(0, 1);
            // un_scr_data->fillAll(*un_data);
            // us_scr_data->fillAll(*us_data);
            // p_scr_data->fillAll(*p_data);

            MatrixXd A_box(9, 9); // 9 x 9 Matrix
            VectorXd b;           // 9 x 1 RHS vector
            VectorXd sol;         // 9 x 1 solution vector

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++) // cell-centers
            {
                un_scr_data->fillAll(un_data);
                us_scr_data->fillAll(us_data);
                p_scr_data->fillAll(p_data);

                const CellIndex<NDIM>& idx = ci(); // (i,j)
                VectorNd x;                        // <- Eigen3 vector
                for (int d = 0; d < NDIM; ++d)
                    x[d] = xlow[d] + dx[d] * (idx(d) - idx_low(d) + 0.5); // Get's physical location of idx.

                SideIndex<NDIM> lower_x_idx(idx, 0, 0); // (i-1/2,j)
                SideIndex<NDIM> upper_x_idx(idx, 0, 1); // (i+1/2,j)
                SideIndex<NDIM> lower_y_idx(idx, 1, 0); // (i,j-1/2)
                SideIndex<NDIM> upper_y_idx(idx, 1, 1); // (i,j+1/2)

                // thn at sidess
                double thn_lower_x = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx - xp)); // thn(i-1/2,j)
                double thn_upper_x = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx + xp)); // thn(i+1/2,j)
                double thn_lower_y = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx - yp)); // thn(i,j-1/2)
                double thn_upper_y = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx + yp)); // thn(i,j+1/2)

                // thn at corners
                double thn_imhalf_jphalf = 0.25 * ((*thn_data)(idx - xp) + (*thn_data)(idx) + (*thn_data)(idx + yp) +
                                                   (*thn_data)(idx - xp + yp)); // thn(i-1/2,j+1/2)
                double thn_imhalf_jmhalf = 0.25 * ((*thn_data)(idx) + (*thn_data)(idx - xp) + (*thn_data)(idx - yp) +
                                                   (*thn_data)(idx - xp - yp)); // thn(i-1/2,j-1/2)
                double thn_iphalf_jphalf = 0.25 * ((*thn_data)(idx + xp) + (*thn_data)(idx) + (*thn_data)(idx + yp) +
                                                   (*thn_data)(idx + xp + yp)); // thn(i+1/2,j+1/2)
                double thn_iphalf_jmhalf = 0.25 * ((*thn_data)(idx + xp) + (*thn_data)(idx) + (*thn_data)(idx - yp) +
                                                   (*thn_data)(idx + xp - yp)); // thn(i+1/2,j-1/2)

                // network at west edge
                A_box(0, 0) = eta_n / (dx[0] * dx[0]) * (-(*thn_data)(idx) - (*thn_data)(idx - xp)) -
                              eta_n / (dx[1] * dx[1]) * (thn_imhalf_jmhalf + thn_imhalf_jphalf) -
                              xi / nu_n * thn_lower_x * convertToThs(thn_lower_x);
                A_box(0, 1) = eta_n / (dx[0] * dx[0]) * ((*thn_data)(idx));
                A_box(0, 2) = eta_n / (dx[0] * dx[1]) * ((*thn_data)(idx)-thn_imhalf_jmhalf);
                A_box(0, 3) = eta_n / (dx[0] * dx[1]) * (thn_imhalf_jphalf - (*thn_data)(idx));
                A_box(0, 4) = xi / nu_n * thn_lower_x * convertToThs(thn_lower_x);
                A_box(0, 5) = A_box(0, 6) = A_box(0, 7) = 0.0;
                A_box(0, 8) = -thn_lower_x / dx[0];

                A_box(4, 4) =
                    eta_s / (dx[0] * dx[0]) * (-convertToThs((*thn_data)(idx)) - convertToThs((*thn_data)(idx - xp))) -
                    eta_s / (dx[1] * dx[1]) * (convertToThs(thn_imhalf_jmhalf) + convertToThs(thn_imhalf_jphalf)) -
                    xi / nu_n * thn_lower_x * convertToThs(thn_lower_x);
                A_box(4, 5) = eta_s / (dx[0] * dx[0]) * (convertToThs((*thn_data)(idx)));
                A_box(4, 6) =
                    eta_s / (dx[0] * dx[1]) * (convertToThs((*thn_data)(idx)) - convertToThs(thn_imhalf_jmhalf));
                A_box(4, 7) =
                    eta_s / (dx[0] * dx[1]) * (convertToThs(thn_imhalf_jphalf) - convertToThs((*thn_data)(idx)));
                A_box(4, 0) = xi / nu_s * thn_lower_x * convertToThs(thn_lower_x);
                A_box(4, 1) = A_box(4, 2) = A_box(4, 3) = 0.0;
                A_box(4, 8) = -convertToThs(thn_lower_x) / dx[0];

                // network at east edge
                A_box(1, 0) = eta_n / (dx[0] * dx[0]) * (*thn_data)(idx);
                A_box(1, 1) = eta_n / (dx[0] * dx[0]) * (-(*thn_data)(idx + xp) - (*thn_data)(idx)) -
                              eta_n / (dx[1] * dx[1]) * (thn_iphalf_jphalf + thn_iphalf_jmhalf) -
                              xi / nu_n * thn_upper_x * convertToThs(thn_upper_x);
                A_box(1, 2) = eta_n / (dx[1] * dx[0]) * (thn_iphalf_jmhalf - (*thn_data)(idx));
                A_box(1, 3) = eta_n / (dx[1] * dx[0]) * ((*thn_data)(idx)-thn_iphalf_jphalf);
                A_box(1, 4) = A_box(1, 6) = A_box(1, 7) = 0.0;
                A_box(1, 5) = xi / nu_n * thn_upper_x * convertToThs(thn_upper_x);
                A_box(1, 8) = thn_upper_x / dx[0];

                A_box(5, 4) = eta_s / (dx[0] * dx[0]) * convertToThs((*thn_data)(idx));
                A_box(5, 5) =
                    eta_s / (dx[0] * dx[0]) * (-convertToThs((*thn_data)(idx + xp)) - convertToThs((*thn_data)(idx))) -
                    eta_s / (dx[1] * dx[1]) * (convertToThs(thn_iphalf_jphalf) + convertToThs(thn_iphalf_jmhalf)) -
                    xi / nu_s * thn_upper_x * convertToThs(thn_upper_x);
                A_box(5, 6) =
                    eta_n / (dx[1] * dx[0]) * (convertToThs(thn_iphalf_jmhalf) - convertToThs((*thn_data)(idx)));
                A_box(5, 7) =
                    eta_n / (dx[1] * dx[0]) * (convertToThs((*thn_data)(idx)) - convertToThs(thn_iphalf_jphalf));
                A_box(5, 0) = A_box(5, 2) = A_box(5, 3) = 0.0;
                A_box(5, 1) = xi / nu_s * thn_upper_x * convertToThs(thn_upper_x);
                A_box(5, 8) = convertToThs(thn_upper_x) / dx[0];

                // network at south edge
                A_box(2, 0) = eta_n / (dx[0] * dx[1]) * ((*thn_data)(idx)-thn_imhalf_jmhalf);
                A_box(2, 1) = eta_n / (dx[0] * dx[1]) * (thn_iphalf_jmhalf - (*thn_data)(idx));
                A_box(2, 2) = eta_n / (dx[1] * dx[1]) * (-(*thn_data)(idx) - (*thn_data)(idx - yp)) -
                              eta_n / (dx[0] * dx[0]) * (thn_iphalf_jmhalf + thn_imhalf_jmhalf) -
                              xi / nu_n * thn_lower_y * convertToThs(thn_lower_y);
                A_box(2, 3) = eta_n / (dx[1] * dx[1]) * ((*thn_data)(idx));
                A_box(2, 4) = A_box(2, 5) = A_box(2, 7) = 0.0;
                A_box(2, 6) = xi / nu_n * thn_lower_y * convertToThs(thn_lower_y);
                A_box(2, 8) = -thn_lower_y / dx[1];

                A_box(6, 4) =
                    eta_s / (dx[0] * dx[1]) * (convertToThs((*thn_data)(idx)) - convertToThs(thn_imhalf_jmhalf));
                A_box(6, 5) =
                    eta_s / (dx[0] * dx[1]) * (convertToThs(thn_iphalf_jmhalf) - convertToThs((*thn_data)(idx)));
                A_box(6, 6) =
                    eta_s / (dx[1] * dx[1]) * (-convertToThs((*thn_data)(idx)) - convertToThs((*thn_data)(idx - yp))) -
                    eta_s / (dx[0] * dx[0]) * (convertToThs(thn_iphalf_jmhalf) + convertToThs(thn_imhalf_jmhalf)) -
                    xi / nu_n * thn_lower_y * convertToThs(thn_lower_y);
                A_box(6, 7) = eta_s / (dx[1] * dx[1]) * (convertToThs((*thn_data)(idx)));
                A_box(6, 0) = A_box(6, 1) = A_box(6, 3) = 0.0;
                A_box(6, 2) = xi / nu_s * thn_lower_y * convertToThs(thn_lower_y);
                A_box(6, 8) = -convertToThs(thn_lower_y) / dx[1];

                // network at north edge
                A_box(3, 0) = eta_n / (dx[0] * dx[1]) * (thn_imhalf_jphalf - (*thn_data)(idx));
                A_box(3, 1) = eta_n / (dx[0] * dx[1]) * ((*thn_data)(idx)-thn_iphalf_jphalf);
                A_box(3, 2) = eta_n / (dx[1] * dx[1]) * ((*thn_data)(idx));
                A_box(3, 3) = eta_n / (dx[1] * dx[1]) * (-(*thn_data)(idx) - (*thn_data)(idx + yp)) -
                              eta_n / (dx[0] * dx[0]) * (thn_iphalf_jphalf - thn_imhalf_jphalf) -
                              xi / nu_n * thn_upper_y * convertToThs(thn_upper_y);
                A_box(3, 4) = A_box(3, 5) = A_box(3, 6) = 0.0;
                A_box(3, 7) = xi / nu_n * thn_upper_y * convertToThs(thn_upper_y);
                A_box(3, 8) = thn_upper_y / dx[1];

                A_box(7, 4) =
                    eta_s / (dx[0] * dx[1]) * (convertToThs(thn_imhalf_jphalf) - convertToThs((*thn_data)(idx)));
                A_box(7, 5) =
                    eta_s / (dx[0] * dx[1]) * (convertToThs((*thn_data)(idx)) - convertToThs(thn_iphalf_jphalf));
                A_box(7, 6) = eta_s / (dx[1] * dx[1]) * (convertToThs((*thn_data)(idx)));
                A_box(7, 7) =
                    eta_s / (dx[1] * dx[1]) * (-convertToThs((*thn_data)(idx)) - convertToThs((*thn_data)(idx + yp))) -
                    eta_s / (dx[0] * dx[0]) * (convertToThs(thn_iphalf_jphalf) - convertToThs(thn_imhalf_jphalf)) -
                    xi / nu_s * thn_upper_y * convertToThs(thn_upper_y);
                A_box(7, 0) = A_box(7, 1) = A_box(7, 2) = 0.0;
                A_box(7, 3) = xi / nu_s * thn_upper_y * convertToThs(thn_upper_y);
                A_box(7, 8) = convertToThs(thn_upper_y) / dx[1];

                // incompressible constrain term at center
                A_box(8, 0) = thn_lower_x / dx[0];
                A_box(8, 1) = -thn_upper_x / dx[0];
                A_box(8, 2) = thn_lower_y / dx[1];
                A_box(8, 3) = -thn_upper_y / dx[1];
                A_box(8, 4) = convertToThs(thn_lower_x) / dx[0];
                A_box(8, 5) = -convertToThs(thn_upper_x) / dx[0];
                A_box(8, 6) = convertToThs(thn_lower_y) / dx[1];
                A_box(8, 7) = -convertToThs(thn_upper_y) / dx[1];
                A_box(8, 8) = 0.0;

                // set-up RHS vector (include terms from residual (f_un,f_us,f_p))
                // should be populated with the values from the previous iteration

                // network at west edge
                b(0) = (*f_un_data)(lower_x_idx)-thn_lower_x / dx[0] * (*p_scr_data)(idx - xp) -
                       eta_n / (dx[0] * dx[0]) * (*thn_data)(idx - xp) * (*un_scr_data)(lower_x_idx - xp) +
                       eta_n / (dx[1] * dx[1]) * thn_imhalf_jphalf * (*un_scr_data)(lower_x_idx + yp) -
                       eta_n / (dx[1] * dx[1]) * thn_imhalf_jmhalf * (*un_scr_data)(lower_x_idx - yp) +
                       eta_n / (dx[0] * dx[1]) * thn_imhalf_jphalf * (*un_scr_data)(upper_y_idx - xp) -
                       eta_n / (dx[0] * dx[1]) * thn_imhalf_jmhalf * (*un_scr_data)(lower_y_idx - xp) -
                       eta_n / (dx[0] * dx[1]) * (*thn_data)(idx - xp) *
                           ((*un_scr_data)(upper_y_idx - xp) - (*un_scr_data)(lower_y_idx - xp));

                // solvent at west edge
                b(4) =
                    (*f_us_data)(lower_x_idx)-convertToThs(thn_lower_x) / dx[0] * (*p_scr_data)(idx - xp) -
                    eta_s / (dx[0] * dx[0]) * convertToThs((*thn_data)(idx - xp)) * (*us_scr_data)(lower_x_idx - xp) +
                    eta_s / (dx[1] * dx[1]) * convertToThs(thn_imhalf_jphalf) * (*us_scr_data)(lower_x_idx + yp) -
                    eta_s / (dx[1] * dx[1]) * convertToThs(thn_imhalf_jmhalf) * (*us_scr_data)(lower_x_idx - yp) +
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_imhalf_jphalf) * (*us_scr_data)(upper_y_idx - xp) -
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_imhalf_jmhalf) * (*us_scr_data)(lower_y_idx - xp) -
                    eta_s / (dx[0] * dx[1]) * convertToThs((*thn_data)(idx - xp)) *
                        ((*us_scr_data)(upper_y_idx - xp) - (*us_scr_data)(lower_y_idx - xp));

                // network at east edge
                b(1) = (*f_un_data)(upper_x_idx) + thn_upper_x / dx[0] * (*p_scr_data)(idx + xp) -
                       eta_n / (dx[0] * dx[0]) * (*thn_data)(idx + xp) * (*un_scr_data)(upper_x_idx + xp) -
                       eta_n / (dx[1] * dx[1]) * thn_iphalf_jphalf * (*un_scr_data)(upper_x_idx + yp) -
                       eta_n / (dx[1] * dx[1]) * thn_iphalf_jmhalf * (*un_scr_data)(upper_x_idx - yp) -
                       eta_n / (dx[0] * dx[1]) * thn_iphalf_jphalf * (*un_scr_data)(upper_y_idx + xp) +
                       eta_n / (dx[0] * dx[1]) * thn_iphalf_jmhalf * (*un_scr_data)(lower_y_idx + xp) +
                       eta_n / (dx[0] * dx[1]) * (*thn_data)(idx + xp) *
                           ((*un_scr_data)(upper_y_idx + xp) - (*un_scr_data)(lower_y_idx + xp));

                // solvent at east edge
                b(5) =
                    (*f_us_data)(upper_x_idx) + convertToThs(thn_upper_x) / dx[0] * (*p_scr_data)(idx + xp) -
                    eta_s / (dx[0] * dx[0]) * convertToThs((*thn_data)(idx + xp)) * (*us_scr_data)(upper_x_idx + xp) -
                    eta_s / (dx[1] * dx[1]) * convertToThs(thn_iphalf_jphalf) * (*us_scr_data)(upper_x_idx + yp) -
                    eta_s / (dx[1] * dx[1]) * convertToThs(thn_iphalf_jmhalf) * (*us_scr_data)(upper_x_idx - yp) -
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_iphalf_jphalf) * (*us_scr_data)(upper_y_idx + xp) +
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_iphalf_jmhalf) * (*us_scr_data)(lower_y_idx + xp) +
                    eta_s / (dx[0] * dx[1]) * convertToThs((*thn_data)(idx + xp)) *
                        ((*us_scr_data)(upper_y_idx + xp) - (*us_scr_data)(lower_y_idx + xp));

                // network at south edge
                b(2) = (*f_un_data)(lower_y_idx)-thn_lower_y / dx[1] * (*p_scr_data)(idx - yp) -
                       eta_n / (dx[1] * dx[1]) * (*thn_data)(idx - yp) * (*un_scr_data)(lower_y_idx - yp) -
                       eta_n / (dx[0] * dx[0]) * thn_iphalf_jmhalf * (*un_scr_data)(lower_y_idx + xp) -
                       eta_n / (dx[0] * dx[0]) * thn_imhalf_jmhalf * (*un_scr_data)(lower_y_idx - xp) +
                       eta_n / (dx[0] * dx[1]) * thn_iphalf_jmhalf * (*un_scr_data)(upper_x_idx - yp) -
                       eta_n / (dx[0] * dx[1]) * thn_imhalf_jmhalf * (*un_scr_data)(lower_x_idx - yp) -
                       eta_n / (dx[0] * dx[1]) * (*thn_data)(idx - yp) *
                           ((*un_scr_data)(upper_x_idx - yp) - (*un_scr_data)(lower_x_idx - yp));

                // solvent at south edge
                b(6) =
                    (*f_us_data)(lower_y_idx)-convertToThs(thn_lower_y) / dx[1] * (*p_scr_data)(idx - yp) -
                    eta_s / (dx[1] * dx[1]) * convertToThs((*thn_data)(idx - yp)) * (*us_scr_data)(lower_y_idx - yp) -
                    eta_s / (dx[0] * dx[0]) * convertToThs(thn_iphalf_jmhalf) * (*us_scr_data)(lower_y_idx + xp) -
                    eta_s / (dx[0] * dx[0]) * convertToThs(thn_imhalf_jmhalf) * (*us_scr_data)(lower_y_idx - xp) +
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_iphalf_jmhalf) * (*us_scr_data)(upper_x_idx - yp) -
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_imhalf_jmhalf) * (*us_scr_data)(lower_x_idx - yp) -
                    eta_s / (dx[0] * dx[1]) * convertToThs((*thn_data)(idx - yp)) *
                        ((*us_scr_data)(upper_x_idx - yp) - (*us_scr_data)(lower_x_idx - yp));

                // network at north edge
                b(3) = (*f_un_data)(upper_y_idx) + thn_upper_y / dx[1] * (*p_scr_data)(idx + yp) -
                       eta_n / (dx[1] * dx[1]) * (*thn_data)(idx + yp) * (*un_scr_data)(upper_y_idx + yp) -
                       eta_n / (dx[0] * dx[0]) * thn_iphalf_jphalf * (*un_scr_data)(upper_y_idx + xp) -
                       eta_n / (dx[0] * dx[0]) * thn_imhalf_jphalf * (*un_scr_data)(upper_y_idx - xp) -
                       eta_n / (dx[0] * dx[1]) * thn_iphalf_jphalf * (*un_scr_data)(upper_x_idx + yp) +
                       eta_n / (dx[0] * dx[1]) * thn_imhalf_jphalf * (*un_scr_data)(lower_x_idx + yp) +
                       eta_n / (dx[0] * dx[1]) * (*thn_data)(idx + yp) *
                           ((*un_scr_data)(upper_x_idx + yp) - (*un_scr_data)(lower_x_idx + yp));

                // solvent at north edge
                b(7) =
                    (*f_us_data)(upper_y_idx) + convertToThs(thn_upper_y) / dx[1] * (*p_scr_data)(idx + yp) -
                    eta_s / (dx[1] * dx[1]) * convertToThs((*thn_data)(idx + yp)) * (*us_scr_data)(upper_y_idx + yp) -
                    eta_s / (dx[0] * dx[0]) * convertToThs(thn_iphalf_jphalf) * (*us_scr_data)(upper_y_idx + xp) -
                    eta_s / (dx[0] * dx[0]) * convertToThs(thn_imhalf_jphalf) * (*us_scr_data)(upper_y_idx - xp) -
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_iphalf_jphalf) * (*us_scr_data)(upper_x_idx + yp) +
                    eta_s / (dx[0] * dx[1]) * convertToThs(thn_imhalf_jphalf) * (*us_scr_data)(lower_x_idx + yp) +
                    eta_s / (dx[0] * dx[1]) * convertToThs((*thn_data)(idx + yp)) *
                        ((*us_scr_data)(upper_x_idx + yp) - (*us_scr_data)(lower_x_idx + yp));

                // pressure at cell center
                b(8) = (*f_p_data)(idx);

                sol = A_box.lu().solve(b); // solve Ax = b per cell
                // sol = A_box.partialPivLu().solve(b);
                (*un_data)(lower_x_idx) = sol(0);
                (*un_data)(upper_x_idx) = sol(1);
                (*un_data)(lower_y_idx) = sol(2);
                (*un_data)(upper_y_idx) = sol(3);
                (*us_data)(lower_x_idx) = sol(4);
                (*us_data)(upper_x_idx) = sol(5);
                (*us_data)(lower_y_idx) = sol(6);
                (*us_data)(upper_y_idx) = sol(7);
                (*p_data)(idx) = sol(8);

            } // cell centers
        }     // patches
    }         // num_sweeps
    IBTK_TIMER_STOP(t_smooth_error);
    return;
}

bool
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                                      const SAMRAIVectorReal<NDIM, double>& residual,
                                                                      int coarsest_ln)
{
    return false;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                                   const SAMRAIVectorReal<NDIM, double>& solution,
                                                                   const SAMRAIVectorReal<NDIM, double>& rhs,
                                                                   int coarsest_level_num,
                                                                   int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    // Get the vector components. These pull out patch data indices
    const int un_idx = solution.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = solution.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int P_idx = solution.getComponentDescriptorIndex(2);  // pressure
    const int rhs_un_idx = rhs.getComponentDescriptorIndex(0);  // RHS Un
    const int rhs_us_idx = rhs.getComponentDescriptorIndex(1);  // RHS Us
    const int rhs_P_idx = rhs.getComponentDescriptorIndex(2);   // RHS P
    const int res_un_idx = residual.getComponentDescriptorIndex(0);
    const int res_us_idx = residual.getComponentDescriptorIndex(1);
    const int res_P_idx = residual.getComponentDescriptorIndex(2);
    const int thn_idx = d_thn_idx;

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(4);
    transaction_comps[0] = InterpolationTransactionComponent(un_idx,
                                                             un_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_un_bc_coefs, // modifiy?
                                                             d_un_fill_pattern);
    transaction_comps[1] = InterpolationTransactionComponent(us_idx,
                                                             us_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_us_bc_coefs, // modify?
                                                             d_us_fill_pattern);
    transaction_comps[2] = InterpolationTransactionComponent(P_idx,
                                                             CC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef,
                                                             d_P_fill_pattern);
    transaction_comps[3] = InterpolationTransactionComponent(thn_idx,
                                                             CC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef); // defaults to fill corner

    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time); // Fills in all of the ghost cells
    const double eta_n = 1.0;
    const double eta_s = 1.0;
    const double xi = 1.0;
    const double nu_n = 1.0;
    const double nu_s = 1.0;

    for (int ln = coarsest_level_num; ln <= finest_level_num; ++ln)
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
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(P_idx);
            Pointer<CellData<NDIM, double>> rhs_P_data =
                patch->getPatchData(rhs_P_idx); // result of applying operator (eqn 3)
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> rhs_un_data =
                patch->getPatchData(rhs_un_idx); // result of applying operator (eqn 1)
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<SideData<NDIM, double>> rhs_us_data =
                patch->getPatchData(rhs_us_idx); // result of applying operator (eqn 2)
            Pointer<SideData<NDIM, double>> res_un_data = patch->getPatchData(res_un_idx);
            Pointer<SideData<NDIM, double>> res_us_data = patch->getPatchData(res_un_idx);
            Pointer<CellData<NDIM, double>> res_P_data = patch->getPatchData(res_P_idx);
            IntVector<NDIM> xp(1, 0), yp(0, 1);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++) // cell-centers
            {
                const CellIndex<NDIM>& idx = ci();
                VectorNd x; // <- Eigen3 vector
                for (int d = 0; d < NDIM; ++d)
                    x[d] = xlow[d] + dx[d] * (idx(d) - idx_low(d) + 0.5); // Get's physical location of idx.

                SideIndex<NDIM> lower_x_idx(idx, 0, 0); // (i-1/2,j)
                SideIndex<NDIM> upper_x_idx(idx, 0, 1); // (i+1/2,j)
                SideIndex<NDIM> lower_y_idx(idx, 1, 0); // (i,j-1/2)
                SideIndex<NDIM> upper_y_idx(idx, 1, 1); // (i,j+1/2)

                // thn at sidess
                double thn_lower_x = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx - xp)); // thn(i-1/2,j)
                double thn_upper_x = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx + xp)); // thn(i+1/2,j)
                double thn_lower_y = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx - yp)); // thn(i,j-1/2)
                double thn_upper_y = 0.5 * ((*thn_data)(idx) + (*thn_data)(idx + yp)); // thn(i,j+1/2)

                // conservation of mass

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
                (*res_P_data)(idx) = (*rhs_P_data)(idx) - (div_un_thn + div_us_ths);
            }

            for (SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
            {
                const SideIndex<NDIM>& idx = si(); // axis = 0, (i-1/2,j)
                // pout << "\n At idx " << idx << " \n ";

                CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
                CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
                SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
                SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)
                SideIndex<NDIM> l_y_idx(idx_c_low, 1, 0);    // (i-1,j-1/2)
                SideIndex<NDIM> u_y_idx(idx_c_low, 1, 1);    // (i-1,j+1/2)

                // thn at sides
                double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i-1/2,j)
                // pout << "thn_lower is " << thn_lower << "\n";
                // thn at corners
                double thn_imhalf_jphalf =
                    0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up + yp) +
                            (*thn_data)(idx_c_low + yp)); // thn(i-1/2,j+1/2)
                double thn_imhalf_jmhalf =
                    0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up - yp) +
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
                double ddx_Thn_dy_vn = -eta_n / (dx[0] * dx[1]) *
                                       ((*thn_data)(idx_c_up) * ((*un_data)(upper_y_idx) - (*un_data)(lower_y_idx)) -
                                        (*thn_data)(idx_c_low) * ((*un_data)(u_y_idx) - (*un_data)(l_y_idx)));

                double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
                double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*res_un_data)(idx) = (*rhs_un_data)(idx) - (ddx_Thn_dx_un + ddy_Thn_dy_un + ddy_Thn_dx_vn +
                                                             ddx_Thn_dy_vn + drag_n + pressure_n);

                // solvent equation
                double ddx_Ths_dx_us =
                    eta_s / (dx[0] * dx[0]) *
                    (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(idx + xp) - (*us_data)(idx)) -
                     convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(idx) - (*us_data)(idx - xp)));
                double ddy_Ths_dy_us = eta_s / (dx[1] * dx[1]) *
                                       (convertToThs(thn_imhalf_jphalf) * ((*us_data)(idx + yp) - (*us_data)(idx)) -
                                        convertToThs(thn_imhalf_jmhalf) * ((*us_data)(idx) - (*us_data)(idx - yp)));
                double ddy_Ths_dx_vs =
                    eta_s / (dx[1] * dx[0]) *
                    (convertToThs(thn_imhalf_jphalf) * ((*us_data)(upper_y_idx) - (*us_data)(u_y_idx)) -
                     convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_y_idx) - (*us_data)(l_y_idx)));
                double ddx_Ths_dy_vs =
                    -eta_s / (dx[0] * dx[1]) *
                    (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_y_idx) - (*us_data)(lower_y_idx)) -
                     convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_y_idx) - (*us_data)(l_y_idx)));

                double drag_s = -xi / nu_s * thn_lower * convertToThs(thn_lower) * ((*us_data)(idx) - (*un_data)(idx));
                double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*res_us_data)(idx) = (*rhs_us_data)(idx) - (ddx_Ths_dx_us + ddy_Ths_dy_us + ddy_Ths_dx_vs +
                                                             ddx_Ths_dy_vs + drag_s + pressure_s);
            }

            // pout << "\n\n Looping over y-dir side-centers \n\n";

            for (SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
            {
                const SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)

                CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
                CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
                SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
                SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
                SideIndex<NDIM> l_x_idx(idx_c_low, 0, 0);    // (i-1/2,j-1)
                SideIndex<NDIM> u_x_idx(idx_c_low, 0, 1);    // (i+1/2,j-1)

                // thn at sides
                double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i,j-1/2)

                // thn at corners
                double thn_imhalf_jmhalf =
                    0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up - xp) +
                            (*thn_data)(idx_c_low - xp)); // thn(i-1/2,j-1/2)
                double thn_iphalf_jmhalf =
                    0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up + xp) +
                            (*thn_data)(idx_c_low + xp)); // thn(i+1/2,j-1/2)

                // components of second row (y-component of network vel) of network equation
                double ddy_Thn_dy_un = eta_n / (dx[0] * dx[0]) *
                                       ((*thn_data)(idx_c_up) * ((*un_data)(idx + yp) - (*un_data)(idx)) -
                                        (*thn_data)(idx_c_low) * ((*un_data)(idx) - (*un_data)(idx - yp)));
                double ddx_Thn_dx_un = eta_n / (dx[1] * dx[1]) *
                                       (thn_iphalf_jmhalf * ((*un_data)(idx + xp) - (*un_data)(idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(idx) - (*un_data)(idx - xp)));
                double ddx_Thn_dy_vn = eta_n / (dx[1] * dx[0]) *
                                       (thn_iphalf_jmhalf * ((*un_data)(upper_x_idx) - (*un_data)(u_x_idx)) -
                                        thn_imhalf_jmhalf * ((*un_data)(lower_x_idx) - (*un_data)(l_x_idx)));
                double ddy_Thn_dx_vn = -eta_n / (dx[0] * dx[1]) *
                                       ((*thn_data)(idx_c_up) * ((*un_data)(upper_x_idx) - (*un_data)(lower_x_idx)) -
                                        (*thn_data)(idx_c_low) * ((*un_data)(u_x_idx) - (*un_data)(l_x_idx)));

                double drag_n = -xi / nu_n * thn_lower * convertToThs(thn_lower) * ((*un_data)(idx) - (*us_data)(idx));
                double pressure_n = -thn_lower / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*res_un_data)(idx) = (*rhs_un_data)(idx) - (ddy_Thn_dy_un + ddx_Thn_dx_un + ddx_Thn_dy_vn +
                                                             ddy_Thn_dx_vn + drag_n + pressure_n);

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
                    (convertToThs(thn_iphalf_jmhalf) * ((*us_data)(upper_x_idx) - (*us_data)(u_x_idx)) -
                     convertToThs(thn_imhalf_jmhalf) * ((*us_data)(lower_x_idx) - (*us_data)(l_x_idx)));
                double ddy_Ths_dx_vs =
                    -eta_s / (dx[0] * dx[1]) *
                    (convertToThs((*thn_data)(idx_c_up)) * ((*us_data)(upper_x_idx) - (*us_data)(lower_x_idx)) -
                     convertToThs((*thn_data)(idx_c_low)) * ((*us_data)(u_x_idx) - (*us_data)(l_x_idx)));

                double drag_s = -xi / nu_s * thn_lower * convertToThs(thn_lower) * ((*us_data)(idx) - (*un_data)(idx));
                double pressure_s = -convertToThs(thn_lower) / dx[0] * ((*p_data)(idx_c_up) - (*p_data)(idx_c_low));
                (*res_us_data)(idx) = (*rhs_us_data)(idx) - (ddy_Ths_dy_us + ddx_Ths_dx_us + ddx_Ths_dy_vs +
                                                             ddy_Ths_dx_vs + drag_s + pressure_s);
            }
        }
    }

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}

double
convertToThs(double Thn)
{
    return 1 - Thn; // Thn+Ths = 1
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::setToZero(SAMRAIVectorReal<NDIM, double>& vec, int level_num)
{
    const int un_idx = vec.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = vec.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int P_idx = vec.getComponentDescriptorIndex(2);  // pressure
    // d_level_data_ops[level_num]->setToScalar(un_idx, 0.0, /*interior_only*/ false);
    // d_level_data_ops[level_num]->setToScalar(us_idx, 0.0, /*interior_only*/ false);
    // d_level_data_ops[level_num]->setToScalar(p_idx, 0.0, /*interior_only*/ false);
    return;
} // setToZero

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM, double>& solution,
    const SAMRAIVectorReal<NDIM, double>& rhs)
{
    // Setup solution and rhs vectors.
    Pointer<CellVariable<NDIM, double>> solution_var = solution.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double>> rhs_var = rhs.getComponentVariable(0);

    Pointer<CellDataFactory<NDIM, double>> solution_pdat_fac = solution_var->getPatchDataFactory();
    Pointer<CellDataFactory<NDIM, double>> rhs_pdat_fac = rhs_var->getPatchDataFactory();
    d_hierarchy = solution.getPatchHierarchy();

#if !defined(NDEBUG)
    TBOX_ASSERT(solution_var);
    TBOX_ASSERT(rhs_var);
    TBOX_ASSERT(solution_pdat_fac);
    TBOX_ASSERT(rhs_pdat_fac);
#endif

    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::deallocateOperatorState()
{
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
