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

#include "ibtk/VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.h"
#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/CartCellDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"

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

#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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

VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator(const std::string& object_name,
                                                                         const Pointer<Database> input_db,
                                                                         const std::string& default_options_prefix)
    : PoissonFACPreconditionerStrategy(
          object_name,
          new CellVariable<NDIM, double>(object_name + "::cell_scratch", DEFAULT_DATA_DEPTH),
          CELLG,
          input_db,
          default_options_prefix),
      d_level_solver_default_options_prefix(default_options_prefix + "level_")
{
    // Set some default values.
    d_prolongation_method = "LINEAR_REFINE";
    d_restriction_method = "CONSERVATIVE_COARSEN";
    d_level_solver_db = new MemoryDatabase(object_name + "::level_solver_db");
    d_coarse_solver_type = CCPoissonSolverManager::PETSC_LEVEL_SOLVER;
    d_coarse_solver_default_options_prefix = default_options_prefix + "level_0_";
    d_coarse_solver_max_iterations = 1;
    d_coarse_solver_db = new MemoryDatabase(object_name + "::coarse_solver_db");

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("smoother_type")) d_smoother_type = input_db->getString("smoother_type");
        if (input_db->keyExists("prolongation_method"))
            d_prolongation_method = input_db->getString("prolongation_method");
        if (input_db->keyExists("restriction_method")) d_restriction_method = input_db->getString("restriction_method");
        if (input_db->keyExists("level_solver_type")) d_level_solver_type = input_db->getString("level_solver_type");
        if (input_db->keyExists("level_solver_rel_residual_tol"))
            d_level_solver_rel_residual_tol = input_db->getDouble("level_solver_rel_residual_tol");
        if (input_db->keyExists("level_solver_abs_residual_tol"))
            d_level_solver_abs_residual_tol = input_db->getDouble("level_solver_abs_residual_tol");
        if (input_db->keyExists("level_solver_max_iterations"))
            d_level_solver_max_iterations = input_db->getInteger("level_solver_max_iterations");
        if (input_db->isDatabase("level_solver_db"))
        {
            d_level_solver_db = input_db->getDatabase("level_solver_db");
        }
        if (input_db->keyExists("coarse_solver_type")) d_coarse_solver_type = input_db->getString("coarse_solver_type");
        if (input_db->keyExists("coarse_solver_rel_residual_tol"))
            d_coarse_solver_rel_residual_tol = input_db->getDouble("coarse_solver_rel_residual_tol");
        if (input_db->keyExists("coarse_solver_abs_residual_tol"))
            d_coarse_solver_abs_residual_tol = input_db->getDouble("coarse_solver_abs_residual_tol");
        if (input_db->keyExists("coarse_solver_max_iterations"))
            d_coarse_solver_max_iterations = input_db->getInteger("coarse_solver_max_iterations");
        if (input_db->isDatabase("coarse_solver_db"))
        {
            d_coarse_solver_db = input_db->getDatabase("coarse_solver_db");
        }
    }

    // Configure the coarse level solver.
    setCoarseSolverType(d_coarse_solver_type);

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error =
                     TimerManager::getManager()->getTimer("IBTK::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::computeResidual()"););
    return;
}

VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::~VCTwoFluidStaggeredStokesBoxRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::setSmootherType(const std::string& level_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherType():\n"
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (d_level_solver_type != level_solver_type)
    {
        d_level_solvers.clear();
    }
    d_level_solver_type = level_solver_type;
    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::setCoarseSolverType(const std::string& coarse_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setCoarseSolverType():\n"
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (d_coarse_solver_type != coarse_solver_type) d_coarse_solver.setNull();
    d_coarse_solver_type = coarse_solver_type;
    if (!d_coarse_solver)
    {
        d_coarse_solver = CCPoissonSolverManager::getManager()->allocateSolver(d_coarse_solver_type,
                                                                               d_object_name + "::coarse_solver",
                                                                               d_coarse_solver_db,
                                                                               d_coarse_solver_default_options_prefix);
    }
    return;
}

// create another member function to set-up Thn
// Thn is defined in the input file and read in using muParserCartGridFunction
void
VCTwoFluidStaggeredStokesOperator::setThnIdx(int thn_idx)
{
    d_thn_idx = thn_idx;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error, // Solution
                                                 const SAMRAIVectorReal<NDIM, double>& residual,      // RHS
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

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);

    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
        const double* const xlow = pgeom->getXLower(); // {xlow[0], xlow[1]} -> physical location of bottom left of box.
        const hier::Index<NDIM>& idx_low = patch->getBox().lower();
        Pointer<CellData<NDIM, double> > thn_data = patch->getPatchData(thn_idx);
        Pointer<SideData<NDIM, double> > un_data = patch->getPatchData(un_idx);
        Pointer<SideData<NDIM, double> > us_data = patch->getPatchData(us_idx);
        Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(P_idx);
        IntVector<NDIM> xp(1, 0), yp(0, 1);

        Matrix<double, 9, 9> A_box;       // 9 x 9 Matrix
        Vector9f b;                       // 9 x 1 RHS vector 
        // Vector9f u_vec;                   // 9 x 1 solution vector

        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++) // cell-centers
        {
            const CellIndex<NDIM>& idx = ci();   // (i,j)
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

            // thn at corners
            double thn_imhalf_jphalf = 0.25 * ((*thn_data)(idx-xp) + (*thn_data)(idx) + (*thn_data)(idx+yp) +
                    (*thn_data)(idx-xp+yp)); // thn(i-1/2,j+1/2)
            double thn_imhalf_jmhalf = 0.25 * ((*thn_data)(idx) + (*thn_data)(idx-xp) + (*thn_data)(idx-yp) +
                    (*thn_data)(idx-xp-yp)); // thn(i-1/2,j-1/2)
            double thn_iphalf_jphalf = 0.25 * ((*thn_data)(idx+xp) + (*thn_data)(idx) + (*thn_data)(idx+yp) +
                    (*thn_data)(idx+xp+yp)); // thn(i+1/2,j+1/2)
            double thn_iphalf_jmhalf = 0.25 * ((*thn_data)(idx+xp) + (*thn_data)(idx) + (*thn_data)(idx-yp) +
                    (*thn_data)(idx+xp-yp)); // thn(i+1/2,j-1/2)

            // network at west edge
            A_box(0,0) = eta_n / (dx[0] * dx[0])*(-(*thn_data)(idx)-(*thn_data)(idx-xp)) - eta_n/(dx[1]*dx[1])*(thn_imhalf_jmhalf+thn_imhalf_jphalf)
                            - xi/nu_n * thn_lower_x * convertToThs(thn_lower_x);
            A_box(0,1) = eta_n / (dx[0] * dx[0])*((*thn_data)(idx));
            A_box(0,2) = eta_n / (dx[0] * dx[1])*((*thn_data)(idx)-thn_imhalf_jmhalf);
            A_box(0,3) = eta_n / (dx[0] * dx[1])*(thn_imhalf_jphalf-(*thn_data)(idx));
            A_box(0,4) = xi/nu_n * thn_lower_x * convertToThs(thn_lower_x);
            A_box(0,5) = A_box(0,6) = A_box(0,7) = 0.0;
            A_box(0,8) = -thn_lower_x/dx[0];

            A_box(4,4) = eta_s / (dx[0] * dx[0])*(-convertToThs((*thn_data)(idx))-convertToThs((*thn_data)(idx-xp))) - eta_s/(dx[1]*dx[1])*(convertToThs(thn_imhalf_jmhalf)
                        + convertToThs(thn_imhalf_jphalf)) - xi/nu_n * thn_lower_x * convertToThs(thn_lower_x);
            A_box(4,5) = eta_s / (dx[0] * dx[0])*(convertToThs((*thn_data)(idx)));
            A_box(4,6) = eta_s / (dx[0] * dx[1])*(convertToThs((*thn_data)(idx))-convertToThs(thn_imhalf_jmhalf));
            A_box(4,7) = eta_s / (dx[0] * dx[1])*(convertToThs(thn_imhalf_jphalf)-convertToThs((*thn_data)(idx)));
            A_box(4,0) = xi/nu_s * thn_lower_x * convertToThs(thn_lower_x);
            A_box(4,1) = A_box(4,2) = A_box(4,3) = 0.0;
            A_box(4,8) = -convertToThs(thn_lower_x)/dx[0];

            // network at east edge
            A_box(1,0) = eta_n / (dx[0] * dx[0])*(*thn_data)(idx);
            A_box(1,1) = eta_n / (dx[0] * dx[0])*(-(*thn_data)(idx+xp)-(*thn_data)(idx)) - eta_n / (dx[1]*dx[1])*(thn_iphalf_jphalf+thn_iphalf_jmhalf)
                        - xi/nu_n * thn_upper_x * convertToThs(thn_upper_x);
            A_box(1,2) = eta_n / (dx[1] * dx[0])*(thn_iphalf_jmhalf-(*thn_data)(idx));
            A_box(1,3) = eta_n / (dx[1] * dx[0])*((*thn_data)(idx)-thn_iphalf_jphalf)
            A_box(1,4) = A_box(1,6) = A_box(1,7) = 0.0;
            A_box(1,5) = xi/nu_n * thn_upper_x * convertToThs(thn_upper_x);
            A_box(1,8) = thn_upper_x/dx[0];

            A_box(5,4) = eta_s / (dx[0] * dx[0])*convertToThs((*thn_data)(idx));
            A_box(5,5) = eta_s / (dx[0] * dx[0])*(-convertToThs((*thn_data)(idx+xp))-convertToThs((*thn_data)(idx)))
                        - eta_s / (dx[1]*dx[1])* (convertToThs(thn_iphalf_jphalf) + convertToThs(thn_iphalf_jmhalf)) - xi/nu_s * thn_upper_x * convertToThs(thn_upper_x);
            A_box(5,6) = eta_n / (dx[1] * dx[0])*(convertToThs(thn_iphalf_jmhalf)-convertToThs((*thn_data)(idx)));
            A_box(5,7) = eta_n / (dx[1] * dx[0])*(convertToThs((*thn_data)(idx))-convertToThs(thn_iphalf_jphalf));
            A_box(5,0) = A_box(5,2) = A_box(5,3) = 0.0;
            A_box(5,1) = xi/nu_s * thn_upper_x * convertToThs(thn_upper_x);
            A_box(5,8) = convertToThs(thn_upper_x)/dx[0];

            // network at south edge
            A_box(2,0) = eta_n / (dx[0] * dx[1])*((*thn_data)(idx)-thn_imhalf_jmhalf);
            A_box(2,1) = eta_n / (dx[0] * dx[1])*(thn_iphalf_jmhalf-(*thn_data)(idx));
            A_box(2,2) = eta_n / (dx[1] * dx[1])*(-(*thn_data)(idx)-(*thn_data)(idx-yp)) - eta_n/(dx[0]*dx[0])*(thn_iphalf_jmhalf+thn_imhalf_jmhalf)
                        -xi / nu_n*thn_lower_y*convertToThs(thn_lower_y));
            A_box(2,3) = eta_n / (dx[1] * dx[1])*((*thn_data)(idx));
            A_box(2,4) = A_box(2,5) = A_box(2,7) = 0.0;
            A_box(2,6) = xi/nu_n * thn_lower_y * convertToThs(thn_lower_y);
            A_box(2,8) = -thn_lower_y/ dx[1];

            A_box(6,4) = eta_s / (dx[0] * dx[1])*(convertToThs((*thn_data)(idx))-convertToThs(thn_imhalf_jmhalf));
            A_box(6,5) = eta_s / (dx[0] * dx[1])*(convertToThs(thn_iphalf_jmhalf)-convertToThs((*thn_data)(idx)));
            A_box(6,6) = eta_s / (dx[1] * dx[1])*(-convertToThs((*thn_data)(idx))-convertToThs((*thn_data)(idx-yp)))
                        - eta_s / (dx[0] * dx[0])*(convertToThs(thn_iphalf_jmhalf)+convertToThs(thn_imhalf_jmhalf))- xi/nu_n*thn_lower_y*convertToThs(thn_lower_y);
            A_box(6,7) = eta_s / (dx[1] * dx[1])*(convertToThs((*thn_data)(idx)));
            A_box(6,0) = A_box(6,1) = A_box(6,3) = 0.0;
            A_box(6,2) = xi/nu_s * thn_lower_y * convertToThs(thn_lower_y);
            A_box(6,8) = -convertToThs(thn_lower_y)/ dx[1];

            // network at north edge
            A_box(3,0) = eta_n / (dx[0] * dx[1])*(thn_imhalf_jphalf-(*thn_data)(idx));
            A_box(3,1) = eta_n / (dx[0] * dx[1])*((*thn_data)(idx)-thn_iphalf_jphalf);
            A_box(3,2) = eta_n / (dx[1] * dx[1])*((*thn_data)(idx));
            A_box(3,3) = eta_n / (dx[1] * dx[1])*(-(*thn_data)(idx)-(*thn_data)(idx+yp)) - eta_n / (dx[0] * dx[0])*(thn_iphalf_jphalf-thn_imhalf_jphalf)
                        -xi/nu_n*thn_upper_y*convertToThs(thn_upper_y);
            A_box(3,4) = A_box(3,5) = A_box(3,6) = 0.0;
            A_box(3,7) = xi/nu_n * thn_upper_y * convertToThs(thn_upper_y);
            A_box(3,8) = thn_upper_y/ dx[1];

            A_box(7,4) = eta_s / (dx[0] * dx[1])*(convertToThs(thn_imhalf_jphalf)-convertToThs((*thn_data)(idx)));
            A_box(7,5) = eta_s / (dx[0] * dx[1])*(convertToThs((*thn_data)(idx))-convertToThs(thn_iphalf_jphalf));
            A_box(7,6) = eta_s / (dx[1] * dx[1])*(convertToThs((*thn_data)(idx)));
            A_box(7,7) = eta_s / (dx[1] * dx[1])*(-convertToThs((*thn_data)(idx))-convertToThs((*thn_data)(idx+yp)))
                        -  eta_s / (dx[0] * dx[0])*(convertToThs(thn_iphalf_jphalf)-convertToThs(thn_imhalf_jphalf))-xi / nu_s*thn_upper_y*convertToThs(thn_upper_y);
            A_box(7,0) = A_box(7,1) = A_box(7,2) = 0.0;
            A_box(7,3) = xi/nu_s * thn_upper_y * convertToThs(thn_upper_y);
            A_box(7,8) = convertToThs(thn_upper_y)/ dx[1];

            // incompressible constrain term at center
            A_box(8,0) = thn_lower_x / dx[0];
            A_box(8,1) = -thn_upper_x / dx[0];
            A_box(8,2) = thn_lower_y / dx[1];
            A_box(8,3) = -thn_upper_y / dx[1];
            A_box(8,4) = convertToThs(thn_lower_x) / dx[0];
            A_box(8,5) = -convertToThs(thn_upper_x) / dx[0];
            A_box(8,6) = convertToThs(thn_lower_y) / dx[1];
            A_box(8,7) = -convertToThs(thn_upper_y) / dx[1];
            A_box(8,8) = 0.0;

            //set-up RHS vector

            // network at west edge
            b(0) = -thn_lower_x/dx[0]*(*P_data)(idx-xp) - eta_n/(dx[0]*dx[0])*(*thn_data)(idx-xp)*(*un_data)(lower_x_idx-xp) + eta_n/(dx[1]*dx[1])*thn_imhalf_jphalf*(*un_data)(lower_x_idx+yp)-eta_n/(dx[1]*dx[1])*thn_imhalf_jmhalf*(*un_data)(lower_x_idx-yp)...
                + eta_n/(dx[0]*dx[1])*thn_imhalf_jphalf*(*un_data)(upper_y_idx-xp) - eta_n/(dx[0]*dx[1])*thn_imhalf_jmhalf*(*un_data)(lower_y_idx-xp) - eta_n/(dx[0]*dx[1])*(*thn_data)(idx-xp)*((*un_data)(upper_y_idx-xp) -(*un_data)(lower_y_idx-xp));
            
            // network at east edge
            b(1) = thn_upper_x/dx[0]*(*P_data)(idx+xp) - eta_n/(dx[0]*dx[0])*(*thn_data)(idx+xp)*(*un_data)(upper_x_idx+xp)-eta_n/(dx[1]*dx[1])*thn_iphalf_jphalf*(*un_data)(upper_x_idx+yp) - eta_n/(dx[1]*dx[1])*thn_iphalf_jmhalf*(*un_data)(upper_x_idx-yp)...
             - eta_n/(dx[0]*dx[1])*thn_iphalf_jphalf*(*un_data)(upper_y_idx+xp) + eta_n/(dx[0]*dx[1])*thn_iphalf_jmhalf*(*un_data)(lower_y_idx+xp) + eta_n/(dx[0]*dx[1])*(*thn_data)(idx+xp)*((*un_data)(upper_y_idx+xp) -(*un_data)(lower_y_idx+xp));

            // solvent at west edge
            b(4) = -convertToThs(thn_lower_x)/dx[0]*(*P_data)(idx-xp) - eta_s/(dx[0]*dx[0])*convertToThs((*thn_data)(idx-xp))*(*us_data)(lower_x_idx-xp) + eta_s/(dx[1]*dx[1])*convertToThs(thn_imhalf_jphalf)*(*us_data)(lower_x_idx+yp)-eta_s/(dx[1]*dx[1])*convertToThs(thn_imhalf_jmhalf)*(*us_data)(lower_x_idx-yp)...
                + eta_s/(dx[0]*dx[1])*convertToThs(thn_imhalf_jphalf)*(*us_data)(upper_y_idx-xp) - eta_s/(dx[0]*dx[1])*convertToThs(thn_imhalf_jmhalf)*(*us_data)(lower_y_idx-xp) - eta_s/(dx[0]*dx[1])*convertToThs((*thn_data)(idx-xp))*((*us_data)(upper_y_idx-xp) -(*us_data)(lower_y_idx-xp));

            // solvent at east edge
            b(5) = convertToThs(thn_upper_x)/dx[0]*(*P_data)(idx+xp) - eta_s/(dx[0]*dx[0])*convertToThs((*thn_data)(idx+xp))*(*us_data)(upper_x_idx+xp)-eta_s/(dx[1]*dx[1])*convertToThs(thn_iphalf_jphalf)*(*us_data)(upper_x_idx+yp) - eta_s/(dx[1]*dx[1])*convertToThs(thn_iphalf_jmhalf)*(*us_data)(upper_x_idx-yp)...
             - eta_s/(dx[0]*dx[1])*convertToThs(thn_iphalf_jphalf)*(*us_data)(upper_y_idx+xp) + eta_s/(dx[0]*dx[1])*convertToThs(thn_iphalf_jmhalf)*(*us_data)(lower_y_idx+xp) + eta_s/(dx[0]*dx[1])*convertToThs((*thn_data)(idx+xp))*((*us_data)(upper_y_idx+xp) -(*us_data)(lower_y_idx+xp)); 

            // network at south edge
            b(2) = -thn_lower_y/dx[1]*(*P_data)(idx-yp)-eta_n/(dx[1]*dx[1])*(*thn_data)(idx-yp)*(*un_data)(lower_y_idx-yp) - eta_n/(dx[0]*dx[0])*thn_iphalf_jmhalf*(*un_data)(lower_y_idx+xp)-eta_n/(dx[0]*dx[0])*thn_imhalf_jmhalf*(*un_data)(lower_y_idx-xp)...
            + eta_n/(dx[0]*dx[1])*thn_iphalf_jmhalf*(*un_data)(upper_x_idx-yp) - eta_n/(dx[0]*dx[1])*thn_imhalf_jmhalf*(*un_data)(lower_x_idx-yp) - eta_n/(dx[0]*dx[1])*(*thn_data)(idx-yp)*((*un_data)(upper_x_idx-yp)-(*un_data)(lower_x_idx-yp));

            // solvent at south edge
            b(6) = -convertToThs(thn_lower_y)/dx[1]*(*P_data)(idx-yp)-eta_s/(dx[1]*dx[1])*convertToThs((*thn_data)(idx-yp))*(*us_data)(lower_y_idx-yp) - eta_s/(dx[0]*dx[0])*convertToThs(thn_iphalf_jmhalf)*(*us_data)(lower_y_idx+xp)-eta_s/(dx[0]*dx[0])*convertToThs(thn_imhalf_jmhalf)*(*us_data)(lower_y_idx-xp)...
            + eta_s/(dx[0]*dx[1])*convertToThs(thn_iphalf_jmhalf)*(*us_data)(upper_x_idx-yp) - eta_s/(dx[0]*dx[1])*convertToThs(thn_imhalf_jmhalf)*(*us_data)(lower_x_idx-yp) - eta_s/(dx[0]*dx[1])*convertToThs((*thn_data)(idx-yp))*((*us_data)(upper_x_idx-yp)-(*us_data)(lower_x_idx-yp));

            // network at north edge
            b(3) = thn_upper_y/dx[1]*(*P_data)(idx+yp)-eta_n/(dx[1]*dx[1])*(*thn_data)(idx+yp)*(*un_data)(upper_y_idx+yp) - eta_n/(dx[0]*dx[0])*thn_iphalf_jphalf*(*un_data)(upper_y_idx+xp)-eta_n/(dx[0]*dx[0])*thn_imhalf_jphalf*(*un_data)(upper_y_idx-xp)...
            - eta_n/(dx[0]*dx[1])*thn_iphalf_jphalf*(*un_data)(upper_x_idx+yp) + eta_n/(dx[0]*dx[1])*thn_imhalf_jphalf*(*un_data)(lower_x_idx+yp) + eta_n/(dx[0]*dx[1])*(*thn_data)(idx+yp)*((*un_data)(upper_x_idx+yp)-(*un_data)(lower_x_idx+yp));

            // solvent at north edge
            b(7) = convertToThs(thn_upper_y)/dx[1]*(*P_data)(idx+yp)-eta_s/(dx[1]*dx[1])*convertToThs((*thn_data)(idx+yp))*(*us_data)(upper_y_idx+yp) - eta_s/(dx[0]*dx[0])*convertToThs(thn_iphalf_jphalf)*(*us_data)(upper_y_idx+xp)-eta_s/(dx[0]*dx[0])*convertToThs(thn_imhalf_jphalf)*(*us_data)(upper_y_idx-xp)...
            - eta_s/(dx[0]*dx[1])*convertToThs(thn_iphalf_jphalf)*(*us_data)(upper_x_idx+yp) + eta_s/(dx[0]*dx[1])*convertToThs(thn_imhalf_jphalf)*(*us_data)(lower_x_idx+yp) + eta_s/(dx[0]*dx[1])*convertToThs((*thn_data)(idx+yp))*((*us_data)(upper_x_idx+yp)-(*us_data)(lower_x_idx+yp));
        
            b(8) = 0.0;

            error = A_box.lu().solve(b); // solve Ax = b per cell 
        }

    }
    IBTK_TIMER_STOP(t_smooth_error);
    return;
}

bool
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                        const SAMRAIVectorReal<NDIM, double>& residual,
                                                        int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);
#if !defined(NDEBUG)
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
    TBOX_ASSERT(d_coarse_solver);
#endif
    Pointer<SAMRAIVectorReal<NDIM, double> > e_level = getLevelSAMRAIVectorReal(error, coarsest_ln);
    Pointer<SAMRAIVectorReal<NDIM, double> > r_level = getLevelSAMRAIVectorReal(residual, coarsest_ln);
    d_coarse_solver->setSolutionTime(d_solution_time);
    d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
    auto p_coarse_solver = dynamic_cast<LinearSolver*>(d_coarse_solver.getPointer());
    if (p_coarse_solver)
    {
        bool initial_guess_nonzero = true;
        auto p_petsc_solver = dynamic_cast<PETScKrylovLinearSolver*>(p_coarse_solver);
        auto p_petsc_level_solver = dynamic_cast<PETScLevelSolver*>(p_coarse_solver);
        if (p_petsc_solver || p_petsc_level_solver)
        {
            const KSP& petsc_ksp = p_petsc_solver ? p_petsc_solver->getPETScKSP() : p_petsc_level_solver->getPETScKSP();
            KSPType ksp_type;
            KSPGetType(petsc_ksp, &ksp_type);

            if (!std::strcmp(ksp_type, "preonly")) initial_guess_nonzero = false;
        }
        p_coarse_solver->setInitialGuessNonzero(initial_guess_nonzero);
    }

    d_coarse_solver->solveSystem(*e_level, *r_level);
    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                     const SAMRAIVectorReal<NDIM, double>& solution,
                                                     const SAMRAIVectorReal<NDIM, double>& rhs,
                                                     int coarsest_level_num,
                                                     int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    const int res_idx = residual.getComponentDescriptorIndex(0);
    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<CellVariable<NDIM, double> > res_var = residual.getComponentVariable(0);
    const Pointer<CellVariable<NDIM, double> > sol_var = solution.getComponentVariable(0);

    // Fill ghost-cell values.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    Pointer<CellNoCornersFillPattern> fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    InterpolationTransactionComponent transaction_comp(sol_idx,
                                                       DATA_REFINE_TYPE,
                                                       USE_CF_INTERPOLATION,
                                                       DATA_COARSEN_TYPE,
                                                       BDRY_EXTRAP_TYPE,
                                                       CONSISTENT_TYPE_2_BDRY,
                                                       d_bc_coefs,
                                                       fill_pattern);
    if (d_level_bdry_fill_ops[finest_level_num])
    {
        d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponent(transaction_comp);
    }
    else
    {
        d_level_bdry_fill_ops[finest_level_num] = new HierarchyGhostCellInterpolation();
        d_level_bdry_fill_ops[finest_level_num]->initializeOperatorState(
            transaction_comp, d_hierarchy, coarsest_level_num, finest_level_num);
    }
    d_level_bdry_fill_ops[finest_level_num]->setHomogeneousBc(true);
    d_level_bdry_fill_ops[finest_level_num]->fillData(d_solution_time);
    InterpolationTransactionComponent default_transaction_comp(d_solution->getComponentDescriptorIndex(0),
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_bc_coefs,
                                                               fill_pattern);
    d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponent(default_transaction_comp);

    // Compute the residual, r = f - A*u.
    if (!d_level_math_ops[finest_level_num])
    {
        d_level_math_ops[finest_level_num] =
            new HierarchyMathOps(d_object_name + "::hier_math_ops_" + std::to_string(finest_level_num),
                                 d_hierarchy,
                                 coarsest_level_num,
                                 finest_level_num);
    }
    d_level_math_ops[finest_level_num]->laplace(
        res_idx, res_var, d_poisson_spec, sol_idx, sol_var, nullptr, d_solution_time);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_cc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}

double
convertToThs(double Thn)
{
    return 1 - Thn; // Thn+Ths = 1
}

void
PoissonFACPreconditionerStrategy::setToZero(SAMRAIVectorReal<NDIM, double>& vec, int level_num)
{
    const int un_idx = vec.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = vec.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int P_idx = vec.getComponentDescriptorIndex(2);  // pressure
    d_level_data_ops[level_num]->setToScalar(un_idx, 0.0, /*interior_only*/ false);
    d_level_data_ops[level_num]->setToScalar(us_idx, 0.0, /*interior_only*/ false);
    d_level_data_ops[level_num]->setToScalar(p_idx, 0.0, /*interior_only*/ false);
    return;
} // setToZero

/////////////////////////////// PROTECTED ////////////////////////////////////
void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
                                                                        const SAMRAIVectorReal<NDIM, double>& rhs,
                                                                        const int coarsest_reset_ln,
                                                                        const int finest_reset_ln)
{
    // Setup solution and rhs vectors.
    Pointer<CellVariable<NDIM, double> > solution_var = solution.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > rhs_var = rhs.getComponentVariable(0);

    Pointer<CellDataFactory<NDIM, double> > solution_pdat_fac = solution_var->getPatchDataFactory();
    Pointer<CellDataFactory<NDIM, double> > rhs_pdat_fac = rhs_var->getPatchDataFactory();

#if !defined(NDEBUG)
    TBOX_ASSERT(solution_var);
    TBOX_ASSERT(rhs_var);
    TBOX_ASSERT(solution_pdat_fac);
    TBOX_ASSERT(rhs_pdat_fac);
#endif

    if (solution_pdat_fac->getDefaultDepth() != rhs_pdat_fac->getDefaultDepth())
    {
        TBOX_ERROR("VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = " << solution_pdat_fac->getDefaultDepth() << "\n"
                   << "  rhs      data depth = " << rhs_pdat_fac->getDefaultDepth() << std::endl);
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellDataFactory<NDIM, double> > scratch_pdat_fac =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
    scratch_pdat_fac->setDefaultDepth(solution_pdat_fac->getDefaultDepth());

    // Initialize the fine level solvers when needed.
    d_level_solvers.resize(d_finest_ln + 1);
    for (int ln = std::max(0, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PoissonSolver>& level_solver = d_level_solvers[ln];
        if (!level_solver)
        {
            level_solver = CCPoissonSolverManager::getManager()->allocateSolver(d_level_solver_type,
                                                                                d_object_name + "::level_solver",
                                                                                d_level_solver_db,
                                                                                d_level_solver_default_options_prefix +
                                                                                    std::to_string(ln) + "_");
        }

        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setPoissonSpecifications(d_poisson_spec);
        level_solver->setPhysicalBcCoefs(d_bc_coefs);
        level_solver->setMaxIterations(d_level_solver_max_iterations);
        level_solver->setAbsoluteTolerance(d_level_solver_abs_residual_tol);
        level_solver->setRelativeTolerance(d_level_solver_rel_residual_tol);
        level_solver->setHomogeneousBc(true);
        level_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, ln),
                                            *getLevelSAMRAIVectorReal(*d_rhs, ln));
    }

    // Initialize the coarse level solvers when needed.
    if (coarsest_reset_ln == d_coarsest_ln)
    {
        // Note that since the coarse level solver is solving for the error, it
        // must always employ homogeneous boundary conditions.
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setPoissonSpecifications(d_poisson_spec);
        d_coarse_solver->setPhysicalBcCoefs(d_bc_coefs);
        d_coarse_solver->setMaxIterations(d_coarse_solver_max_iterations);
        d_coarse_solver->setAbsoluteTolerance(d_coarse_solver_abs_residual_tol);
        d_coarse_solver->setRelativeTolerance(d_coarse_solver_rel_residual_tol);
        d_coarse_solver->setHomogeneousBc(true);
        d_coarse_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, d_coarsest_ln),
                                               *getLevelSAMRAIVectorReal(*d_rhs, d_coarsest_ln));
    }

    // Setup specialized transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    IBTK_DO_ONCE(geometry->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen()););

    // Setup coarse-fine interface and physical boundary operators.
    d_cf_bdry_op = new CartCellDoubleQuadraticCFInterpolation();
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);
    d_bc_op = new CartCellRobinPhysBdryOp(d_scratch_idx, d_bc_coefs, false);

    // Setup fill pattern spec objects.
    if (d_poisson_spec.dIsConstant())
    {
        d_op_stencil_fill_pattern = new CellNoCornersFillPattern(CELLG, true, false, false);
    }
    else
    {
        d_op_stencil_fill_pattern = nullptr;
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, 1);
            d_patch_bc_box_overlap[ln][patch_counter] = BoxList<NDIM>(ghost_box);
            d_patch_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }
    return;
}

void
VCTwoFluidStaggeredStokesBoxRelaxationFACOperator::deallocateOperatorStateSpecialized(const int /*coarsest_reset_ln*/,
                                                                        const int /*finest_reset_ln*/)
{
    if (!d_is_initialized) return;
    if (!d_in_initialize_operator_state)
    {
        d_patch_bc_box_overlap.clear();
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            if (d_level_solvers[ln]) d_level_solvers[ln]->deallocateSolverState();
        }
        if (d_coarse_solver) d_coarse_solver->deallocateSolverState();
    }
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

////////////////////////////////////////////////////////////////////////////// 
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////