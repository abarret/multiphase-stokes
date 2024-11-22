/////////////////////////////// INCLUDES /////////////////////////////////////

#include "multiphase/MultiphaseStaggeredStokesBlockFACOperator.h"
#include "multiphase/fd_operators.h"
#include "multiphase/utility_functions.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

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
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/RefinePatchStrategySet.h>

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

#include <LocationIndexRobinBcCoefs.h>

#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////
namespace multiphase
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
static const int SIDEG = 1;

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
static const std::string DATA_COARSEN_TYPE = "NONE";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;
} // namespace
/////////////////////////////// PUBLIC ///////////////////////////////////////

MultiphaseStaggeredStokesBlockFACOperator::MultiphaseStaggeredStokesBlockFACOperator(
    const std::string& object_name,
    const std::string& default_options_prefix,
    const MultiphaseParameters& params)
    : FACPreconditionerStrategy(object_name),
      d_default_un_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_un_bc_coef", Pointer<Database>(nullptr))),
      d_default_us_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_us_bc_coef", Pointer<Database>(nullptr))),
      d_un_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef.get())),
      d_us_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef.get())),
      d_default_thn_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_thn_bc_coef", Pointer<Database>(nullptr))),
      d_thn_bc_coef(d_default_thn_bc_coef.get()),
      d_mask_var(new SideVariable<NDIM, int>(d_object_name + "::mask_var")),
      d_params(params)
{
    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_default_un_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_un_bc_coef.get());
        auto p_default_us_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_us_bc_coef.get());
        p_default_un_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_un_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        p_default_us_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_us_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        auto p_default_thn_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_thn_bc_coef.get());
        p_default_thn_bc_coef->setBoundarySlope(2 * d, 0.0);
        p_default_thn_bc_coef->setBoundarySlope(2 * d + 1, 0.0);
    }
    // Create variables and register them with the variable database.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> d_ctx = var_db->getContext(d_object_name + "::context");

    // State scratch variables: Velocity and pressure.
    // Prepend with d_object_name to ensure there are no conflicts.
    Pointer<SideVariable<NDIM, double>> un_scr_var = new SideVariable<NDIM, double>(d_object_name + "::un_scr");
    Pointer<SideVariable<NDIM, double>> us_scr_var = new SideVariable<NDIM, double>(d_object_name + "::us_scr");

    // Register patch data indices
    if (var_db->checkVariableExists(d_object_name + "::un_scr"))
    {
        un_scr_var = var_db->getVariable(d_object_name + "::un_scr");
        d_un_scr_idx = var_db->mapVariableAndContextToIndex(un_scr_var, d_ctx);
        var_db->removePatchDataIndex(d_un_scr_idx);
    }
    if (var_db->checkVariableExists(d_object_name + "::us_scr"))
    {
        us_scr_var = var_db->getVariable(d_object_name + "::us_scr");
        d_us_scr_idx = var_db->mapVariableAndContextToIndex(us_scr_var, d_ctx);
        var_db->removePatchDataIndex(d_us_scr_idx);
    }
    d_un_scr_idx = var_db->registerVariableAndContext(un_scr_var, d_ctx, IntVector<NDIM>(1));
    d_us_scr_idx = var_db->registerVariableAndContext(us_scr_var, d_ctx, IntVector<NDIM>(1));

    if (var_db->checkVariableExists(d_mask_var->getName()))
    {
        d_mask_var = var_db->getVariable(d_mask_var->getName());
        d_mask_idx = var_db->mapVariableAndContextToIndex(d_mask_var, d_ctx);
        var_db->removePatchDataIndex(d_mask_idx);
    }
    d_mask_idx = var_db->registerVariableAndContext(d_mask_var, d_ctx, IntVector<NDIM>(0));

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error = TimerManager::getManager()->getTimer(
                     "IBTK::MultiphaseStaggeredStokesBlockFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::MultiphaseStaggeredStokesBlockFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::MultiphaseStaggeredStokesBlockFACOperator::computeResidual()"););
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
    RobinBcCoefStrategy<NDIM>* thn_bc_coef)
{
#ifndef NDEBUG
    TBOX_ASSERT(un_bc_coefs.size() == NDIM);
    TBOX_ASSERT(us_bc_coefs.size() == NDIM);
#endif
    // Only replace the boundary conditions that actually exist
    for (int d = 0; d < NDIM; ++d)
    {
        if (un_bc_coefs[d])
            d_un_bc_coefs[d] = un_bc_coefs[d];
        else
            d_un_bc_coefs[d] = d_default_un_bc_coef.get();
        if (us_bc_coefs[d])
            d_us_bc_coefs[d] = us_bc_coefs[d];
        else
            d_us_bc_coefs[d] = d_default_us_bc_coef.get();
    }
    if (thn_bc_coef)
        d_thn_bc_coef = thn_bc_coef;
    else
        d_thn_bc_coef = d_default_thn_bc_coef.get();
}

MultiphaseStaggeredStokesBlockFACOperator::~MultiphaseStaggeredStokesBlockFACOperator()
{
    // Dallocate operator state first
    deallocateOperatorState();
    return;
}

// create another member function to set-up Thn
// Thn is defined in the input file and read in using muParserCartGridFunction
void
MultiphaseStaggeredStokesBlockFACOperator::setThnIdx(int thn_idx)
{
    d_thn_idx = thn_idx;
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::restrictResidual(const SAMRAIVectorReal<NDIM, double>& src,
                                                                    SAMRAIVectorReal<NDIM, double>& dst,
                                                                    const int dst_ln)
{
    // Pull out patch indices.
    const int dst_un_idx = dst.getComponentDescriptorIndex(0);
    const int dst_us_idx = dst.getComponentDescriptorIndex(1);
    std::array<int, 2> dst_idxs = { dst_un_idx, dst_us_idx };

    const int src_un_idx = src.getComponentDescriptorIndex(0);
    const int src_us_idx = src.getComponentDescriptorIndex(1);
    std::array<int, 2> src_idxs = { src_un_idx, src_us_idx };

    // SAMRAI's refine operators will copy data from patch interiors if the patch indices are different.
    // Therefore, I don't think we need to do this
    // TODO: test if this is necessary.
    if (dst_un_idx != src_un_idx)
    {
        HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        level_sc_data_ops.copyData(dst_un_idx, src_un_idx, false /*interior_only*/);
    }
    if (dst_us_idx != src_us_idx)
    {
        HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        level_sc_data_ops.copyData(dst_us_idx, src_us_idx, false /*interior_only*/);
    }

    // Now perform restriction
    performRestriction(dst_idxs, src_idxs, dst_ln);
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::prolongError(const SAMRAIVectorReal<NDIM, double>& src,
                                                                SAMRAIVectorReal<NDIM, double>& dst,
                                                                const int dst_ln)
{
    // Pull out patch indices.
    const int dst_un_idx = dst.getComponentDescriptorIndex(0);
    const int dst_us_idx = dst.getComponentDescriptorIndex(1);
    std::array<int, 2> dst_idxs = { dst_un_idx, dst_us_idx };

    const int src_un_idx = src.getComponentDescriptorIndex(0);
    const int src_us_idx = src.getComponentDescriptorIndex(1);
    std::array<int, 2> src_idxs = { src_un_idx, src_us_idx };

    // Now perform prolongation
    performProlongation(dst_idxs, src_idxs, dst_ln);
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::prolongErrorAndCorrect(const SAMRAIVectorReal<NDIM, double>& src,
                                                                          SAMRAIVectorReal<NDIM, double>& dst,
                                                                          const int dst_ln)
{
    // Pull out patch indices.
    const int dst_un_idx = dst.getComponentDescriptorIndex(0);
    const int dst_us_idx = dst.getComponentDescriptorIndex(1);
    std::array<int, 2> dst_idxs = { d_un_scr_idx, d_us_scr_idx };

    const int src_un_idx = src.getComponentDescriptorIndex(0);
    const int src_us_idx = src.getComponentDescriptorIndex(1);
    std::array<int, 2> src_idxs = { src_un_idx, src_us_idx };

    // I'm not sure why/if we need to do this, but this was done in other implementations...
    // TODO: Test if this is necessary.
    if (dst_un_idx != src_un_idx)
    {
        HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, dst_ln - 1, dst_ln - 1);
        level_sc_data_ops.add(dst_un_idx, dst_un_idx, src_un_idx, false /*interior_only*/);
    }
    if (dst_us_idx != src_us_idx)
    {
        HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, dst_ln - 1, dst_ln - 1);
        level_sc_data_ops.add(dst_us_idx, dst_us_idx, src_us_idx, false /*interior_only*/);
    }

    // Now prolong and correct data
    performProlongation(dst_idxs, src_idxs, dst_ln);
    HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
    level_sc_data_ops.add(dst_un_idx, dst_un_idx, d_un_scr_idx, false /*interior_only*/);
    level_sc_data_ops.add(dst_us_idx, dst_us_idx, d_us_scr_idx, false /*interior_only*/);
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::smoothError(
    SAMRAIVectorReal<NDIM, double>& error,          // Solution
    const SAMRAIVectorReal<NDIM, double>& residual, // RHS - A*error on the entire patch level
    int level_num,
    int num_sweeps,
    bool /*performing_pre_sweeps*/,
    bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

    d_hierarchy = error.getPatchHierarchy();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);

    // Get the vector components. These pull out patch data indices
    const int un_idx = error.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = error.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int thn_idx = d_thn_idx;
    const int f_un_idx = residual.getComponentDescriptorIndex(0); // RHS_Un
    const int f_us_idx = residual.getComponentDescriptorIndex(1); // RHS_Us

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > 0 && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);

            Pointer<SideData<NDIM, double>> un_scr_data = patch->getPatchData(d_un_scr_idx);
            Pointer<SideData<NDIM, double>> us_scr_data = patch->getPatchData(d_us_scr_idx);

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                un_scr_data->getArrayData(axis).copy(un_data->getArrayData(axis),
                                                     d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                     IntVector<NDIM>(0));
                us_scr_data->getArrayData(axis).copy(us_data->getArrayData(axis),
                                                     d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                     IntVector<NDIM>(0));
            }
        }
    }

    IntVector<NDIM> xp(1, 0), yp(0, 1);

    // outer for loop for number of sweeps
    for (int sweep = 0; sweep < num_sweeps; sweep++)
    {   
        if (level_num > 0)
        {
            if (sweep > 0)
            {
                // Copy coarse-fine interface ghost cell values which are cached in the scratch data.
                int patch_counter = 0;
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());

                    Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
                    Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);

                    Pointer<SideData<NDIM, double>> un_scr_data = patch->getPatchData(d_un_scr_idx);
                    Pointer<SideData<NDIM, double>> us_scr_data = patch->getPatchData(d_us_scr_idx);

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        un_data->getArrayData(axis).copy(un_scr_data->getArrayData(axis),
                                                            d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                            IntVector<NDIM>(0));
                        us_data->getArrayData(axis).copy(us_scr_data->getArrayData(axis),
                                                            d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                            IntVector<NDIM>(0));
                    }

                }
            }
            // Fill in ghost cells. We only want to use values on our current level to fill in ghost cells.
            performGhostFilling({ un_idx, us_idx }, level_num);

            // Compute the normal extension of the solution at coarse fine interfaces if we are not on the coarsest
            // level
            // TODO: 0 here is coarsest level number, we should set this to a variable.
            d_sc_bdry_op->setPatchDataIndices({ un_idx, us_idx });
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                // TODO: 1 here is the ghost width to fill, this should be variable.
                const IntVector<NDIM>& gcw_to_fill = 1;
                d_sc_bdry_op->computeNormalExtension(*patch, ratio, gcw_to_fill);
            }
        }
        else if (sweep > 0)
        {
                performGhostFilling({ un_idx, us_idx }, level_num);
        }

        // loop through all patches on this level
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<SideData<NDIM, double>> f_un_data = patch->getPatchData(f_un_idx);
            Pointer<SideData<NDIM, double>> f_us_data = patch->getPatchData(f_us_idx);
            Pointer<SideData<NDIM, int>> mask_data = patch->getPatchData(d_mask_idx);

            // Enforce any Dirichlet boundary conditions.
            if (d_bc_un_helper->patchTouchesDirichletBoundary(patch))
            {
                d_bc_un_helper->copyDataAtDirichletBoundaries(un_data, f_un_data, patch);
            }
            if (d_bc_us_helper->patchTouchesDirichletBoundary(patch))
            {
                d_bc_us_helper->copyDataAtDirichletBoundaries(us_data, f_us_data, patch);
            }

            if (d_params.isVariableDrag())
            {
                // Place var xi implementation here
                Pointer<SideData<NDIM, double>> xi_data = patch->getPatchData(d_params.xi_idx);

            }// end variable drag
            else // use constant xi coefficient
            {
                // Gauss-Seidel iterations with successive under relaxation
                for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& idx = si();           // axis = 0, (i-1/2,j)

                    SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
                    SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
                    SAMRAI::pdat::SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
                    SAMRAI::pdat::SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
                    SAMRAI::pdat::SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
                    SAMRAI::pdat::SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)

                    // thn at sides
                    double thn_lower_x = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up));     // thn(i-1/2,j)

                    // thn at corners
                    double thn_imh_jph = 0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + yp) +
                                                    (*thn_data)(idx_c_low + yp)); // thn(i-1/2,j+1/2)
                    double thn_imh_jmh = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up - yp) +
                                                    (*thn_data)(idx_c_low - yp)); // thn(i-1/2,j-1/2)
                    double thn_iph_jph = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + xp) + (*thn_data)(idx_c_up + yp) +
                                                    (*thn_data)(idx_c_up + xp + yp)); // thn(i+1/2,j+1/2)
                    double thn_iph_jmh = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + xp) + (*thn_data)(idx_c_up - yp) +
                                                    (*thn_data)(idx_c_up + xp - yp)); // thn(i+1/2,j-1/2)

                    double dx_dx = (dx[0] * dx[0]);
                    double dy_dy = (dx[1] * dx[1]);
                    double dx_dy = (dx[0] * dx[1]);

                    // Solve for network velocities
                    double un_imhalf_j = (*f_un_data)(lower_x_idx);   
                    un_imhalf_j -= d_D * d_params.eta_n * ((thn_imh_jph-(*thn_data)(idx))/(dx_dy) * (*un_data)(upper_y_idx) 
                                    + ((*thn_data)(idx-xp))/(dx_dx) * (*un_data)(lower_x_idx - xp) 
                                    + ((*thn_data)(idx)/dx_dx * (*un_data)(upper_x_idx)) 
                                    + (thn_imh_jph/dy_dy * (*un_data)(lower_x_idx + yp)) 
                                    + (thn_imh_jmh/dy_dy * (*un_data)(lower_x_idx - yp))
                                    + ((*thn_data)(idx - xp) - thn_imh_jph)/(dx_dy) * (*un_data)(upper_y_idx - xp)
                                    + (thn_imh_jmh - (*thn_data)(idx - xp))/(dx_dy) * (*un_data)(lower_y_idx - xp)
                                    + ((*thn_data)(idx)-thn_imh_jmh)/(dx_dy) * (*un_data)(lower_y_idx)) 
                                    - d_D * d_params.xi * (*us_data)(lower_x_idx) * thn_lower_x * convertToThs(thn_lower_x);
                    un_imhalf_j /= d_D * d_params.eta_n *
                                       (-((*thn_data)(idx) + (*thn_data)(idx - xp)) / (dx_dx) -
                                        (thn_imh_jph + thn_imh_jmh) / (dy_dy)) +
                                   d_D * d_params.xi * thn_lower_x * convertToThs(thn_lower_x) + d_C * thn_lower_x;
                    (*un_data)(lower_x_idx) = (1.0 - d_w) * (*un_data)(lower_x_idx) + d_w * un_imhalf_j;

                    // Solve for solvent velocities
                    double us_imhalf_j = (*f_us_data)(lower_x_idx);   
                    us_imhalf_j -= d_D * d_params.eta_s * ((convertToThs(thn_imh_jph)-convertToThs((*thn_data)(idx)))/(dx_dy) * (*us_data)(upper_y_idx) 
                                + (convertToThs((*thn_data)(idx-xp)))/(dx_dx) * (*us_data)(lower_x_idx - xp) 
                                + (convertToThs((*thn_data)(idx))/dx_dx * (*us_data)(upper_x_idx)) 
                                + (convertToThs(thn_imh_jph)/dy_dy * (*us_data)(lower_x_idx + yp)) 
                                + (convertToThs(thn_imh_jmh)/dy_dy * (*us_data)(lower_x_idx - yp))
                                + (convertToThs((*thn_data)(idx - xp)) - convertToThs(thn_imh_jph))/(dx_dy) * (*us_data)(upper_y_idx - xp)
                                + (convertToThs(thn_imh_jmh) - convertToThs((*thn_data)(idx - xp)))/(dx_dy) * (*us_data)(lower_y_idx - xp)
                                + (convertToThs((*thn_data)(idx))-convertToThs(thn_imh_jmh))/(dx_dy) * (*us_data)(lower_y_idx)) 
                                - d_D * d_params.xi * (*un_data)(lower_x_idx) * thn_lower_x * convertToThs(thn_lower_x);
                    us_imhalf_j /=
                        d_D * d_params.eta_s *
                            (-(convertToThs((*thn_data)(idx)) + convertToThs((*thn_data)(idx - xp))) / (dx_dx) -
                             (convertToThs(thn_imh_jph) + convertToThs(thn_imh_jmh)) / (dy_dy)) +
                        d_D * d_params.xi * thn_lower_x * convertToThs(thn_lower_x) + d_C * convertToThs(thn_lower_x);
                    (*us_data)(lower_x_idx) = (1.0 - d_w) * (*us_data)(lower_x_idx) + d_w * us_imhalf_j;
                } // side-centers in x-dir

                for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // axis = 1, (i,j-1/2)

                    SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i,j-1)
                    SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
                    SAMRAI::pdat::SideIndex<NDIM> lower_x_idx(idx_c_up, 0, 0); // (i-1/2,j)
                    SAMRAI::pdat::SideIndex<NDIM> upper_x_idx(idx_c_up, 0, 1); // (i+1/2,j)
                    SAMRAI::pdat::SideIndex<NDIM> lower_y_idx(idx_c_up, 1, 0); // (i,j-1/2)
                    SAMRAI::pdat::SideIndex<NDIM> upper_y_idx(idx_c_up, 1, 1); // (i,j+1/2)

                    // thn at sides
                    double thn_lower_y = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i,j-1/2)

                    // thn at corners
                    double thn_imh_jmh = 0.25 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up) + (*thn_data)(idx_c_up - xp) +
                                                    (*thn_data)(idx_c_low - xp)); // thn(i-1/2,j-1/2)
                    double thn_iph_jmh = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_low) + (*thn_data)(idx_c_up + xp) +
                                                    (*thn_data)(idx_c_low + xp)); // thn(i+1/2,j-1/2)
                    double thn_imh_jph = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_up - xp) + (*thn_data)(idx_c_up + yp) +
                                                    (*thn_data)(idx_c_up - xp + yp)); // thn(i-1/2,j+1/2)
                    double thn_iph_jph = 0.25 * ((*thn_data)(idx_c_up) + (*thn_data)(idx_c_up + xp) + (*thn_data)(idx_c_up + yp) +
                                                    (*thn_data)(idx_c_up + xp + yp)); // thn(i+1/2,j+1/2)

                    double dx_dx = (dx[0] * dx[0]);
                    double dy_dy = (dx[1] * dx[1]);
                    double dx_dy = (dx[0] * dx[1]);
                    
                    double un_i_jmhalf = (*f_un_data)(lower_y_idx);
                    un_i_jmhalf -= d_D * d_params.eta_n * ((thn_imh_jmh/dx_dx) * (*un_data)(lower_y_idx - xp)
                                    + (thn_iph_jmh/dx_dx) * (*un_data)(lower_y_idx + xp)
                                    + ((*thn_data)(idx)/dy_dy) * (*un_data)(upper_y_idx)
                                    + ((*thn_data)(idx - yp)/dy_dy) * (*un_data)(lower_y_idx - yp)
                                    + (thn_imh_jmh - (*thn_data)(idx - yp))/(dx_dy) * (*un_data)(lower_x_idx - yp)
                                    + ((*thn_data)(idx - yp) - thn_iph_jmh)/(dx_dy) * (*un_data)(upper_x_idx - yp)
                                    + ((*thn_data)(idx) - thn_imh_jmh)/(dx_dy) * (*un_data)(lower_x_idx)
                                    + (thn_iph_jmh - (*thn_data)(idx))/(dx_dy) * (*un_data)(upper_x_idx)) 
                                    - d_D * d_params.xi * (*us_data)(lower_y_idx)  * thn_lower_y * convertToThs(thn_lower_y);
                    un_i_jmhalf /= d_D * d_params.eta_n *
                                    (-((*thn_data)(idx) + (*thn_data)(idx - yp)) / (dy_dy) -
                                        (thn_iph_jmh + thn_imh_jmh) / (dx_dx)) +
                                d_D * d_params.xi * thn_lower_y * convertToThs(thn_lower_y) + d_C * thn_lower_y;
                    (*un_data)(lower_y_idx) = (1.0 - d_w) * (*un_data)(lower_y_idx) + d_w * un_i_jmhalf;

                    double us_i_jmhalf = (*f_us_data)(lower_y_idx);
                    us_i_jmhalf -= d_D * d_params.eta_s * ((convertToThs(thn_imh_jmh)/dx_dx) * (*us_data)(lower_y_idx - xp)
                                + (convertToThs(thn_iph_jmh)/dx_dx) * (*us_data)(lower_y_idx + xp)
                                + (convertToThs((*thn_data)(idx))/dy_dy) * (*us_data)(upper_y_idx)
                                + (convertToThs((*thn_data)(idx - yp))/dy_dy) * (*us_data)(lower_y_idx - yp)
                                + (convertToThs(thn_imh_jmh) - convertToThs((*thn_data)(idx - yp)))/(dx_dy) * (*us_data)(lower_x_idx - yp)
                                + (convertToThs((*thn_data)(idx - yp)) - convertToThs(thn_iph_jmh))/(dx_dy) * (*us_data)(upper_x_idx - yp)
                                + (convertToThs((*thn_data)(idx)) - convertToThs(thn_imh_jmh))/(dx_dy) * (*us_data)(lower_x_idx)
                                + (convertToThs(thn_iph_jmh) - convertToThs((*thn_data)(idx)))/(dx_dy) * (*us_data)(upper_x_idx)) 
                                - d_D * d_params.xi * (*un_data)(lower_y_idx) * thn_lower_y * convertToThs(thn_lower_y);
                    us_i_jmhalf /=
                        d_D * d_params.eta_s *
                            (-(convertToThs((*thn_data)(idx)) + convertToThs((*thn_data)(idx - yp))) / (dy_dy) -
                            (convertToThs(thn_iph_jmh) + convertToThs(thn_imh_jmh)) / (dx_dx)) +
                        d_D * d_params.xi * thn_lower_y * convertToThs(thn_lower_y) + d_C * convertToThs(thn_lower_y);
                    (*us_data)(lower_y_idx) = (1.0 - d_w) * (*us_data)(lower_y_idx) + d_w * us_i_jmhalf;
                } // side-centers in y-dir
            } // end constant xi coefficient
        } // patches
    }// sweeps
    performGhostFilling({ un_idx, us_idx }, level_num); // Synchronization at patch boundaries
   
    IBTK_TIMER_STOP(t_smooth_error);
    return;
}

bool
MultiphaseStaggeredStokesBlockFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                                      const SAMRAIVectorReal<NDIM, double>& residual,
                                                                      int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);

    smoothError(error, residual, coarsest_ln, 10, false, false);

    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
} // solveCoarsestLevel

void
MultiphaseStaggeredStokesBlockFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                                   const SAMRAIVectorReal<NDIM, double>& solution,
                                                                   const SAMRAIVectorReal<NDIM, double>& rhs,
                                                                   int coarsest_level_num,
                                                                   int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    // Get the vector components. These pull out patch data indices
    const int un_idx = solution.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = solution.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int rhs_un_idx = rhs.getComponentDescriptorIndex(0);  // RHS Un
    const int rhs_us_idx = rhs.getComponentDescriptorIndex(1);  // RHS Us
    const int res_un_idx = residual.getComponentDescriptorIndex(0);
    const int res_us_idx = residual.getComponentDescriptorIndex(1);
    const int thn_idx = d_thn_idx;

    d_un_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    d_us_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
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
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(transaction_comps, d_hierarchy, coarsest_level_num, finest_level_num);
    d_hier_bdry_fill->setHomogeneousBc(true);
    d_hier_bdry_fill->fillData(d_solution_time); // Fills in all of the ghost cells

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
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> rhs_un_data =
                patch->getPatchData(rhs_un_idx); // result of applying operator (eqn 1)
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<SideData<NDIM, double>> rhs_us_data =
                patch->getPatchData(rhs_us_idx); // result of applying operator (eqn 2)
            Pointer<SideData<NDIM, double>> res_un_data = patch->getPatchData(res_un_idx);
            Pointer<SideData<NDIM, double>> res_us_data = patch->getPatchData(res_us_idx);
            IntVector<NDIM> xp(1, 0), yp(0, 1);

            // if (d_params.isVariableDrag())
            //     //accumulateMomentumWithoutPressureOnPatchVariableDrag(
            //         //patch, res_un_idx, res_us_idx, un_idx, us_idx, thn_idx, d_params, d_C, d_D, d_D);
            // else
            accumulateMomentumWithoutPressureOnPatchConstantCoefficient(
                patch, res_un_idx, res_us_idx, un_idx, us_idx, thn_idx, d_params, d_C, d_D);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SideIndex<NDIM>& idx = si();
                    (*res_un_data)(idx) = (*rhs_un_data)(idx) - (*res_un_data)(idx);
                    (*res_us_data)(idx) = (*rhs_us_data)(idx) - (*res_us_data)(idx);
                }
            }

            if (d_bc_un_helper->patchTouchesDirichletBoundary(patch) ||
                d_bc_us_helper->patchTouchesDirichletBoundary(patch))
            {
                Pointer<SideData<NDIM, int>> mask_data = patch->getPatchData(d_mask_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        if ((*mask_data)(idx) == 1)
                        {
                            (*res_un_data)(idx) = (*rhs_un_data)(idx) - (*un_data)(idx);
                            (*res_us_data)(idx) = (*rhs_us_data)(idx) - (*us_data)(idx);
                        }
                    }
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::setToZero(SAMRAIVectorReal<NDIM, double>& vec, int level_num)
{
    const int un_idx = vec.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = vec.getComponentDescriptorIndex(1); // solvent velocity, Us
    d_hierarchy = vec.getPatchHierarchy();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
        Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
        un_data->fillAll(0.0);
        us_data->fillAll(0.0);
    }
    return;
} // setToZero

void
MultiphaseStaggeredStokesBlockFACOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& sol,
                                                                           const SAMRAIVectorReal<NDIM, double>& rhs)
{
    if (d_is_initialized) deallocateOperatorState();
    // Setup solution and rhs vectors.
    Pointer<SideVariable<NDIM, double>> un_sol_var = sol.getComponentVariable(0);
    Pointer<SideVariable<NDIM, double>> us_sol_var = sol.getComponentVariable(1);

    Pointer<SideVariable<NDIM, double>> un_rhs_var = rhs.getComponentVariable(0);
    Pointer<SideVariable<NDIM, double>> us_rhs_var = rhs.getComponentVariable(1);

    d_hierarchy = sol.getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Rudimentary error checking.
#if !defined(NDEBUG)
    TBOX_ASSERT(un_sol_var && us_sol_var);
    TBOX_ASSERT(un_rhs_var && us_rhs_var);
    TBOX_ASSERT(d_hierarchy == rhs.getPatchHierarchy());
#endif

    // Allocate scratch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_un_scr_idx, d_new_time);
        level->allocatePatchData(d_us_scr_idx, d_new_time);
    }

    // Set up physical boundary condition helpers
    d_bc_un_helper = new StaggeredPhysicalBoundaryHelper();
    d_bc_un_helper->cacheBcCoefData(d_un_bc_coefs, d_solution_time, d_hierarchy);
    d_bc_us_helper = new StaggeredPhysicalBoundaryHelper();
    d_bc_us_helper->cacheBcCoefData(d_us_bc_coefs, d_solution_time, d_hierarchy);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_mask_idx)) level->allocatePatchData(d_mask_idx);
    }
    // Note that this assumes that we are using dirichlet conditions at the SAME locations for network and solvent.
    // To remove this restriction, we are going to want to do something more complicated than
    // StaggeredPhysicalBoundaryHelper.
    d_bc_un_helper->setupMaskingFunction(d_mask_idx);

    // Set up boundary condition operator
    d_un_bc_op = new CartSideRobinPhysBdryOp(d_un_scr_idx, d_un_bc_coefs, false);
    d_us_bc_op = new CartSideRobinPhysBdryOp(d_us_scr_idx, d_us_bc_coefs, false);
    std::vector<RefinePatchStrategy<NDIM>*> bc_op_ptrs(2);
    bc_op_ptrs[0] = d_un_bc_op;
    bc_op_ptrs[1] = d_us_bc_op;
    d_vel_P_bc_op = std::make_unique<RefinePatchStrategySet>(bc_op_ptrs.begin(), bc_op_ptrs.end(), false);

    // Cache prolongation operators. Creating refinement schedules can be expensive for hierarchies with many levels. We
    // create the schedules here, and SAMRAI will determine whether they need to be regenerated whenever we switch patch
    // indices.
    // TODO: Set the refine type via input database or setter.
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    d_un_prolong_op = grid_geom->lookupRefineOperator(un_sol_var, "CONSERVATIVE_LINEAR_REFINE");
    d_us_prolong_op = grid_geom->lookupRefineOperator(us_sol_var, "CONSERVATIVE_LINEAR_REFINE");
    const int un_sol_idx = sol.getComponentDescriptorIndex(0);
    const int us_sol_idx = sol.getComponentDescriptorIndex(1);

    d_prolong_alg = new RefineAlgorithm<NDIM>();
    d_prolong_alg->registerRefine(d_un_scr_idx, un_sol_idx, d_un_scr_idx, d_un_prolong_op);
    d_prolong_alg->registerRefine(d_us_scr_idx, us_sol_idx, d_us_scr_idx, d_us_prolong_op);

    d_prolong_scheds.resize(finest_ln - coarsest_ln + 1);
    // Note start from zero because you can't prolong to the coarsest level.
    for (int dst_ln = coarsest_ln + 1; dst_ln <= finest_ln; ++dst_ln)
    {
        d_prolong_scheds[dst_ln] = d_prolong_alg->createSchedule(d_hierarchy->getPatchLevel(dst_ln),
                                                                 Pointer<PatchLevel<NDIM>>(),
                                                                 dst_ln - 1,
                                                                 d_hierarchy,
                                                                 d_vel_P_bc_op.get());
    }

    // Cache restriction operators.
    // TODO: Set the coarsen type via input database or setter.
    d_un_restrict_op = grid_geom->lookupCoarsenOperator(un_sol_var, "CONSERVATIVE_COARSEN");
    d_us_restrict_op = grid_geom->lookupCoarsenOperator(us_sol_var, "CONSERVATIVE_COARSEN");

    d_restrict_alg = new CoarsenAlgorithm<NDIM>();
    d_restrict_alg->registerCoarsen(d_un_scr_idx, un_sol_idx, d_un_restrict_op);
    d_restrict_alg->registerCoarsen(d_us_scr_idx, us_sol_idx, d_us_restrict_op);

    d_restrict_scheds.resize(finest_ln - coarsest_ln);
    // Note don't create one for finest_ln because you can't restrict to the finest level.
    for (int dst_ln = coarsest_ln; dst_ln < finest_ln; ++dst_ln)
    {
        d_restrict_scheds[dst_ln] =
            d_restrict_alg->createSchedule(d_hierarchy->getPatchLevel(dst_ln), d_hierarchy->getPatchLevel(dst_ln + 1));
    }

    // Create operators for only filling ghost cels.
    d_ghostfill_no_restrict_alg = new RefineAlgorithm<NDIM>();
    d_ghostfill_no_restrict_alg->registerRefine(
        d_un_scr_idx, un_sol_idx, d_un_scr_idx, Pointer<RefineOperator<NDIM>>());
    d_ghostfill_no_restrict_alg->registerRefine(
        d_us_scr_idx, us_sol_idx, d_us_scr_idx, Pointer<RefineOperator<NDIM>>());
    d_ghostfill_no_restrict_scheds.resize(finest_ln - coarsest_ln + 1);
    for (int dst_ln = coarsest_ln; dst_ln <= finest_ln; ++dst_ln)
    {
        // We only want to fill in ghost cells from the current level
        // TODO: the second argument here should fill in physical boundary conditions. This only works for periodic
        // conditions.
        d_ghostfill_no_restrict_scheds[dst_ln] =
            d_ghostfill_no_restrict_alg->createSchedule(d_hierarchy->getPatchLevel(dst_ln), d_vel_P_bc_op.get());
    }

    // Coarse-fine boundary operators
    d_sc_bdry_op = new CartSideDoubleQuadraticCFInterpolation();
    d_sc_bdry_op->setConsistentInterpolationScheme(false);
    d_sc_bdry_op->setPatchHierarchy(d_hierarchy);

    // Get overlap information for setting patch boundary conditions.
    d_patch_side_bc_box_overlap.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, 1);

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxList<NDIM>(ghost_box);
            d_patch_cell_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }

    d_is_initialized = true;
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::deallocateOperatorState()
{
    // Delete the operators, schedules, and algorithms
    d_un_restrict_op.setNull();
    d_us_restrict_op.setNull();
    d_restrict_scheds.resize(0);
    d_restrict_alg.setNull();

    d_un_prolong_op.setNull();
    d_us_prolong_op.setNull();
    d_prolong_scheds.resize(0);
    d_prolong_alg.setNull();

    d_ghostfill_no_restrict_scheds.resize(0);
    d_ghostfill_no_restrict_alg.setNull();

    // Deallocate any data associated with the operator
    int coarsest_ln = 0;
    int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_un_scr_idx);
        level->deallocatePatchData(d_us_scr_idx);
        level->deallocatePatchData(d_mask_idx);
    }

    d_is_initialized = false;
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::setUnderRelaxationParamater(const double w)
{
    d_w = w;
}

void
MultiphaseStaggeredStokesBlockFACOperator::setCandDCoefficients(const double C, const double D)
{
    d_C = C;
    d_D = D;
}

void
MultiphaseStaggeredStokesBlockFACOperator::performProlongation(const std::array<int, 2>& dst_idxs,
                                                                       const std::array<int, 2>& src_idxs,
                                                                       const int dst_ln)
{
    d_un_bc_op->setPatchDataIndex(dst_idxs[0]);
    d_un_bc_op->setHomogeneousBc(true);

    d_us_bc_op->setPatchDataIndex(dst_idxs[1]);
    d_us_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(dst_idxs[0], src_idxs[0], dst_idxs[0], d_un_prolong_op);
    refine_alg.registerRefine(dst_idxs[1], src_idxs[1], dst_idxs[1], d_us_prolong_op);

    refine_alg.resetSchedule(d_prolong_scheds[dst_ln]);
    d_prolong_scheds[dst_ln]->fillData(d_new_time);
    d_prolong_alg->resetSchedule(d_prolong_scheds[dst_ln]);
    return;
}

void
MultiphaseStaggeredStokesBlockFACOperator::performRestriction(const std::array<int, 2>& dst_idxs,
                                                                      const std::array<int, 2>& src_idxs,
                                                                      const int dst_ln)
{
    CoarsenAlgorithm<NDIM> coarsen_alg;
    coarsen_alg.registerCoarsen(dst_idxs[0], src_idxs[0], d_un_restrict_op);
    coarsen_alg.registerCoarsen(dst_idxs[1], src_idxs[1], d_us_restrict_op);

    coarsen_alg.resetSchedule(d_restrict_scheds[dst_ln]);
    d_restrict_scheds[dst_ln]->coarsenData();
    d_restrict_alg->resetSchedule(d_restrict_scheds[dst_ln]);
}

void
MultiphaseStaggeredStokesBlockFACOperator::performGhostFilling(const std::array<int, 2>& dst_idxs,
                                                                       const int dst_ln)
{
    d_un_bc_op->setPatchDataIndex(dst_idxs[0]);
    d_un_bc_op->setHomogeneousBc(true);

    d_us_bc_op->setPatchDataIndex(dst_idxs[1]);
    d_us_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(dst_idxs[0], dst_idxs[0], dst_idxs[0], Pointer<RefineOperator<NDIM>>());
    refine_alg.registerRefine(dst_idxs[1], dst_idxs[1], dst_idxs[1], Pointer<RefineOperator<NDIM>>());

    refine_alg.resetSchedule(d_ghostfill_no_restrict_scheds[dst_ln]);
    d_ghostfill_no_restrict_scheds[dst_ln]->fillData(d_new_time);
    d_ghostfill_no_restrict_alg->resetSchedule(d_ghostfill_no_restrict_scheds[dst_ln]);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
