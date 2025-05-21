/////////////////////////////// INCLUDES /////////////////////////////////////

#include "multiphase/MultiphaseStaggeredStokesBoxRelaxationFACOperator.h"
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

// FORTRAN ROUTINES
#define R_B_G_S IBTK_FC_FUNC_(rbgs, RBGS)
#define R_B_G_S_mask IBTK_FC_FUNC_(rbgs_mask, RBGS)
#define R_B_G_S_var_xi IBTK_FC_FUNC_(rbgs_var_xi, RBGS)
#define R_B_G_S_var_xi_mask IBTK_FC_FUNC_(rbgs_mask_var_xi, RBGS)

extern "C"
{
    void R_B_G_S(const double*, // dx
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
                 double* const, // p_data_
                 const int&,    // p_gcw
                 double* const, // f_p_data
                 const int&,    // f_p_gcw
                 double* const, // f_un_data_0
                 double* const, // f_un_data_1
                 const int&,    // f_un_gcw
                 double* const, // f_us_data_0
                 double* const, // f_us_data_1
                 const int&,    // f_us_gcw
                 double* const, // thn_data
                 const int&,    // thn_gcw
                 double* const, // thn_nc_data
                 const int&,
                 double* const, // thn_sc_0
                 double* const, // thn_sc_1
                 const int&,
                 const double&, // eta_n
                 const double&, // eta_s
                 const double&, // lambda_n
                 const double&, // lambda_s
                 const double&, // nu_n
                 const double&, // nu_s
                 const double&, // xi
                 const double&, // w = under relaxation factor
                 const double&, // C in C*u term
                 const double&, // D
                 const int&);   // red_or_black

    void R_B_G_S_mask(const double*, // dx
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
                      double* const, // p_data_
                      const int&,    // p_gcw
                      double* const, // f_p_data
                      const int&,    // f_p_gcw
                      double* const, // f_un_data_0
                      double* const, // f_un_data_1
                      const int&,    // f_un_gcw
                      double* const, // f_us_data_0
                      double* const, // f_us_data_1
                      const int&,    // f_us_gcw
                      double* const, // thn_data
                      const int&,    // thn_gcw
                      double* const, // thn_nc_data
                      const int&,
                      double* const, // thn_sc_0
                      double* const, // thn_sc_1
                      const int&,
                      const double&, // eta_n
                      const double&, // eta_s
                      const double&, // lambda_n
                      const double&, // lambda_s
                      const double&, // nu_n
                      const double&, // nu_s
                      const double&, // xi
                      const double&, // w = under relaxation factor
                      const double&, // C in C*u term
                      const double&, // D
                      const int&,    // red_or_black
                      const int*,    // mask_0
                      const int*,    // mask_1
                      const int&);   // mask_gcw

    void R_B_G_S_var_xi(const double*, // dx
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
                        double* const, // p_data_
                        const int&,    // p_gcw
                        double* const, // f_p_data
                        const int&,    // f_p_gcw
                        double* const, // f_un_data_0
                        double* const, // f_un_data_1
                        const int&,    // f_un_gcw
                        double* const, // f_us_data_0
                        double* const, // f_us_data_1
                        const int&,    // f_us_gcw
                        double* const, // thn_data
                        const int&,    // thn_gcw
                        double* const, // thn_nc_data
                        const int&,
                        double* const, // thn_sc_0
                        double* const, // thn_sc_1
                        const int&,
                        const double&, // eta_n
                        const double&, // eta_s
                        const double&, // lambda_n
                        const double&, // lambda_s
                        double* const, // xi_0
                        double* const, // xi_1
                        const int&,    // xi_gcw
                        const double&, // w = under relaxation factor
                        const double&, // C in C*u term
                        const double&, // D
                        const int&);   // red_or_black

    void R_B_G_S_var_xi_mask(const double*, // dx
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
                             double* const, // p_data_
                             const int&,    // p_gcw
                             double* const, // f_p_data
                             const int&,    // f_p_gcw
                             double* const, // f_un_data_0
                             double* const, // f_un_data_1
                             const int&,    // f_un_gcw
                             double* const, // f_us_data_0
                             double* const, // f_us_data_1
                             const int&,    // f_us_gcw
                             double* const, // thn_data
                             const int&,    // thn_gcw
                             double* const, // thn_nc_data
                             const int&,
                             double* const, // thn_sc_0
                             double* const, // thn_sc_1
                             const int&,
                             const double&, // eta_n
                             const double&, // eta_s
                             const double&, // lambda_n
                             const double&, // lambda_s
                             double* const, // xi_0
                             double* const, // xi_1
                             const int&,    // xi_gcw
                             const double&, // w = under relaxation factor
                             const double&, // C in C*u term
                             const double&, // D
                             const int&,    // red_or_black
                             const int*,    // mask_0
                             const int*,    // mask_1
                             const int&);   // mask_gcw
}
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

MultiphaseStaggeredStokesBoxRelaxationFACOperator::MultiphaseStaggeredStokesBoxRelaxationFACOperator(
    const std::string& object_name,
    const std::string& default_options_prefix,
    const MultiphaseParameters& params,
    const std::unique_ptr<VolumeFractionDataManager>& thn_manager)
    : FACPreconditionerStrategy(object_name),
      d_default_un_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_un_bc_coef", Pointer<Database>(nullptr))),
      d_default_us_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_us_bc_coef", Pointer<Database>(nullptr))),
      d_un_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_un_bc_coef.get())),
      d_us_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_us_bc_coef.get())),
      d_default_P_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(nullptr))),
      d_P_bc_coef(d_default_P_bc_coef.get()),
      d_mask_var(new SideVariable<NDIM, int>(d_object_name + "::mask_var")),
      d_params(params),
      d_thn_manager(thn_manager)
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
        auto p_default_P_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_P_bc_coef.get());
        p_default_P_bc_coef->setBoundarySlope(2 * d, 0.0);
        p_default_P_bc_coef->setBoundarySlope(2 * d + 1, 0.0);
    }
    // Create variables and register them with the variable database.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> d_ctx = var_db->getContext(d_object_name + "::context");

    // State scratch variables: Velocity and pressure.
    // Prepend with d_object_name to ensure there are no conflicts.
    Pointer<SideVariable<NDIM, double>> un_scr_var = new SideVariable<NDIM, double>(d_object_name + "::un_scr");
    Pointer<SideVariable<NDIM, double>> us_scr_var = new SideVariable<NDIM, double>(d_object_name + "::us_scr");
    Pointer<CellVariable<NDIM, double>> p_scr_var = new CellVariable<NDIM, double>(d_object_name + "::p_scr");

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
    if (var_db->checkVariableExists(d_object_name + "::p_scr"))
    {
        p_scr_var = var_db->getVariable(d_object_name + "::p_scr");
        d_p_scr_idx = var_db->mapVariableAndContextToIndex(p_scr_var, d_ctx);
        var_db->removePatchDataIndex(d_p_scr_idx);
    }
    d_un_scr_idx = var_db->registerVariableAndContext(un_scr_var, d_ctx, IntVector<NDIM>(1));
    d_us_scr_idx = var_db->registerVariableAndContext(us_scr_var, d_ctx, IntVector<NDIM>(1));
    d_p_scr_idx = var_db->registerVariableAndContext(p_scr_var, d_ctx, IntVector<NDIM>(1));

    if (var_db->checkVariableExists(d_mask_var->getName()))
    {
        d_mask_var = var_db->getVariable(d_mask_var->getName());
        d_mask_idx = var_db->mapVariableAndContextToIndex(d_mask_var, d_ctx);
        var_db->removePatchDataIndex(d_mask_idx);
    }
    d_mask_idx = var_db->registerVariableAndContext(d_mask_var, d_ctx, IntVector<NDIM>(0));

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStaggeredStokesBoxRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStaggeredStokesBoxRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "multiphase::MultiphaseStaggeredStokesBoxRelaxationFACOperator::computeResidual()"););
    return;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
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
    if (P_bc_coef)
        d_P_bc_coef = P_bc_coef;
    else
        d_P_bc_coef = d_default_P_bc_coef.get();
}

MultiphaseStaggeredStokesBoxRelaxationFACOperator::~MultiphaseStaggeredStokesBoxRelaxationFACOperator()
{
    // Dallocate operator state first
    deallocateOperatorState();
    return;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::restrictResidual(const SAMRAIVectorReal<NDIM, double>& src,
                                                                    SAMRAIVectorReal<NDIM, double>& dst,
                                                                    const int dst_ln)
{
    // Pull out patch indices.
    const int dst_un_idx = dst.getComponentDescriptorIndex(0);
    const int dst_us_idx = dst.getComponentDescriptorIndex(1);
    const int dst_p_idx = dst.getComponentDescriptorIndex(2);
    std::array<int, 3> dst_idxs = { dst_un_idx, dst_us_idx, dst_p_idx };

    const int src_un_idx = src.getComponentDescriptorIndex(0);
    const int src_us_idx = src.getComponentDescriptorIndex(1);
    const int src_p_idx = src.getComponentDescriptorIndex(2);
    std::array<int, 3> src_idxs = { src_un_idx, src_us_idx, src_p_idx };

    // Now perform restriction
    performRestriction(dst_idxs, src_idxs, dst_ln);
    return;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::prolongError(const SAMRAIVectorReal<NDIM, double>& src,
                                                                SAMRAIVectorReal<NDIM, double>& dst,
                                                                const int dst_ln)
{
    // Pull out patch indices.
    const int dst_un_idx = dst.getComponentDescriptorIndex(0);
    const int dst_us_idx = dst.getComponentDescriptorIndex(1);
    const int dst_p_idx = dst.getComponentDescriptorIndex(2);
    std::array<int, 3> dst_idxs = { dst_un_idx, dst_us_idx, dst_p_idx };

    const int src_un_idx = src.getComponentDescriptorIndex(0);
    const int src_us_idx = src.getComponentDescriptorIndex(1);
    const int src_p_idx = src.getComponentDescriptorIndex(2);
    std::array<int, 3> src_idxs = { src_un_idx, src_us_idx, src_p_idx };

    // Now perform prolongation
    performProlongation(dst_idxs, src_idxs, dst_ln);
    return;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::prolongErrorAndCorrect(const SAMRAIVectorReal<NDIM, double>& src,
                                                                          SAMRAIVectorReal<NDIM, double>& dst,
                                                                          const int dst_ln)
{
    // Pull out patch indices.
    const int dst_un_idx = dst.getComponentDescriptorIndex(0);
    const int dst_us_idx = dst.getComponentDescriptorIndex(1);
    const int dst_p_idx = dst.getComponentDescriptorIndex(2);
    std::array<int, 3> dst_idxs = { d_un_scr_idx, d_us_scr_idx, d_p_scr_idx };

    const int src_un_idx = src.getComponentDescriptorIndex(0);
    const int src_us_idx = src.getComponentDescriptorIndex(1);
    const int src_p_idx = src.getComponentDescriptorIndex(2);
    std::array<int, 3> src_idxs = { src_un_idx, src_us_idx, src_p_idx };

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
    if (dst_p_idx != src_p_idx)
    {
        HierarchyCellDataOpsReal<NDIM, double> level_cc_data_ops(d_hierarchy, dst_ln - 1, dst_ln - 1);
        level_cc_data_ops.add(dst_p_idx, dst_p_idx, src_p_idx, false /*interior_only*/);
    }

    // Now prolong and correct data
    performProlongation(dst_idxs, src_idxs, dst_ln);
    HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
    level_sc_data_ops.add(dst_un_idx, dst_un_idx, d_un_scr_idx, false /*interior_only*/);
    level_sc_data_ops.add(dst_us_idx, dst_us_idx, d_us_scr_idx, false /*interior_only*/);
    HierarchyCellDataOpsReal<NDIM, double> level_cc_data_ops(d_hierarchy, dst_ln, dst_ln);
    level_cc_data_ops.add(dst_p_idx, dst_p_idx, d_p_scr_idx, false /*interior_only*/);
    return;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::smoothError(
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
    const int thn_cc_idx = d_thn_manager->getCellIndex();
    const int thn_sc_idx = d_thn_manager->getSideIndex();
    const int thn_nc_idx = d_thn_manager->getNodeIndex();
    const int f_un_idx = residual.getComponentDescriptorIndex(0); // RHS_Un
    const int f_us_idx = residual.getComponentDescriptorIndex(1); // RHS_Us
    const int f_P_idx = residual.getComponentDescriptorIndex(2);  // RHS_pressure

    d_hierarchy = error.getPatchHierarchy();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > 0 && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(P_idx);

            Pointer<SideData<NDIM, double>> un_scr_data = patch->getPatchData(d_un_scr_idx);
            Pointer<SideData<NDIM, double>> us_scr_data = patch->getPatchData(d_us_scr_idx);
            Pointer<CellData<NDIM, double>> p_scr_data = patch->getPatchData(d_p_scr_idx);

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                un_scr_data->getArrayData(axis).copy(un_data->getArrayData(axis),
                                                     d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                     IntVector<NDIM>(0));
                us_scr_data->getArrayData(axis).copy(us_data->getArrayData(axis),
                                                     d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                     IntVector<NDIM>(0));
            }

            p_scr_data->getArrayData().copy(
                p_data->getArrayData(), d_patch_cell_bc_box_overlap[level_num][patch_counter], IntVector<NDIM>(0));
        }
    }

    // outer for loop for number of sweeps
    for (int sweep = 0; sweep < 2 * num_sweeps; sweep++)
    {
        if (level_num > 0)
        {
            // 1 here because we do red-black ordering.
            if (sweep > 1)
            {
                // Copy coarse-fine interface ghost cell values which are cached in the scratch data.
                int patch_counter = 0;
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());

                    Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
                    Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
                    Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(P_idx);

                    Pointer<SideData<NDIM, double>> un_scr_data = patch->getPatchData(d_un_scr_idx);
                    Pointer<SideData<NDIM, double>> us_scr_data = patch->getPatchData(d_us_scr_idx);
                    Pointer<CellData<NDIM, double>> p_scr_data = patch->getPatchData(d_p_scr_idx);

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        un_data->getArrayData(axis).copy(un_scr_data->getArrayData(axis),
                                                         d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                         IntVector<NDIM>(0));
                        us_data->getArrayData(axis).copy(us_scr_data->getArrayData(axis),
                                                         d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                         IntVector<NDIM>(0));
                    }

                    p_data->getArrayData().copy(p_scr_data->getArrayData(),
                                                d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                IntVector<NDIM>(0));
                }
            }
            // Fill in ghost cells. We only want to use values on our current level to fill in ghost cells.
            performGhostFilling({ un_idx, us_idx, P_idx }, level_num);

            // Compute the normal extension of the solution at coarse fine interfaces if we are not on the coarsest
            // level
            // TODO: 0 here is coarsest level number, we should set this to a variable.
            d_cc_bdry_op->setPatchDataIndex(P_idx);
            d_sc_bdry_op->setPatchDataIndices({ un_idx, us_idx });
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                // TODO: 1 here is the ghost width to fill, this should be variable.
                const IntVector<NDIM>& gcw_to_fill = 1;
                d_sc_bdry_op->computeNormalExtension(*patch, ratio, gcw_to_fill);
                d_cc_bdry_op->computeNormalExtension(*patch, ratio, gcw_to_fill);
            }
        }
        else if (sweep > 0)
        {
            performGhostFilling({ un_idx, us_idx, P_idx }, level_num);
        }

        IntVector<NDIM> xp(1, 0), yp(0, 1);

        // loop through all patches on this level
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_cc_idx);
            Pointer<NodeData<NDIM, double>> thn_nc_data = patch->getPatchData(thn_nc_idx);
            Pointer<SideData<NDIM, double>> thn_sc_data = patch->getPatchData(thn_sc_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(P_idx);
            Pointer<SideData<NDIM, double>> f_un_data = patch->getPatchData(f_un_idx);
            Pointer<SideData<NDIM, double>> f_us_data = patch->getPatchData(f_us_idx);
            Pointer<CellData<NDIM, double>> f_p_data = patch->getPatchData(f_P_idx);
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

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            int red_or_black = sweep % 2; // red = 0 and black = 1
            if (d_bc_un_helper->patchTouchesDirichletBoundary(patch) ||
                d_bc_us_helper->patchTouchesDirichletBoundary(patch))
            {
                if (d_params.isVariableDrag())
                {
                    Pointer<SideData<NDIM, double>> xi_data = patch->getPatchData(d_params.xi_idx);
                    R_B_G_S_var_xi_mask(dx,
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
                                        p_data->getPointer(),
                                        p_data->getGhostCellWidth().min(),
                                        f_p_data->getPointer(),
                                        f_p_data->getGhostCellWidth().min(),
                                        f_un_data->getPointer(0),
                                        f_un_data->getPointer(1),
                                        f_un_data->getGhostCellWidth().min(),
                                        f_us_data->getPointer(0),
                                        f_us_data->getPointer(1),
                                        f_us_data->getGhostCellWidth().min(),
                                        thn_data->getPointer(),
                                        thn_data->getGhostCellWidth().min(),
                                        thn_nc_data->getPointer(),
                                        thn_nc_data->getGhostCellWidth().min(),
                                        thn_sc_data->getPointer(0),
                                        thn_sc_data->getPointer(1),
                                        thn_sc_data->getGhostCellWidth().min(),
                                        d_params.eta_n,
                                        d_params.eta_s,
                                        d_params.lambda_n,
                                        d_params.lambda_s,
                                        xi_data->getPointer(0),
                                        xi_data->getPointer(1),
                                        xi_data->getGhostCellWidth().min(),
                                        d_w,
                                        d_C,
                                        d_D,
                                        red_or_black,
                                        mask_data->getPointer(0),
                                        mask_data->getPointer(1),
                                        mask_data->getGhostCellWidth().max());
                }
                else
                {
                    R_B_G_S_mask(dx,
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
                                 p_data->getPointer(),
                                 p_data->getGhostCellWidth().min(),
                                 f_p_data->getPointer(),
                                 f_p_data->getGhostCellWidth().min(),
                                 f_un_data->getPointer(0),
                                 f_un_data->getPointer(1),
                                 f_un_data->getGhostCellWidth().min(),
                                 f_us_data->getPointer(0),
                                 f_us_data->getPointer(1),
                                 f_us_data->getGhostCellWidth().min(),
                                 thn_data->getPointer(),
                                 thn_data->getGhostCellWidth().min(),
                                 thn_nc_data->getPointer(),
                                 thn_nc_data->getGhostCellWidth().min(),
                                 thn_sc_data->getPointer(0),
                                 thn_sc_data->getPointer(1),
                                 thn_sc_data->getGhostCellWidth().min(),
                                 d_params.eta_n,
                                 d_params.eta_s,
                                 d_params.lambda_n,
                                 d_params.lambda_s,
                                 d_params.nu_n,
                                 d_params.nu_s,
                                 d_params.xi,
                                 d_w,
                                 d_C,
                                 d_D,
                                 red_or_black,
                                 mask_data->getPointer(0),
                                 mask_data->getPointer(1),
                                 mask_data->getGhostCellWidth().max());
                }
            }
            else
            {
                if (d_params.isVariableDrag())
                {
                    Pointer<SideData<NDIM, double>> xi_data = patch->getPatchData(d_params.xi_idx);
                    R_B_G_S_var_xi(dx,
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
                                   p_data->getPointer(),
                                   p_data->getGhostCellWidth().min(),
                                   f_p_data->getPointer(),
                                   f_p_data->getGhostCellWidth().min(),
                                   f_un_data->getPointer(0),
                                   f_un_data->getPointer(1),
                                   f_un_data->getGhostCellWidth().min(),
                                   f_us_data->getPointer(0),
                                   f_us_data->getPointer(1),
                                   f_us_data->getGhostCellWidth().min(),
                                   thn_data->getPointer(),
                                   thn_data->getGhostCellWidth().min(),
                                   thn_nc_data->getPointer(),
                                   thn_nc_data->getGhostCellWidth().min(),
                                   thn_sc_data->getPointer(0),
                                   thn_sc_data->getPointer(1),
                                   thn_sc_data->getGhostCellWidth().min(),
                                   d_params.eta_n,
                                   d_params.eta_s,
                                   d_params.lambda_n,
                                   d_params.lambda_s,
                                   xi_data->getPointer(0),
                                   xi_data->getPointer(1),
                                   xi_data->getGhostCellWidth().min(),
                                   d_w,
                                   d_C,
                                   d_D,
                                   red_or_black);
                }
                else
                {
                    R_B_G_S(dx,
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
                            p_data->getPointer(),
                            p_data->getGhostCellWidth().min(),
                            f_p_data->getPointer(),
                            f_p_data->getGhostCellWidth().min(),
                            f_un_data->getPointer(0),
                            f_un_data->getPointer(1),
                            f_un_data->getGhostCellWidth().min(),
                            f_us_data->getPointer(0),
                            f_us_data->getPointer(1),
                            f_us_data->getGhostCellWidth().min(),
                            thn_data->getPointer(),
                            thn_data->getGhostCellWidth().min(),
                            thn_nc_data->getPointer(),
                            thn_nc_data->getGhostCellWidth().min(),
                            thn_sc_data->getPointer(0),
                            thn_sc_data->getPointer(1),
                            thn_sc_data->getGhostCellWidth().min(),
                            d_params.eta_n,
                            d_params.eta_s,
                            d_params.lambda_n,
                            d_params.lambda_s,
                            d_params.nu_n,
                            d_params.nu_s,
                            d_params.xi,
                            d_w,
                            d_C,
                            d_D,
                            red_or_black);
                }
            }
        } // patchess
    }     // num_sweeps
    performGhostFilling({ un_idx, us_idx, P_idx }, level_num);
    IBTK_TIMER_STOP(t_smooth_error);
    return;
}

bool
MultiphaseStaggeredStokesBoxRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                                      const SAMRAIVectorReal<NDIM, double>& residual,
                                                                      int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);

    smoothError(error, residual, coarsest_ln, 10, false, false);

    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
} // solveCoarsestLevel

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
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
    const int thn_cc_idx = d_thn_manager->getCellIndex();
    const int thn_nc_idx = d_thn_manager->getNodeIndex();
    const int thn_sc_idx = d_thn_manager->getSideIndex();

    d_un_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    d_us_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    d_P_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(3);
    transaction_comps[0] = InterpolationTransactionComponent(un_idx,
                                                             SC_DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_un_bc_coefs, // modifiy?
                                                             d_un_fill_pattern);
    transaction_comps[1] = InterpolationTransactionComponent(us_idx,
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
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(P_idx);
            Pointer<CellData<NDIM, double>> rhs_P_data =
                patch->getPatchData(rhs_P_idx); // result of applying operator (eqn 3)
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_cc_idx);
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            Pointer<SideData<NDIM, double>> rhs_un_data =
                patch->getPatchData(rhs_un_idx); // result of applying operator (eqn 1)
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            Pointer<SideData<NDIM, double>> rhs_us_data =
                patch->getPatchData(rhs_us_idx); // result of applying operator (eqn 2)
            Pointer<SideData<NDIM, double>> res_un_data = patch->getPatchData(res_un_idx);
            Pointer<SideData<NDIM, double>> res_us_data = patch->getPatchData(res_us_idx);
            Pointer<CellData<NDIM, double>> res_P_data = patch->getPatchData(res_P_idx);
            IntVector<NDIM> xp(1, 0), yp(0, 1);

            applyCoincompressibility(patch, res_P_idx, un_idx, us_idx, thn_sc_idx, 1.0);
            if (d_params.isVariableDrag())
                accumulateMomentumForcesOnPatchVariableDrag(patch,
                                                            res_un_idx,
                                                            res_us_idx,
                                                            P_idx,
                                                            un_idx,
                                                            us_idx,
                                                            thn_cc_idx,
                                                            thn_nc_idx,
                                                            thn_sc_idx,
                                                            d_params,
                                                            d_C,
                                                            d_D,
                                                            d_D);
            else
                accumulateMomentumForcesOnPatchConstantCoefficient(patch,
                                                                   res_un_idx,
                                                                   res_us_idx,
                                                                   P_idx,
                                                                   un_idx,
                                                                   us_idx,
                                                                   thn_cc_idx,
                                                                   thn_nc_idx,
                                                                   thn_sc_idx,
                                                                   d_params,
                                                                   d_C,
                                                                   d_D,
                                                                   d_D);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                (*res_P_data)(idx) = (*rhs_P_data)(idx) - (*res_P_data)(idx);
            }

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
MultiphaseStaggeredStokesBoxRelaxationFACOperator::setToZero(SAMRAIVectorReal<NDIM, double>& vec, int level_num)
{
    const int un_idx = vec.getComponentDescriptorIndex(0); // network velocity, Un
    const int us_idx = vec.getComponentDescriptorIndex(1); // solvent velocity, Us
    const int P_idx = vec.getComponentDescriptorIndex(2);  // pressure
    d_hierarchy = vec.getPatchHierarchy();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
        Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(P_idx);
        un_data->fillAll(0.0);
        us_data->fillAll(0.0);
        p_data->fillAll(0.0);
    }
    return;
} // setToZero

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& sol,
                                                                           const SAMRAIVectorReal<NDIM, double>& rhs)
{
    if (d_is_initialized) deallocateOperatorState();
    // Setup solution and rhs vectors.
    Pointer<SideVariable<NDIM, double>> un_sol_var = sol.getComponentVariable(0);
    Pointer<SideVariable<NDIM, double>> us_sol_var = sol.getComponentVariable(1);
    Pointer<CellVariable<NDIM, double>> p_sol_var = sol.getComponentVariable(2);

    Pointer<SideVariable<NDIM, double>> un_rhs_var = rhs.getComponentVariable(0);
    Pointer<SideVariable<NDIM, double>> us_rhs_var = rhs.getComponentVariable(1);
    Pointer<CellVariable<NDIM, double>> p_rhs_var = rhs.getComponentVariable(2);

    d_hierarchy = sol.getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Rudimentary error checking.
#if !defined(NDEBUG)
    TBOX_ASSERT(un_sol_var && us_sol_var && p_sol_var);
    TBOX_ASSERT(un_rhs_var && us_rhs_var && p_rhs_var);
    TBOX_ASSERT(d_hierarchy == rhs.getPatchHierarchy());
#endif

    // Allocate scratch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_un_scr_idx, d_solution_time);
        level->allocatePatchData(d_us_scr_idx, d_solution_time);
        level->allocatePatchData(d_p_scr_idx, d_solution_time);
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
    d_P_bc_op = new CartCellRobinPhysBdryOp(d_p_scr_idx, d_P_bc_coef, false);
    std::vector<RefinePatchStrategy<NDIM>*> bc_op_ptrs(3);
    bc_op_ptrs[0] = d_un_bc_op;
    bc_op_ptrs[1] = d_us_bc_op;
    bc_op_ptrs[2] = d_P_bc_op;
    d_vel_P_bc_op = std::make_unique<RefinePatchStrategySet>(bc_op_ptrs.begin(), bc_op_ptrs.end(), false);

    // Cache prolongation operators. Creating refinement schedules can be expensive for hierarchies with many levels. We
    // create the schedules here, and SAMRAI will determine whether they need to be regenerated whenever we switch patch
    // indices.
    // TODO: Set the refine type via input database or setter.
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    d_un_prolong_op = grid_geom->lookupRefineOperator(un_sol_var, "CONSERVATIVE_LINEAR_REFINE");
    d_us_prolong_op = grid_geom->lookupRefineOperator(us_sol_var, "CONSERVATIVE_LINEAR_REFINE");
    d_p_prolong_op = grid_geom->lookupRefineOperator(p_sol_var, "CONSERVATIVE_LINEAR_REFINE");
    const int un_sol_idx = sol.getComponentDescriptorIndex(0);
    const int us_sol_idx = sol.getComponentDescriptorIndex(1);
    const int p_sol_idx = sol.getComponentDescriptorIndex(2);

    d_prolong_alg = new RefineAlgorithm<NDIM>();
    d_prolong_alg->registerRefine(d_un_scr_idx, un_sol_idx, d_un_scr_idx, d_un_prolong_op);
    d_prolong_alg->registerRefine(d_us_scr_idx, us_sol_idx, d_us_scr_idx, d_us_prolong_op);
    d_prolong_alg->registerRefine(d_p_scr_idx, p_sol_idx, d_p_scr_idx, d_p_prolong_op);

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
    d_p_restrict_op = grid_geom->lookupCoarsenOperator(p_sol_var, "CONSERVATIVE_COARSEN");

    d_restrict_alg = new CoarsenAlgorithm<NDIM>();
    d_restrict_alg->registerCoarsen(d_un_scr_idx, un_sol_idx, d_un_restrict_op);
    d_restrict_alg->registerCoarsen(d_us_scr_idx, us_sol_idx, d_us_restrict_op);
    d_restrict_alg->registerCoarsen(d_p_scr_idx, p_sol_idx, d_p_restrict_op);

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
    d_ghostfill_no_restrict_alg->registerRefine(d_p_scr_idx, p_sol_idx, d_p_scr_idx, Pointer<RefineOperator<NDIM>>());
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
    d_cc_bdry_op = new CartCellDoubleQuadraticCFInterpolation();
    d_cc_bdry_op->setConsistentInterpolationScheme(false);
    d_cc_bdry_op->setPatchHierarchy(d_hierarchy);

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
MultiphaseStaggeredStokesBoxRelaxationFACOperator::deallocateOperatorState()
{
    // Delete the operators, schedules, and algorithms
    d_un_restrict_op.setNull();
    d_us_restrict_op.setNull();
    d_p_restrict_op.setNull();
    d_restrict_scheds.resize(0);
    d_restrict_alg.setNull();

    d_un_prolong_op.setNull();
    d_us_prolong_op.setNull();
    d_p_prolong_op.setNull();
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
        level->deallocatePatchData(d_p_scr_idx);
        level->deallocatePatchData(d_mask_idx);
    }

    d_is_initialized = false;
    return;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::setUnderRelaxationParamater(const double w)
{
    d_w = w;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::setCandDCoefficients(const double C, const double D)
{
    d_C = C;
    d_D = D;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::performProlongation(const std::array<int, 3>& dst_idxs,
                                                                       const std::array<int, 3>& src_idxs,
                                                                       const int dst_ln)
{
    d_un_bc_op->setPatchDataIndex(dst_idxs[0]);
    d_un_bc_op->setHomogeneousBc(true);

    d_us_bc_op->setPatchDataIndex(dst_idxs[1]);
    d_us_bc_op->setHomogeneousBc(true);

    d_P_bc_op->setPatchDataIndex(dst_idxs[2]);
    d_P_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(dst_idxs[0], src_idxs[0], dst_idxs[0], d_un_prolong_op);
    refine_alg.registerRefine(dst_idxs[1], src_idxs[1], dst_idxs[1], d_us_prolong_op);
    refine_alg.registerRefine(dst_idxs[2], src_idxs[2], dst_idxs[2], d_p_prolong_op);

    refine_alg.resetSchedule(d_prolong_scheds[dst_ln]);
    d_prolong_scheds[dst_ln]->fillData(d_solution_time);
    d_prolong_alg->resetSchedule(d_prolong_scheds[dst_ln]);
    return;
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::performRestriction(const std::array<int, 3>& dst_idxs,
                                                                      const std::array<int, 3>& src_idxs,
                                                                      const int dst_ln)
{
    CoarsenAlgorithm<NDIM> coarsen_alg;
    coarsen_alg.registerCoarsen(dst_idxs[0], src_idxs[0], d_un_restrict_op);
    coarsen_alg.registerCoarsen(dst_idxs[1], src_idxs[1], d_us_restrict_op);
    coarsen_alg.registerCoarsen(dst_idxs[2], src_idxs[2], d_p_restrict_op);

    coarsen_alg.resetSchedule(d_restrict_scheds[dst_ln]);
    d_restrict_scheds[dst_ln]->coarsenData();
    d_restrict_alg->resetSchedule(d_restrict_scheds[dst_ln]);
}

void
MultiphaseStaggeredStokesBoxRelaxationFACOperator::performGhostFilling(const std::array<int, 3>& dst_idxs,
                                                                       const int dst_ln)
{
    d_un_bc_op->setPatchDataIndex(dst_idxs[0]);
    d_un_bc_op->setHomogeneousBc(true);

    d_us_bc_op->setPatchDataIndex(dst_idxs[1]);
    d_us_bc_op->setHomogeneousBc(true);

    d_P_bc_op->setPatchDataIndex(dst_idxs[2]);
    d_P_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(dst_idxs[0], dst_idxs[0], dst_idxs[0], Pointer<RefineOperator<NDIM>>());
    refine_alg.registerRefine(dst_idxs[1], dst_idxs[1], dst_idxs[1], Pointer<RefineOperator<NDIM>>());
    refine_alg.registerRefine(dst_idxs[2], dst_idxs[2], dst_idxs[2], Pointer<RefineOperator<NDIM>>());

    refine_alg.resetSchedule(d_ghostfill_no_restrict_scheds[dst_ln]);
    d_ghostfill_no_restrict_scheds[dst_ln]->fillData(d_solution_time);
    d_ghostfill_no_restrict_alg->resetSchedule(d_ghostfill_no_restrict_scheds[dst_ln]);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
