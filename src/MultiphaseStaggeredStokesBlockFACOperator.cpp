#include "multiphase/MultiphaseStaggeredStokesBlockFACOperator.h"
#include "multiphase/MultiphaseStaggeredStokesBoxRelaxationFACOperator.h"
#include "multiphase/MultiphaseStaggeredVelocityAsymmetricFACOperator.h"
#include "multiphase/MultiphaseStaggeredVelocityBlockFACOperator.h"
#include "multiphase/fd_operators.h"
#include "multiphase/utility_functions.h"

#include "ibtk/CCPoissonBoxRelaxationFACOperator.h"
#include "ibtk/CartSideDoubleQuadraticCFInterpolation.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/app_namespaces.h"

#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "tbox/MemoryDatabase.h"

namespace multiphase
{
using namespace SAMRAI;

namespace
{
constexpr int UN_COMP = 0;
constexpr int US_COMP = 1;
constexpr int P_COMP = 2;
} // namespace

MultiphaseStaggeredStokesBlockFACOperator::MultiphaseStaggeredStokesBlockFACOperator(
    const std::string& object_name,
    const std::string& default_options_prefix,
    const MultiphaseParameters& params,
    const std::unique_ptr<VolumeFractionDataManager>& thn_manager)
    : FACPreconditionerStrategy(object_name),
      d_pressure_coefs(object_name + "::pressure_coefs"),
      d_params(params),
      d_thn_manager(thn_manager)
{
    // Reuse the existing FAC transfer operators and composite residual code.
    d_delegate = new MultiphaseStaggeredStokesBoxRelaxationFACOperator(
        object_name + "::Delegate", default_options_prefix, params, thn_manager);

    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    Pointer<hier::VariableContext> ctx = var_db->getContext(object_name + "::ctx");

    d_un_corr_var = new pdat::SideVariable<NDIM, double>(object_name + "::un_corr");
    d_us_corr_var = new pdat::SideVariable<NDIM, double>(object_name + "::us_corr");
    d_p_corr_var = new pdat::CellVariable<NDIM, double>(object_name + "::p_corr");
    d_p_rhs_var = new pdat::CellVariable<NDIM, double>(object_name + "::p_rhs");
    d_un_rhs_var = new pdat::SideVariable<NDIM, double>(object_name + "::un_rhs");
    d_us_rhs_var = new pdat::SideVariable<NDIM, double>(object_name + "::us_rhs");
    d_un_scr_var = new pdat::SideVariable<NDIM, double>(object_name + "::un_scr");
    d_us_scr_var = new pdat::SideVariable<NDIM, double>(object_name + "::us_scr");
    d_fn_scr_var = new pdat::SideVariable<NDIM, double>(object_name + "::fn_scr");
    d_fs_scr_var = new pdat::SideVariable<NDIM, double>(object_name + "::fs_scr");
    d_thn_ths_sq_var = new pdat::SideVariable<NDIM, double>(object_name + "::thn_ths_sq");

    d_un_corr_idx = var_db->registerVariableAndContext(d_un_corr_var, ctx, hier::IntVector<NDIM>(1));
    d_us_corr_idx = var_db->registerVariableAndContext(d_us_corr_var, ctx, hier::IntVector<NDIM>(1));
    d_p_corr_idx = var_db->registerVariableAndContext(d_p_corr_var, ctx, hier::IntVector<NDIM>(1));
    d_p_rhs_idx = var_db->registerVariableAndContext(d_p_rhs_var, ctx, hier::IntVector<NDIM>(1));
    d_un_rhs_idx = var_db->registerVariableAndContext(d_un_rhs_var, ctx, hier::IntVector<NDIM>(1));
    d_us_rhs_idx = var_db->registerVariableAndContext(d_us_rhs_var, ctx, hier::IntVector<NDIM>(1));
    d_un_scr_idx = var_db->registerVariableAndContext(d_un_scr_var, ctx, hier::IntVector<NDIM>(1));
    d_us_scr_idx = var_db->registerVariableAndContext(d_us_scr_var, ctx, hier::IntVector<NDIM>(1));
    d_fn_scr_idx = var_db->registerVariableAndContext(d_fn_scr_var, ctx, hier::IntVector<NDIM>(1));
    d_fs_scr_idx = var_db->registerVariableAndContext(d_fs_scr_var, ctx, hier::IntVector<NDIM>(1));
    d_thn_ths_sq_idx = var_db->registerVariableAndContext(d_thn_ths_sq_var, ctx);

    d_pressure_solver_db = new tbox::MemoryDatabase(object_name + "::PressureLevelSolverDB");
    // Match the existing FAC-style smoothers: use a local box-relaxation
    // Poisson smoother instead of a Krylov level solve so refined levels see
    // proper coarse-fine interface treatment.
    d_pressure_solver_db->putString("smoother_type", "PATCH_GAUSS_SEIDEL");
    d_pressure_solver_db->putString("coarse_solver_type", "PATCH_GAUSS_SEIDEL");
    d_pressure_solver_db->putInteger("coarse_solver_max_iterations", 2);
}

MultiphaseStaggeredStokesBlockFACOperator::~MultiphaseStaggeredStokesBlockFACOperator()
{
    deallocateOperatorState();
}

void
MultiphaseStaggeredStokesBlockFACOperator::setPhysicalBcCoefs(
    const std::vector<solv::RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
    const std::vector<solv::RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
    solv::RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
    d_un_bc_coefs = un_bc_coefs;
    d_us_bc_coefs = us_bc_coefs;
    d_P_bc_coef = P_bc_coef;
    d_delegate->setPhysicalBcCoefs(un_bc_coefs, us_bc_coefs, P_bc_coef);
    if (d_velocity_block_smoother) d_velocity_block_smoother->setPhysicalBcCoefs(un_bc_coefs, us_bc_coefs);
    if (d_velocity_asymmetric_smoother) d_velocity_asymmetric_smoother->setPhysicalBcCoefs(un_bc_coefs, us_bc_coefs);
}

void
MultiphaseStaggeredStokesBlockFACOperator::setToZero(solv::SAMRAIVectorReal<NDIM, double>& error, const int level_num)
{
    d_delegate->setToZero(error, level_num);
}

void
MultiphaseStaggeredStokesBlockFACOperator::restrictResidual(const solv::SAMRAIVectorReal<NDIM, double>& source,
                                                            solv::SAMRAIVectorReal<NDIM, double>& dest,
                                                            const int dest_level_num)
{
    d_delegate->restrictResidual(source, dest, dest_level_num);
}

void
MultiphaseStaggeredStokesBlockFACOperator::prolongError(const solv::SAMRAIVectorReal<NDIM, double>& source,
                                                        solv::SAMRAIVectorReal<NDIM, double>& dest,
                                                        const int dest_level_num)
{
    d_delegate->prolongError(source, dest, dest_level_num);
}

void
MultiphaseStaggeredStokesBlockFACOperator::prolongErrorAndCorrect(const solv::SAMRAIVectorReal<NDIM, double>& source,
                                                                  solv::SAMRAIVectorReal<NDIM, double>& dest,
                                                                  const int dest_level_num)
{
    d_delegate->prolongErrorAndCorrect(source, dest, dest_level_num);
}

void
MultiphaseStaggeredStokesBlockFACOperator::smoothError(solv::SAMRAIVectorReal<NDIM, double>& error,
                                                       const solv::SAMRAIVectorReal<NDIM, double>& residual,
                                                       const int level_num,
                                                       const int num_sweeps,
                                                       bool /*performing_pre_sweeps*/,
                                                       bool /*performing_post_sweeps*/)
{
    if (num_sweeps <= 0) return;

    // The pressure smoother depends on dense-hierarchy `thn` data, which is only
    // guaranteed to exist after the outer FullFACPreconditioner copies it over.
    Pointer<math::HierarchySideDataOpsReal<NDIM, double>> sc_ops =
        Pointer<math::HierarchySideDataOpsReal<NDIM, double>>(
            new math::HierarchySideDataOpsReal<NDIM, double>(d_hierarchy, level_num, level_num));
    Pointer<math::HierarchyCellDataOpsReal<NDIM, double>> cc_ops =
        Pointer<math::HierarchyCellDataOpsReal<NDIM, double>>(
            new math::HierarchyCellDataOpsReal<NDIM, double>(d_hierarchy, level_num, level_num));

    solv::SAMRAIVectorReal<NDIM, double> corr_vec(d_object_name + "::corr", d_hierarchy, level_num, level_num);
    corr_vec.addComponent(d_un_corr_var, d_un_corr_idx, error.getControlVolumeIndex(UN_COMP), sc_ops);
    corr_vec.addComponent(d_us_corr_var, d_us_corr_idx, error.getControlVolumeIndex(US_COMP), sc_ops);
    corr_vec.addComponent(d_p_corr_var, d_p_corr_idx, error.getControlVolumeIndex(P_COMP), cc_ops);

    solv::SAMRAIVectorReal<NDIM, double> rhs_vec(d_object_name + "::rhs", d_hierarchy, level_num, level_num);
    rhs_vec.addComponent(d_un_rhs_var, d_un_rhs_idx, residual.getControlVolumeIndex(UN_COMP), sc_ops);
    rhs_vec.addComponent(d_us_rhs_var, d_us_rhs_idx, residual.getControlVolumeIndex(US_COMP), sc_ops);
    rhs_vec.addComponent(d_p_rhs_var, d_p_rhs_idx, residual.getControlVolumeIndex(P_COMP), cc_ops);

    solv::SAMRAIVectorReal<NDIM, double> u_corr_vec(d_object_name + "::u_corr", d_hierarchy, level_num, level_num);
    u_corr_vec.addComponent(d_un_corr_var, d_un_corr_idx, error.getControlVolumeIndex(UN_COMP), sc_ops);
    u_corr_vec.addComponent(d_us_corr_var, d_us_corr_idx, error.getControlVolumeIndex(US_COMP), sc_ops);

    solv::SAMRAIVectorReal<NDIM, double> u_rhs_vec(d_object_name + "::u_rhs", d_hierarchy, level_num, level_num);
    u_rhs_vec.addComponent(d_fn_scr_var, d_fn_scr_idx, residual.getControlVolumeIndex(UN_COMP), sc_ops);
    u_rhs_vec.addComponent(d_fs_scr_var, d_fs_scr_idx, residual.getControlVolumeIndex(US_COMP), sc_ops);

    solv::SAMRAIVectorReal<NDIM, double> p_corr_vec(d_object_name + "::p_corr", d_hierarchy, level_num, level_num);
    p_corr_vec.addComponent(d_p_corr_var, d_p_corr_idx, error.getControlVolumeIndex(P_COMP), cc_ops);

    solv::SAMRAIVectorReal<NDIM, double> p_rhs_vec(d_object_name + "::p_rhs", d_hierarchy, level_num, level_num);
    p_rhs_vec.addComponent(d_p_rhs_var, d_p_rhs_idx, residual.getControlVolumeIndex(P_COMP), cc_ops);

    if (!d_pressure_smoother)
    {
        updateThnThsSqData();
        initializePressureSmoother();
    }

    const int thn_cc_idx = d_thn_manager->getCellIndex();
    const int thn_sc_idx = d_thn_manager->getSideIndex();
    const int thn_nc_idx = d_thn_manager->getNodeIndex();

    IBTK::FACPreconditionerStrategy* velocity_smoother = nullptr;
    if (d_using_symmetric)
    {
        velocity_smoother = static_cast<IBTK::FACPreconditionerStrategy*>(d_velocity_block_smoother.getPointer());
    }
    else
    {
        velocity_smoother = static_cast<IBTK::FACPreconditionerStrategy*>(d_velocity_asymmetric_smoother.getPointer());
    }

    for (int sweep = 0; sweep < num_sweeps; ++sweep)
    {
        // Compute the single-level residual equation that this smoother
        // approximately inverts.
        d_delegate->computeResidual(rhs_vec, error, residual, level_num, level_num);

        sc_ops->setToScalar(d_un_corr_idx, 0.0, false);
        sc_ops->setToScalar(d_us_corr_idx, 0.0, false);
        cc_ops->setToScalar(d_p_corr_idx, 0.0, false);
        sc_ops->setToScalar(d_un_scr_idx, 0.0, false);
        sc_ops->setToScalar(d_us_scr_idx, 0.0, false);
        sc_ops->setToScalar(d_fn_scr_idx, 0.0, false);
        sc_ops->setToScalar(d_fs_scr_idx, 0.0, false);

        // First half of the Schur approximation: apply an in-level Poisson
        // smoother to approximate (G^T G)^{-1} b_p while respecting
        // coarse-fine interfaces.
        d_pressure_smoother->setToZero(p_corr_vec, level_num);
        if (level_num == d_coarsest_ln)
            d_pressure_smoother->solveCoarsestLevel(p_corr_vec, p_rhs_vec, level_num);
        else
            d_pressure_smoother->smoothError(p_corr_vec, p_rhs_vec, level_num, 2, false, false);

        multiphase_grad_on_hierarchy(
            *d_hierarchy, d_un_scr_idx, d_us_scr_idx, thn_sc_idx, d_p_corr_idx, -1.0, false, level_num, level_num);

        // The momentum block uses these side-centered fields with one-cell
        // stencils, so they also need coarse-fine normal extensions.
        fillVelocityScratchData(level_num);

        if (d_params.isVariableDrag())
        {
            TBOX_ERROR(d_object_name << "::smoothError()\n"
                                     << "  BLOCK_FAC does not yet implement the variable-drag velocity block.\n");
        }
        accumulateMomentumWithoutPressureConstantCoefficient(*d_hierarchy,
                                                             d_fn_scr_idx,
                                                             d_fs_scr_idx,
                                                             d_un_scr_idx,
                                                             d_us_scr_idx,
                                                             thn_cc_idx,
                                                             thn_nc_idx,
                                                             thn_sc_idx,
                                                             d_params,
                                                             d_C,
                                                             d_D,
                                                             level_num,
                                                             level_num);

        applyCoincompressibility(
            *d_hierarchy, d_p_rhs_idx, d_fn_scr_idx, d_fs_scr_idx, thn_sc_idx, 1.0, level_num, level_num);

        // Second half of the Schur approximation: apply the same pressure
        // smoother to approximate (G^T G)^{-1} G^T A G p_1.
        d_pressure_smoother->setToZero(p_corr_vec, level_num);
        if (level_num == d_coarsest_ln)
            d_pressure_smoother->solveCoarsestLevel(p_corr_vec, p_rhs_vec, level_num);
        else
            d_pressure_smoother->smoothError(p_corr_vec, p_rhs_vec, level_num, 2, false, false);

        sc_ops->copyData(d_fn_scr_idx, d_un_rhs_idx, false);
        sc_ops->copyData(d_fs_scr_idx, d_us_rhs_idx, false);
        multiphase_grad_on_hierarchy(
            *d_hierarchy, d_fn_scr_idx, d_fs_scr_idx, thn_sc_idx, d_p_corr_idx, -1.0, true, level_num, level_num);

        if (!d_using_symmetric)
        {
            sc_ops->add(d_fs_scr_idx, d_fs_scr_idx, d_fn_scr_idx, false);
        }

        // Finish the block application with the level-local velocity solve.
        velocity_smoother->solveCoarsestLevel(u_corr_vec, u_rhs_vec, level_num);

        sc_ops->add(error.getComponentDescriptorIndex(UN_COMP),
                    error.getComponentDescriptorIndex(UN_COMP),
                    d_un_corr_idx,
                    false);
        sc_ops->add(error.getComponentDescriptorIndex(US_COMP),
                    error.getComponentDescriptorIndex(US_COMP),
                    d_us_corr_idx,
                    false);
        cc_ops->add(
            error.getComponentDescriptorIndex(P_COMP), error.getComponentDescriptorIndex(P_COMP), d_p_corr_idx, false);
    }
}

bool
MultiphaseStaggeredStokesBlockFACOperator::solveCoarsestLevel(solv::SAMRAIVectorReal<NDIM, double>& error,
                                                              const solv::SAMRAIVectorReal<NDIM, double>& residual,
                                                              const int coarsest_ln)
{
    smoothError(error, residual, coarsest_ln, 2, false, false);
    return true;
}

void
MultiphaseStaggeredStokesBlockFACOperator::computeResidual(solv::SAMRAIVectorReal<NDIM, double>& residual,
                                                           const solv::SAMRAIVectorReal<NDIM, double>& solution,
                                                           const solv::SAMRAIVectorReal<NDIM, double>& rhs,
                                                           const int coarsest_level_num,
                                                           const int finest_level_num)
{
    d_delegate->computeResidual(residual, solution, rhs, coarsest_level_num, finest_level_num);
}

void
MultiphaseStaggeredStokesBlockFACOperator::initializeOperatorState(const solv::SAMRAIVectorReal<NDIM, double>& solution,
                                                                   const solv::SAMRAIVectorReal<NDIM, double>& rhs)
{
    if (d_hierarchy) deallocateOperatorState();

    // Cache the hierarchy range used by the outer FAC preconditioner. The
    // extracted block smoother itself always acts on one FAC level at a time.
    d_hierarchy = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln = solution.getFinestLevelNumber();

    d_delegate->initializeOperatorState(solution, rhs);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<hier::PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_un_corr_idx, d_solution_time);
        level->allocatePatchData(d_us_corr_idx, d_solution_time);
        level->allocatePatchData(d_p_corr_idx, d_solution_time);
        level->allocatePatchData(d_p_rhs_idx, d_solution_time);
        level->allocatePatchData(d_un_rhs_idx, d_solution_time);
        level->allocatePatchData(d_us_rhs_idx, d_solution_time);
        level->allocatePatchData(d_un_scr_idx, d_solution_time);
        level->allocatePatchData(d_us_scr_idx, d_solution_time);
        level->allocatePatchData(d_fn_scr_idx, d_solution_time);
        level->allocatePatchData(d_fs_scr_idx, d_solution_time);
        level->allocatePatchData(d_thn_ths_sq_idx, d_solution_time);
    }

    initializeVelocitySmoother(solution, rhs);
}

void
MultiphaseStaggeredStokesBlockFACOperator::deallocateOperatorState()
{
    if (d_velocity_block_smoother) d_velocity_block_smoother->deallocateOperatorState();
    if (d_velocity_asymmetric_smoother) d_velocity_asymmetric_smoother->deallocateOperatorState();
    if (d_pressure_smoother) d_pressure_smoother->deallocateOperatorState();
    d_pressure_smoother.setNull();
    if (d_delegate) d_delegate->deallocateOperatorState();

    if (d_hierarchy)
    {
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            Pointer<hier::PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(d_un_corr_idx);
            level->deallocatePatchData(d_us_corr_idx);
            level->deallocatePatchData(d_p_corr_idx);
            level->deallocatePatchData(d_p_rhs_idx);
            level->deallocatePatchData(d_un_rhs_idx);
            level->deallocatePatchData(d_us_rhs_idx);
            level->deallocatePatchData(d_un_scr_idx);
            level->deallocatePatchData(d_us_scr_idx);
            level->deallocatePatchData(d_fn_scr_idx);
            level->deallocatePatchData(d_fs_scr_idx);
            level->deallocatePatchData(d_thn_ths_sq_idx);
        }
    }

    d_hierarchy.setNull();
    d_coarsest_ln = IBTK::invalid_level_number;
    d_finest_ln = IBTK::invalid_level_number;
}

void
MultiphaseStaggeredStokesBlockFACOperator::setCandDCoefficients(const double C, const double D)
{
    d_C = C;
    d_D = D;
    d_delegate->setCandDCoefficients(C, D);
    // Keep the already-built velocity smoother in sync when coefficients are
    // updated after construction.
    if (d_velocity_block_smoother) d_velocity_block_smoother->setCandDCoefficients(C, D);
    if (d_velocity_asymmetric_smoother) d_velocity_asymmetric_smoother->setCandDCoefficients(C, D);
}

void
MultiphaseStaggeredStokesBlockFACOperator::setUsingSymmetric(const bool using_symmetric)
{
    if (d_using_symmetric == using_symmetric) return;
    d_using_symmetric = using_symmetric;
    if (d_velocity_block_smoother) d_velocity_block_smoother->deallocateOperatorState();
    if (d_velocity_asymmetric_smoother) d_velocity_asymmetric_smoother->deallocateOperatorState();
    d_velocity_block_smoother.setNull();
    d_velocity_asymmetric_smoother.setNull();
}

void
MultiphaseStaggeredStokesBlockFACOperator::initializePressureSmoother()
{
    if (!d_pressure_smoother)
    {
        d_pressure_smoother = new IBTK::CCPoissonBoxRelaxationFACOperator(
            d_object_name + "::PressureFAC", d_pressure_solver_db, "stokes_block_fac_pressure_");
    }

    d_pressure_coefs.setCConstant(0.0);
    d_pressure_coefs.setDPatchDataId(d_thn_ths_sq_idx);
    d_pressure_smoother->setPoissonSpecifications(d_pressure_coefs);
    d_pressure_smoother->setPhysicalBcCoef(d_P_bc_coef);
    d_pressure_smoother->setHomogeneousBc(true);

    Pointer<math::HierarchyCellDataOpsReal<NDIM, double>> cc_ops =
        Pointer<math::HierarchyCellDataOpsReal<NDIM, double>>(
            new math::HierarchyCellDataOpsReal<NDIM, double>(d_hierarchy, d_coarsest_ln, d_finest_ln));
    solv::SAMRAIVectorReal<NDIM, double> p_corr_vec(
        d_object_name + "::PressureCorr", d_hierarchy, d_coarsest_ln, d_finest_ln);
    p_corr_vec.addComponent(d_p_corr_var, d_p_corr_idx, -1, cc_ops);
    solv::SAMRAIVectorReal<NDIM, double> p_rhs_vec(
        d_object_name + "::PressureRhs", d_hierarchy, d_coarsest_ln, d_finest_ln);
    p_rhs_vec.addComponent(d_p_rhs_var, d_p_rhs_idx, -1, cc_ops);

    d_pressure_smoother->initializeOperatorState(p_corr_vec, p_rhs_vec);
}

void
MultiphaseStaggeredStokesBlockFACOperator::initializeVelocitySmoother(
    const solv::SAMRAIVectorReal<NDIM, double>& solution,
    const solv::SAMRAIVectorReal<NDIM, double>& rhs)
{
    if (d_un_bc_coefs.empty()) d_un_bc_coefs.resize(NDIM, nullptr);
    if (d_us_bc_coefs.empty()) d_us_bc_coefs.resize(NDIM, nullptr);

    Pointer<math::HierarchySideDataOpsReal<NDIM, double>> sc_ops =
        Pointer<math::HierarchySideDataOpsReal<NDIM, double>>(
            new math::HierarchySideDataOpsReal<NDIM, double>(d_hierarchy, d_coarsest_ln, d_finest_ln));

    solv::SAMRAIVectorReal<NDIM, double> u_corr_vec(d_object_name + "::U", d_hierarchy, d_coarsest_ln, d_finest_ln);
    u_corr_vec.addComponent(d_un_corr_var, d_un_corr_idx, solution.getControlVolumeIndex(UN_COMP), sc_ops);
    u_corr_vec.addComponent(d_us_corr_var, d_us_corr_idx, solution.getControlVolumeIndex(US_COMP), sc_ops);

    solv::SAMRAIVectorReal<NDIM, double> u_rhs_vec(d_object_name + "::bU", d_hierarchy, d_coarsest_ln, d_finest_ln);
    u_rhs_vec.addComponent(d_fn_scr_var, d_fn_scr_idx, rhs.getControlVolumeIndex(UN_COMP), sc_ops);
    u_rhs_vec.addComponent(d_fs_scr_var, d_fs_scr_idx, rhs.getControlVolumeIndex(US_COMP), sc_ops);

    if (d_using_symmetric)
    {
        d_velocity_block_smoother = new MultiphaseStaggeredVelocityBlockFACOperator(
            d_object_name + "::VelocityBlockFAC", "stokes_block_fac_velocity_", d_params, d_thn_manager);
        d_velocity_block_smoother->setPhysicalBcCoefs(d_un_bc_coefs, d_us_bc_coefs);
        // The extracted block smoother must propagate the same velocity
        // coefficients used by the outer Stokes operator.
        d_velocity_block_smoother->setCandDCoefficients(d_C, d_D);
        d_velocity_block_smoother->initializeOperatorState(u_corr_vec, u_rhs_vec);
    }
    else
    {
        d_velocity_asymmetric_smoother = new MultiphaseStaggeredVelocityAsymmetricFACOperator(
            d_object_name + "::VelocityAsymFAC", "stokes_block_fac_velocity_", d_params, d_thn_manager);
        d_velocity_asymmetric_smoother->setPhysicalBcCoefs(d_un_bc_coefs, d_us_bc_coefs);
        d_velocity_asymmetric_smoother->setCandDCoefficients(d_C, d_D);
        d_velocity_asymmetric_smoother->initializeOperatorState(u_corr_vec, u_rhs_vec);
    }
}

void
MultiphaseStaggeredStokesBlockFACOperator::updateThnThsSqData()
{
    const int thn_sc_idx = d_thn_manager->getSideIndex();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<hier::PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
            Pointer<pdat::SideData<NDIM, double>> thn_sc_data = patch->getPatchData(thn_sc_idx);
            Pointer<pdat::SideData<NDIM, double>> thn_ths_sq_data = patch->getPatchData(d_thn_ths_sq_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (pdat::SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const pdat::SideIndex<NDIM>& idx = si();
                    const double thn = (*thn_sc_data)(idx);
                    const double ths = convertToThs(thn);
                    (*thn_ths_sq_data)(idx) = thn * thn + ths * ths;
                }
            }
        }
    }
}

void
MultiphaseStaggeredStokesBlockFACOperator::fillVelocityScratchData(const int level_num)
{
    using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    Pointer<xfer::VariableFillPattern<NDIM>> un_fill_pattern = new IBTK::SideNoCornersFillPattern(1, false);
    Pointer<xfer::VariableFillPattern<NDIM>> us_fill_pattern = new IBTK::SideNoCornersFillPattern(1, false);
    IBTK::HierarchyGhostCellInterpolation ghost_fill;
    std::vector<ITC> ghost_comps{
        ITC(d_un_scr_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE", "LINEAR", false, d_un_bc_coefs, un_fill_pattern),
        ITC(d_us_scr_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE", "LINEAR", false, d_us_bc_coefs, us_fill_pattern)
    };
    ghost_fill.initializeOperatorState(ghost_comps, d_hierarchy, level_num, level_num);
    ghost_fill.setHomogeneousBc(true);
    ghost_fill.fillData(d_solution_time);

    if (level_num == d_coarsest_ln) return;

    Pointer<hier::PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);
    IBTK::CartSideDoubleQuadraticCFInterpolation sc_bdry_op;
    sc_bdry_op.setConsistentInterpolationScheme(false);
    sc_bdry_op.setPatchHierarchy(d_hierarchy);
    sc_bdry_op.setPatchDataIndices({ d_un_scr_idx, d_us_scr_idx });

    const hier::IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
    const hier::IntVector<NDIM> gcw_to_fill(1);
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        sc_bdry_op.computeNormalExtension(*patch, ratio, gcw_to_fill);
    }
}
} // namespace multiphase
