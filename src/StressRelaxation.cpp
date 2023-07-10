#include "ibamr/app_namespaces.h" // IWYU pragma: keep

#include "ibtk/ibtk_utilities.h"

#include "CellData.h"
#include "CellIterator.h"
#include "Patch.h"
#include "tbox/Database.h"

// Local includes
#include "StressRelaxation.h"

StressRelaxation::StressRelaxation(const std::string& object_name,
                                   Pointer<Database> input_db,
                                   Pointer<SideVariable<NDIM, double>> u_var,
                                   Pointer<HierarchyIntegrator> u_integrator,
                                   Pointer<CellVariable<NDIM, double>> thn_var,
                                   Pointer<HierarchyIntegrator> thn_integrator)
    : CFRelaxationOperator(object_name, input_db),
      d_EE_var(new CellVariable<NDIM, double>(d_object_name + "::EE_var", 3)),
      d_u_var(u_var),
      d_u_integrator(u_integrator),
      d_thn_var(thn_var),
      d_thn_integrator(thn_integrator)
{
    d_lambda = input_db->getDouble("relaxation_time");
    d_mu = input_db->getDouble("viscosity");

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_EE_idx = var_db->registerVariableAndContext(d_EE_var, var_db->getContext(d_object_name));

    return;
} // Constructor

void
StressRelaxation::setDataOnPatchHierarchy(const int data_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                          const double data_time,
                                          const bool initial_time,
                                          int coarsest_ln,
                                          int finest_ln)
{
    coarsest_ln = coarsest_ln == IBTK::invalid_level_number ? 0 : coarsest_ln;
    finest_ln = finest_ln == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln;

    // Find thn_scr_idx.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_cur_idx = var_db->mapVariableAndContextToIndex(d_u_var, d_u_integrator->getCurrentContext());
    const int u_scr_idx = var_db->mapVariableAndContextToIndex(d_u_var, d_u_integrator->getScratchContext());
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_var, d_u_integrator->getNewContext());
    bool deallocate_after = !d_u_integrator->isAllocatedPatchData(u_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_after) d_u_integrator->allocatePatchData(u_scr_idx, data_time, coarsest_ln, finest_ln);
    d_u_integrator->allocatePatchData(d_EE_idx, data_time, coarsest_ln, finest_ln);

    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_u_integrator->isAllocatedPatchData(u_new_idx, coarsest_ln, finest_ln))
        hier_sc_data_ops.linearSum(u_scr_idx, 0.5, u_cur_idx, 0.5, u_new_idx);
    else
        hier_sc_data_ops.copyData(u_scr_idx, u_cur_idx);

    const int thn_cur_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_thn_integrator->getCurrentContext());
    const int thn_scr_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_thn_integrator->getScratchContext());
    const int thn_new_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_thn_integrator->getNewContext());
    bool deallocate_thn_after = !d_thn_integrator->isAllocatedPatchData(thn_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_thn_after) d_thn_integrator->allocatePatchData(thn_scr_idx, data_time, coarsest_ln, finest_ln);

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_thn_integrator->isAllocatedPatchData(thn_new_idx, coarsest_ln, finest_ln))
        hier_cc_data_ops.linearSum(thn_scr_idx, 0.5, thn_cur_idx, 0.5, thn_new_idx);
    else
        hier_cc_data_ops.copyData(thn_scr_idx, thn_cur_idx);

    // Fill in ghost cells for u
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comps = { ITC(u_scr_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE", "LINEAR", true) };
    Pointer<HierarchyGhostCellInterpolation> hier_ghost_fill = new HierarchyGhostCellInterpolation();
    hier_ghost_fill->initializeOperatorState(ghost_cell_comps, hierarchy, coarsest_ln, finest_ln);

    HierarchyMathOps hier_math_ops("hier_math_ops", hierarchy);
    hier_math_ops.setPatchHierarchy(hierarchy);
    hier_math_ops.resetLevels(coarsest_ln, finest_ln);
    hier_math_ops.strain_rate(d_EE_idx, d_EE_var, u_scr_idx, d_u_var, hier_ghost_fill, data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(ln), data_time, initial_time);

    if (deallocate_after) d_u_integrator->deallocatePatchData(u_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_thn_after) d_thn_integrator->deallocatePatchData(thn_scr_idx, coarsest_ln, finest_ln);
    d_u_integrator->deallocatePatchData(d_EE_idx, coarsest_ln, finest_ln);
}

void
StressRelaxation::setDataOnPatch(const int data_idx,
                                 Pointer<Variable<NDIM>> /*var*/,
                                 Pointer<Patch<NDIM>> patch,
                                 const double /*data_time*/,
                                 const bool initial_time,
                                 Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CellData<NDIM, double>> ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> in_data = patch->getPatchData(d_W_cc_idx);
    Pointer<CellData<NDIM, double>> EE_data = patch->getPatchData(d_EE_idx);
    Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(d_thn_var, d_thn_integrator->getScratchContext());
    ret_data->fillAll(0.0);
    if (initial_time) return;
    const double l_inv = 1.0 / d_lambda;
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        const CellIndex<NDIM>& idx = i();
        MatrixNd mat;
#if (NDIM == 2)
        mat(0, 0) = (*in_data)(idx, 0);
        mat(1, 1) = (*in_data)(idx, 1);
        mat(0, 1) = mat(1, 0) = (*in_data)(idx, 2);

        (*ret_data)(idx, 0) = l_inv * (-mat(0, 0)) + d_mu * l_inv * /*(*thn_data)(idx) **/ (*EE_data)(idx, 0);
        (*ret_data)(idx, 1) = l_inv * (-mat(1, 1)) + d_mu * l_inv * /*(*thn_data)(idx) **/ (*EE_data)(idx, 1);
        (*ret_data)(idx, 2) = l_inv * (-mat(1, 0)) + d_mu * l_inv * /*(*thn_data)(idx) **/ (*EE_data)(idx, 2);
#endif
#if (NDIM == 3)
        mat(0, 0) = (*in_data)(idx, 0);
        mat(1, 1) = (*in_data)(idx, 1);
        mat(2, 2) = (*in_data)(idx, 2);
        mat(1, 2) = mat(2, 1) = (*in_data)(idx, 3);
        mat(0, 2) = mat(2, 0) = (*in_data)(idx, 4);
        mat(0, 1) = mat(1, 0) = (*in_data)(idx, 5);

        (*ret_data)(idx, 0) = l_inv * (-mat(0, 0));
        (*ret_data)(idx, 1) = l_inv * (-mat(1, 1));
        (*ret_data)(idx, 2) = l_inv * (-mat(2, 2));
        (*ret_data)(idx, 3) = l_inv * (-mat(1, 2));
        (*ret_data)(idx, 4) = l_inv * (-mat(0, 2));
        (*ret_data)(idx, 5) = l_inv * (-mat(0, 1));
#endif
    }
} // setDataOnPatch
