#include "multiphase/CFMultiphaseOldroydB.h"

#include <ibamr/app_namespaces.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>

#include <HierarchyCellDataOpsReal.h>
#include <SAMRAI_config.h>

#include <array>

namespace multiphase
{
CFMultiphaseOldroydB::CFMultiphaseOldroydB(std::string object_name,
                                           Pointer<CellVariable<NDIM, double>> thn_var,
                                           Pointer<CellVariable<NDIM, double>> z_var,
                                           Pointer<AdvDiffHierarchyIntegrator> integrator,
                                           Pointer<Database> input_db)
    : CFStrategy(std::move(object_name)), d_thn_var(thn_var), d_z_var(z_var), d_integrator(integrator)
{
    d_relaxation_time = input_db->getDouble("relaxation_time");
    d_alpha = input_db->getDouble("alpha");
} // CFMultiphaseOldroydB

void
CFMultiphaseOldroydB::computeRelaxation(const int R_idx,
                                        Pointer<CellVariable<NDIM, double>> R_var,
                                        const int C_idx,
                                        Pointer<CellVariable<NDIM, double>> C_var,
                                        TensorEvolutionType evolve_type,
                                        Pointer<PatchHierarchy<NDIM>> hierarchy,
                                        const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const double l_inv = 1.0 / d_relaxation_time;

    // Find thn_scr_idx.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int z_cur_idx = var_db->mapVariableAndContextToIndex(d_z_var, d_integrator->getCurrentContext());
    const int z_scr_idx = var_db->mapVariableAndContextToIndex(d_z_var, d_integrator->getScratchContext());
    const int z_new_idx = var_db->mapVariableAndContextToIndex(d_z_var, d_integrator->getNewContext());
    bool deallocate_after = !d_integrator->isAllocatedPatchData(z_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_after) d_integrator->allocatePatchData(z_scr_idx, data_time, coarsest_ln, finest_ln);

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_integrator->isAllocatedPatchData(z_new_idx, coarsest_ln, finest_ln))
        hier_cc_data_ops.linearSum(z_scr_idx, 0.5, z_cur_idx, 0.5, z_new_idx);
    else
        hier_cc_data_ops.copyData(z_scr_idx, z_cur_idx);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double>> C_data = patch->getPatchData(C_idx);
            Pointer<CellData<NDIM, double>> z_data = patch->getPatchData(z_scr_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                const double z = d_alpha * (*z_data)(idx);
                // We could be using a different evolution type.
                MatrixNd mat = convert_to_conformation_tensor(*C_data, idx, evolve_type);
                (*R_data)(idx, 0) = l_inv * (z - mat(0, 0));
                (*R_data)(idx, 1) = l_inv * (z - mat(1, 1));
                (*R_data)(idx, 2) = l_inv * (-mat(0, 1));
            }
        }
    }

    if (deallocate_after) d_integrator->deallocatePatchData(z_scr_idx, coarsest_ln, finest_ln);
}

void
CFMultiphaseOldroydB::computeStress(const int sig_idx,
                                    Pointer<CellVariable<NDIM, double>> sig_var,
                                    Pointer<PatchHierarchy<NDIM>> hierarchy,
                                    const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Find thn_scr_idx.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int thn_cur_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_integrator->getCurrentContext());
    const int thn_scr_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_integrator->getScratchContext());
    const int thn_new_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_integrator->getNewContext());
    bool deallocate_thn_after = !d_integrator->isAllocatedPatchData(thn_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_thn_after) d_integrator->allocatePatchData(thn_scr_idx, data_time, coarsest_ln, finest_ln);

    const int z_cur_idx = var_db->mapVariableAndContextToIndex(d_z_var, d_integrator->getCurrentContext());
    const int z_scr_idx = var_db->mapVariableAndContextToIndex(d_z_var, d_integrator->getScratchContext());
    const int z_new_idx = var_db->mapVariableAndContextToIndex(d_z_var, d_integrator->getNewContext());
    bool deallocate_z_after = !d_integrator->isAllocatedPatchData(z_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_z_after) d_integrator->allocatePatchData(z_scr_idx, data_time, coarsest_ln, finest_ln);

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_integrator->isAllocatedPatchData(z_new_idx, coarsest_ln, finest_ln))
        hier_cc_data_ops.linearSum(z_scr_idx, 0.5, z_cur_idx, 0.5, z_new_idx);
    else
        hier_cc_data_ops.copyData(z_scr_idx, z_cur_idx);

    if (d_integrator->isAllocatedPatchData(thn_new_idx, coarsest_ln, finest_ln))
        hier_cc_data_ops.linearSum(thn_scr_idx, 0.5, thn_cur_idx, 0.5, thn_new_idx);
    else
        hier_cc_data_ops.copyData(thn_scr_idx, thn_cur_idx);

    // Fill thn ghost cells
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp{ ITC(thn_scr_idx,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          true,
                                          "NONE",
                                          "LINEAR",
                                          true,
                                          d_integrator->getPhysicalBcCoefs(d_thn_var)),
                                      ITC(z_scr_idx,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          true,
                                          "NONE",
                                          "LINEAR",
                                          true,
                                          d_integrator->getPhysicalBcCoefs(d_z_var)) };
    HierarchyGhostCellInterpolation ghost_cell_fill;
    ghost_cell_fill.initializeOperatorState(ghost_cell_comp, hierarchy, coarsest_ln, finest_ln);
    ghost_cell_fill.fillData(data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> C_data = patch->getPatchData(sig_idx);
            Pointer<CellData<NDIM, double>> thn_data =
                patch->getPatchData(d_thn_var, d_integrator->getScratchContext());
            Pointer<CellData<NDIM, double>> z_data = patch->getPatchData(d_z_var, d_integrator->getScratchContext());
#if !defined(NDEBUG)
            TBOX_ASSERT(C_data);
            TBOX_ASSERT(thn_data);
#endif
            // Note we only need to fill in one level of ghost cells.
            // For some reason C_data has a width of three ghost cells, so we use thn_data.
            for (CellIterator<NDIM> ci(thn_data->getGhostBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                const double thn = (*thn_data)(idx);
                // Convert conformation to stress in place.
                // Note this acts on C, not a decomposition of C.
                (*C_data)(idx, 0) = thn * ((*C_data)(idx, 0) - d_alpha * (*z_data)(idx));
                (*C_data)(idx, 1) = thn * ((*C_data)(idx, 1) - d_alpha * (*z_data)(idx));
                (*C_data)(idx, 2) = thn * ((*C_data)(idx, 2));
            }
        }
    }

    if (deallocate_thn_after) d_integrator->deallocatePatchData(thn_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_z_after) d_integrator->deallocatePatchData(z_scr_idx, coarsest_ln, finest_ln);
}

} // namespace multiphase
//////////////////////////////////////////////////////////////////////////////
