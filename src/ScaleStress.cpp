#include <ibamr/app_namespaces.h>

#include <HierarchyCellDataOpsReal.h>
#include <SAMRAI_config.h>

#include <array>

// Local includes
#include "ScaleStress.h"

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ScaleStress::ScaleStress(std::string object_name,
                         Pointer<CellVariable<NDIM, double>> thn_var,
                         Pointer<HierarchyIntegrator> thn_integrator)
    : CartGridFunction(std::move(object_name)), d_thn_var(thn_var), d_thn_integrator(thn_integrator)
{
    // intentionally blank
} // ScaleStress

void
ScaleStress::setDataOnPatchHierarchy(const int data_idx,
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
    const int thn_cur_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_thn_integrator->getCurrentContext());
    const int thn_scr_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_thn_integrator->getScratchContext());
    const int thn_new_idx = var_db->mapVariableAndContextToIndex(d_thn_var, d_thn_integrator->getNewContext());
    bool deallocate_after = !d_thn_integrator->isAllocatedPatchData(thn_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_after) d_thn_integrator->allocatePatchData(thn_scr_idx, data_time, coarsest_ln, finest_ln);

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_thn_integrator->isAllocatedPatchData(thn_new_idx, coarsest_ln, finest_ln))
    {
        pout << " Computing linear sum\n";
        hier_cc_data_ops.linearSum(thn_scr_idx, 0.5, thn_cur_idx, 0.5, thn_new_idx);
    }
    else
    {
        hier_cc_data_ops.copyData(thn_scr_idx, thn_cur_idx);
        pout << "Copying data\n";
    }

    pout << "Here\n";
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(ln), data_time, initial_time);

    if (deallocate_after) d_thn_integrator->deallocatePatchData(thn_scr_idx, coarsest_ln, finest_ln);
}

void
ScaleStress::setDataOnPatch(const int data_idx,
                            Pointer<Variable<NDIM>> /*var*/,
                            Pointer<Patch<NDIM>> patch,
                            const double data_time,
                            const bool initial_time,
                            Pointer<PatchLevel<NDIM>> /*level*/)
{
    if (initial_time) return;
    Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(d_thn_var, d_thn_integrator->getScratchContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
    TBOX_ASSERT(thn_data);
#endif
    for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        for (int d = 0; d < Q_data->getDepth(); ++d) (*Q_data)(idx, d) = (*Q_data)(idx, d) * (*thn_data)(idx);
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
