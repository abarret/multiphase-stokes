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

#include "FRMThnForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

double
fx(const VectorNd& x, double t)
{
    return 2.0 * M_PI * std::sin(2.0 * M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]) +
           8.0 * M_PI * M_PI * std::sin(2.0 * M_PI * x[0]) * std::cos(2.0 * M_PI * x[1]);
}

double
fy(const VectorNd& x, double t)
{
    return 2.0 * M_PI * std::cos(2.0 * M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]) -
           8.0 * M_PI * M_PI * std::cos(2.0 * M_PI * x[0]) * std::sin(2.0 * M_PI * x[1]);
}

double
ths(const double thn)
{
    return 1.0 - thn;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

FRMThnForcing::FRMThnForcing(Pointer<Variable<NDIM>> thn_var,
                             Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
                             bool using_thn)
    : d_thn_var(thn_var), d_adv_diff_hier_integrator(adv_diff_hier_integrator), d_using_thn(using_thn)
{
    // intentionally blank
    return;
} // FRMThnForcing

FRMThnForcing::~FRMThnForcing()
{
    // intentionally blank
    return;
} // ~FRMThnForcing

bool
FRMThnForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
FRMThnForcing::setDataOnPatchHierarchy(const int data_idx,
                                       Pointer<Variable<NDIM>> var,
                                       Pointer<PatchHierarchy<NDIM>> hierarchy,
                                       const double data_time,
                                       const bool initial_time,
                                       const int coarsest_ln_in,
                                       const int finest_ln_in)
{
    // Pull out thn index at the correct time.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int thn_cur_idx =
        var_db->mapVariableAndContextToIndex(d_thn_var, d_adv_diff_hier_integrator->getCurrentContext());
    const int thn_new_idx =
        var_db->mapVariableAndContextToIndex(d_thn_var, d_adv_diff_hier_integrator->getNewContext());
    const int thn_scr_idx =
        var_db->mapVariableAndContextToIndex(d_thn_var, d_adv_diff_hier_integrator->getScratchContext());
    bool scr_is_allocated = d_adv_diff_hier_integrator->isAllocatedPatchData(thn_scr_idx);
    if (!scr_is_allocated) d_adv_diff_hier_integrator->allocatePatchData(thn_scr_idx, data_time);

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy);

    if (d_adv_diff_hier_integrator->isAllocatedPatchData(thn_new_idx))
    {
        hier_cc_data_ops.linearSum(thn_scr_idx, 0.5, thn_cur_idx, 0.5, thn_new_idx);
    }
    else
    {
        hier_cc_data_ops.copyData(thn_scr_idx, thn_cur_idx);
    }
    // Fill ghost cells for thn_scratch.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp = { ITC(thn_scr_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE") };
    HierarchyGhostCellInterpolation hier_ghost_fill;
    hier_ghost_fill.initializeOperatorState(ghost_cell_comp, hierarchy);
    hier_ghost_fill.fillData(data_time);

    // Compute a staggered-grid approximation to -gamma*T on each patch level.
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    // Deallocate data if we allocated it.
    if (!scr_is_allocated) d_adv_diff_hier_integrator->deallocatePatchData(thn_scr_idx);

    return;
} // setDataOnPatchHierarchy

void
FRMThnForcing::setDataOnPatch(const int data_idx,
                              Pointer<Variable<NDIM>> /*var*/,
                              Pointer<Patch<NDIM>> patch,
                              const double data_time,
                              const bool initial_time,
                              Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    // Loop through all indices on patch and fill in force*thn.
    Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);
    if (initial_time) return;
    Pointer<CellData<NDIM, double>> thn_scr_data =
        patch->getPatchData(d_thn_var, d_adv_diff_hier_integrator->getScratchContext());
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();
    const hier::Index<NDIM>& idx_low = patch_box.lower();
    const int axis = NDIM - 1;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> si(patch_box, axis); si; si++)
        {
            const SideIndex<NDIM>& idx = si();
            VectorNd x;
            for (int d = 0; d < NDIM; ++d)
                x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (axis == d ? 0.0 : 0.5));
            double thn = 0.5 * ((*thn_scr_data)(idx.toCell(0)) + (*thn_scr_data)(idx.toCell(1)));
            double th = d_using_thn ? thn : ths(thn);
            (*F_data)(idx) = th * (axis == 0 ? fx(x, data_time) : fy(x, data_time));
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
