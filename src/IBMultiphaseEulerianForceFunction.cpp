// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/CartGridFunction.h"

#include "CellData.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataBasicOps.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataBasicOps.h"
#include "SideData.h"
#include "Variable.h"
#include "tbox/Pointer.h"

#include <string>

// Local includes
#include "IBMultiphaseHierarchyIntegrator.h"
#include "utility_functions.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBMultiphaseHierarchyIntegrator::IBMultiphaseEulerianForceFunction::IBMultiphaseEulerianForceFunction(
    const IBMultiphaseHierarchyIntegrator* const ib_solver,
    bool solvent_force)
    : CartGridFunction(ib_solver->getName() + "::IBMultiphaseEulerianForceFunction"),
      d_ib_solver(ib_solver),
      d_solvent_force(solvent_force)
{
    // intentionally blank
    return;
} // IBMultiphaseEulerianForceFunction

void
IBMultiphaseHierarchyIntegrator::IBMultiphaseEulerianForceFunction::setNetworkVolumeFractionPatchIndex(
    const int thn_idx)
{
    d_thn_idx = thn_idx;
}

bool
IBMultiphaseHierarchyIntegrator::IBMultiphaseEulerianForceFunction::isTimeDependent() const
{
    return true;
}

void
IBMultiphaseHierarchyIntegrator::IBMultiphaseEulerianForceFunction::setDataOnPatchHierarchy(
    const int data_idx,
    Pointer<Variable<NDIM>> var,
    Pointer<PatchHierarchy<NDIM>> hierarchy,
    const double data_time,
    const bool initial_time,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    if (initial_time)
    {
        d_ib_solver->d_hier_velocity_data_ops->setToScalar(data_idx, 0.0);
        return;
    }
    if (d_ib_solver->d_body_force_fcn)
    {
        d_ib_solver->d_body_force_fcn->setDataOnPatchHierarchy(
            data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
    }
    else
    {
        d_ib_solver->d_hier_velocity_data_ops->setToScalar(data_idx, 0.0);
    }

    // Set network volume fraction
    Pointer<INSVCTwoFluidStaggeredHierarchyIntegrator> ins_integrator = d_ib_solver->d_ins_hier_integrator;
    Pointer<CellVariable<NDIM, double>> thn_var = ins_integrator->getNetworkVolumeFractionVariable();
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_thn_idx = var_db->mapVariableAndContextToIndex(thn_var, ins_integrator->getCurrentContext());

    const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    return;
} // setDataOnPatchHierarchy

void
IBMultiphaseHierarchyIntegrator::IBMultiphaseEulerianForceFunction::setDataOnPatch(const int data_idx,
                                                                                   Pointer<Variable<NDIM>> /*var*/,
                                                                                   Pointer<Patch<NDIM>> patch,
                                                                                   const double /*data_time*/,
                                                                                   const bool initial_time,
                                                                                   Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<SideData<NDIM, double>> f_sc_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_sc_data);
#endif
    if (initial_time)
    {
        f_sc_data->fillAll(0.0);
        return;
    }
    Pointer<SideData<NDIM, double>> f_ib_data =
        patch->getPatchData(d_solvent_force ? d_ib_solver->d_f_idx : d_ib_solver->d_fn_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_ib_data);
#endif
    Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(d_thn_idx);
    Pointer<SideData<NDIM, double>> temp_data = new SideData<NDIM, double>(patch->getBox(), /*depth*/ 1, /*ghosts*/ 0);
    PatchSideDataBasicOps<NDIM, double> patch_ops;
    // Multiply force and volume fraction and add to output.
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
        {
            const SideIndex<NDIM>& idx = si();
            double thn = 0.5 * ((*thn_data)(idx.toCell(0)) + (*thn_data)(idx.toCell(1)));
            (*f_sc_data)(idx) = (*f_sc_data)(idx) + (*f_ib_data)(idx) * (d_solvent_force ? convertToThs(thn) : thn);
        }
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
