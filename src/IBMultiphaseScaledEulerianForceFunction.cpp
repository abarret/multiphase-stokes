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
#include "multiphase/IBMultiphaseHierarchyIntegrator.h"
#include "multiphase/utility_functions.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace multiphase
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBMultiphaseHierarchyIntegrator::IBMultiphaseScaledEulerianForceFunction::IBMultiphaseScaledEulerianForceFunction(
    const IBMultiphaseHierarchyIntegrator* const ib_solver,
    const int ib_idx,
    const bool scale_by_solvent)
    : CartGridFunction(ib_solver->getName() + "::IBMultiphaseScaledEulerianForceFunction"),
      d_ib_solver(ib_solver),
      d_ib_idx(ib_idx),
      d_scale_by_solvent(scale_by_solvent)
{
    // intentionally blank
    return;
} // IBMultiphaseScaledEulerianForceFunction

bool
IBMultiphaseHierarchyIntegrator::IBMultiphaseScaledEulerianForceFunction::isTimeDependent() const
{
    return true;
}

void
IBMultiphaseHierarchyIntegrator::IBMultiphaseScaledEulerianForceFunction::setDataOnPatchHierarchy(
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
IBMultiphaseHierarchyIntegrator::IBMultiphaseScaledEulerianForceFunction::setDataOnPatch(
    const int data_idx,
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
    Pointer<INSVCTwoFluidStaggeredHierarchyIntegrator> ins_integrator = d_ib_solver->d_ins_hier_integrator;
    Pointer<CellData<NDIM, double>> thn_data =
        patch->getPatchData(ins_integrator->getNetworkVolumeFractionVariable(), ins_integrator->getCurrentContext());
    Pointer<SideData<NDIM, double>> f_ib_data = patch->getPatchData(d_ib_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_ib_data);
#endif
    PatchSideDataBasicOps<NDIM, double> patch_ops;
    // Multiply force and volume fraction and add to output.
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
        {
            const SideIndex<NDIM>& idx = si();
            const double thn = 0.5 * ((*thn_data)(idx.toCell(0)) + (*thn_data)(idx.toCell(1)));
            (*f_sc_data)(idx) = (*f_sc_data)(idx) + (*f_ib_data)(idx) * (d_scale_by_solvent ? convertToThs(thn) : thn);
        }
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
