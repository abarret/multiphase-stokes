#include <ibamr/app_namespaces.h>

#include <multiphase/MultiphaseStaggeredHierarchyIntegrator.h>
#include <multiphase/VolumeFractionDataManager.h>
#include <multiphase/utility_functions.h>
#include <multiphase/volume_fraction_functions.h>

namespace multiphase
{
VolumeFractionDataManager::VolumeFractionDataManager(std::string object_name,
                                                     Pointer<CellVariable<NDIM, double>> thn_var,
                                                     Pointer<VariableContext> ctx,
                                                     RobinBcCoefStrategy<NDIM>* thn_bc_coef,
                                                     const double regularize_thn)
    : d_object_name(std::move(object_name)),
      d_thn_cc_var(thn_var),
      d_thn_nc_var(get_var<NodeVariable<NDIM, double>>(d_object_name + "::thn_nc")),
      d_thn_sc_var(get_var<SideVariable<NDIM, double>>(d_object_name + "::thn_sc")),
      d_cc_ndim_var(get_var<CellVariable<NDIM, double>>(d_object_name + "::NDIM", NDIM)),
      d_ctx(ctx),
      d_thn_bc_coef(thn_bc_coef),
      d_regularize_thn(regularize_thn)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    IntVector<NDIM> one_gcw = 1;
    IntVector<NDIM> zero_gcw = 0;
    d_thn_cc_idx = var_db->registerVariableAndContext(d_thn_cc_var, d_ctx, one_gcw);
    d_thn_sc_idx = var_db->registerVariableAndContext(d_thn_sc_var, d_ctx, zero_gcw);
    d_thn_nc_idx = var_db->registerVariableAndContext(d_thn_nc_var, d_ctx, zero_gcw);
    d_cc_ndim_idx = var_db->registerVariableAndContext(d_cc_ndim_var, d_ctx, one_gcw);
}

VolumeFractionDataManager::VolumeFractionDataManager(std::string object_name,
                                                     Pointer<VariableContext> ctx,
                                                     RobinBcCoefStrategy<NDIM>* thn_bc_coef,
                                                     const double regularize_thn)
    : d_object_name(std::move(object_name)),
      d_thn_cc_var(get_var<CellVariable<NDIM, double>>(d_object_name + "::thn_cc")),
      d_thn_nc_var(get_var<NodeVariable<NDIM, double>>(d_object_name + "::thn_nc")),
      d_thn_sc_var(get_var<SideVariable<NDIM, double>>(d_object_name + "::thn_sc")),
      d_cc_ndim_var(get_var<CellVariable<NDIM, double>>(d_object_name + "::NDIM", NDIM)),
      d_ctx(ctx),
      d_thn_bc_coef(thn_bc_coef),
      d_regularize_thn(regularize_thn)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    if (!d_ctx) d_ctx = var_db->getContext(d_object_name + "::CTX");
    IntVector<NDIM> one_gcw = 1;
    IntVector<NDIM> zero_gcw = 0;
    d_thn_cc_idx = var_db->registerVariableAndContext(d_thn_cc_var, d_ctx, one_gcw);
    d_thn_sc_idx = var_db->registerVariableAndContext(d_thn_sc_var, d_ctx, zero_gcw);
    d_thn_nc_idx = var_db->registerVariableAndContext(d_thn_nc_var, d_ctx, zero_gcw);
    d_cc_ndim_idx = var_db->registerVariableAndContext(d_cc_ndim_var, d_ctx, one_gcw);
}

VolumeFractionDataManager::~VolumeFractionDataManager()
{
    deallocateData();
}

void
VolumeFractionDataManager::updateVolumeFraction(const int thn_cc_idx,
                                                Pointer<PatchHierarchy<NDIM>>& hierarchy,
                                                const double time)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Allocate patch data
    allocatePatchData(hierarchy, coarsest_ln, finest_ln, time);

    // Copy data from the provided index
    if (thn_cc_idx != d_thn_cc_idx)
    {
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, false);
        hier_cc_data_ops.copyData(d_thn_cc_idx, thn_cc_idx);
    }

    updateVolumeFractionPrivate(hierarchy, coarsest_ln, finest_ln, time);
}

void
VolumeFractionDataManager::updateVolumeFraction(CartGridFunction& thn_fcn,
                                                Pointer<PatchHierarchy<NDIM>>& hierarchy,
                                                double time,
                                                TimePoint time_pt,
                                                bool initial_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Allocate patch data
    allocatePatchData(hierarchy, coarsest_ln, finest_ln, time);

    // Use the provided function to fill in cell data
    // Check if this is a special CartGridFunction
    if (auto thn_cart_fcn = dynamic_cast<ThnCartGridFunction*>(&thn_fcn))
        thn_cart_fcn->setDataOnPatchHierarchy(
            d_thn_cc_idx, d_thn_cc_var, hierarchy, time, time_pt, false, coarsest_ln, finest_ln);
    else
        thn_fcn.setDataOnPatchHierarchy(d_thn_cc_idx, d_thn_cc_var, hierarchy, time, false, coarsest_ln, finest_ln);

    updateVolumeFractionPrivate(hierarchy, coarsest_ln, finest_ln, time);
}

bool
VolumeFractionDataManager::isAllocated(Pointer<PatchHierarchy<NDIM>>& hierarchy) const
{
    const auto& it = d_hierarchies.find(hierarchy);
    if (it != d_hierarchies.end())
        return true;
    else
        return false;
}

void
VolumeFractionDataManager::deallocateData(Pointer<PatchHierarchy<NDIM>>& hierarchy)
{
    auto it = d_hierarchies.find(hierarchy);
    if (it != d_hierarchies.end()) deallocateDataPrivate(it);
}

void
VolumeFractionDataManager::deallocateData()
{
    for (auto it = d_hierarchies.begin(); it != d_hierarchies.end();) it = deallocateDataPrivate(it);
}

VolumeFractionDataManager::iterator
VolumeFractionDataManager::deallocateDataPrivate(const iterator& it)
{
    deallocate_patch_data(
        { d_thn_cc_idx, d_thn_nc_idx, d_thn_sc_idx, d_cc_ndim_idx }, **it, 0, (*it)->getFinestLevelNumber());
    return d_hierarchies.erase(it);
}

void
VolumeFractionDataManager::allocatePatchData(Pointer<PatchHierarchy<NDIM>>& hierarchy,
                                             int coarsest_ln,
                                             int finest_ln,
                                             double time)
{
    allocate_patch_data(
        { d_thn_cc_idx, d_thn_nc_idx, d_thn_sc_idx, d_cc_ndim_idx }, *hierarchy, time, coarsest_ln, finest_ln);
    d_hierarchies.insert(hierarchy);
}
void
VolumeFractionDataManager::fillNodeAndSideData(Pointer<PatchHierarchy<NDIM>>& hierarchy, double time)
{
    HierarchyMathOps hier_math_ops("Hier_math_ops", hierarchy);

    hier_math_ops.interp(d_thn_nc_idx, d_thn_nc_var, true, d_thn_cc_idx, d_thn_cc_var, nullptr, time);

    convert_to_ndim_cc(d_cc_ndim_idx, d_thn_cc_idx, *hierarchy);
    hier_math_ops.interp(d_thn_sc_idx, d_thn_sc_var, true, d_cc_ndim_idx, d_cc_ndim_var, nullptr, time);
}

void
VolumeFractionDataManager::updateVolumeFractionPrivate(Pointer<PatchHierarchy<NDIM>>& hierarchy,
                                                       const int coarsest_ln,
                                                       const int finest_ln,
                                                       const double time)
{
    // Fill ghost cells
    fill_ghost_cells(d_thn_cc_idx,
                     hierarchy,
                     coarsest_ln,
                     finest_ln,
                     time,
                     "CONSERVATIVE_LINEAR_REFINE",
                     true,
                     "NONE",
                     "LINEAR",
                     true,
                     d_thn_bc_coef);

    // Normalize volume fraction, including ghost cells
    normalize_volume_fraction(d_thn_cc_idx, *hierarchy, d_regularize_thn, coarsest_ln, finest_ln);

    // Fill node and side data
    fillNodeAndSideData(hierarchy, time);
}
} // namespace multiphase