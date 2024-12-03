#ifndef included_multiphase_utility_functions_inc
#define included_multiphase_utility_functions_inc
#include "multiphase/utility_functions.h"
namespace multiphase
{
inline void
allocate_patch_data(const SAMRAI::hier::ComponentSelector& comp,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                    const double time,
                    int coarsest_ln,
                    int finest_ln)
{
    coarsest_ln = coarsest_ln < 0 ? 0 : coarsest_ln;
    finest_ln = finest_ln < 0 ? hierarchy->getFinestLevelNumber() : finest_ln;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(comp, time);
    }
}

inline void
allocate_patch_data(const std::set<int>& idxs,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                    const double time,
                    const int coarsest_ln,
                    const int finest_ln)
{
    SAMRAI::hier::ComponentSelector comp;
    for (const auto& idx : idxs) comp.setFlag(idx);
    allocate_patch_data(comp, hierarchy, time, coarsest_ln, finest_ln);
}

inline void
allocate_patch_data(const int idx,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                    const double time,
                    const int coarsest_ln,
                    const int finest_ln)
{
    SAMRAI::hier::ComponentSelector comp;
    comp.setFlag(idx);
    allocate_patch_data(comp, hierarchy, time, coarsest_ln, finest_ln);
}

inline void
deallocate_patch_data(const SAMRAI::hier::ComponentSelector& comp,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                      int coarsest_ln,
                      int finest_ln)
{
    coarsest_ln = coarsest_ln < 0 ? 0 : coarsest_ln;
    finest_ln = finest_ln < 0 ? hierarchy->getFinestLevelNumber() : finest_ln;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(comp);
    }
}

inline void
deallocate_patch_data(const std::set<int>& idxs,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                      const int coarsest_ln,
                      const int finest_ln)
{
    SAMRAI::hier::ComponentSelector comps;
    for (const auto& idx : idxs) comps.setFlag(idx);

    deallocate_patch_data(comps, hierarchy, coarsest_ln, finest_ln);
}

inline void
deallocate_patch_data(const int idx,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                      const int coarsest_ln,
                      const int finest_ln)
{
    SAMRAI::hier::ComponentSelector comps;
    comps.setFlag(idx);
    deallocate_patch_data(comps, hierarchy, coarsest_ln, finest_ln);
}

inline void
set_valid_level_numbers(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, int& coarsest_ln, int& finest_ln)
{
    coarsest_ln = coarsest_ln == IBTK::invalid_level_number ? 0 : coarsest_ln;
    finest_ln = finest_ln == IBTK::invalid_level_number ? hierarchy.getFinestLevelNumber() : finest_ln;
}

inline void
convert_to_ndim_cc(const int dst_idx, const int cc_idx, SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy)
{
    for (int ln = 0; ln <= hierarchy.getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> cc_data = patch->getPatchData(cc_idx);
            for (SAMRAI::pdat::CellIterator<NDIM> ci(dst_data->getGhostBox()); ci; ci++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& idx = ci();
                for (int d = 0; d < dst_data->getDepth(); ++d) (*dst_data)(idx, d) = (*cc_data)(idx);
            }
        }
    }
}
} // namespace multiphase
#endif
