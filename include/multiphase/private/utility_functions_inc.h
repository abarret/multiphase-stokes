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
} // namespace multiphase
#endif
