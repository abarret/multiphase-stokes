#ifndef included_utility_functions
#define included_utility_functions

#include <CellData.h>
#include <FaceData.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>

namespace IBAMR
{

inline double
convertToThs(double Thn)
{
    return 1.0 - Thn; // Thn+Ths = 1
}

inline void
pre_div_interp(const int dst_idx,
               const int thn_idx,
               const int un_idx,
               const int us_idx,
               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> us_data = patch->getPatchData(us_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& idx = si();
                    double thn = 0.5 * ((*thn_data)(idx.toCell(0)) + (*thn_data)(idx.toCell(1)));
                    double ths = convertToThs(thn);
                    (*dst_data)(idx) = thn * (*un_data)(idx) + ths * (*us_data)(idx);
                }
            }
        }
    }
}

// Copy data from a side-centered variable to a face-centered variable.
inline void
copy_side_to_face(const int u_fc_idx,
                  const int u_sc_idx,
                  SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double>> u_f_data = patch->getPatchData(u_fc_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> u_s_data = patch->getPatchData(u_sc_idx);
            const SAMRAI::hier::Box<NDIM> box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAI::pdat::SideIterator<NDIM> i(box, axis); i; i++)
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& si = i();
                    SAMRAI::pdat::FaceIndex<NDIM> fi(si.toCell(0), axis, 1);
                    (*u_f_data)(fi) = (*u_s_data)(si);
                }
            }
        }
    }

    return;
} // copy_side_to_face
} // namespace IBAMR
#endif
