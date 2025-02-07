#ifndef included_multiphase_utility_functions
#define included_multiphase_utility_functions

#include <CellData.h>
#include <FaceData.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>

#include <set>

namespace multiphase
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

// Multiply thn and side centered quantity and fill in the side centered dst_idx. Note that cell centered data should
// already have ghost cells filled. dst_idx can be the same as sc_idx.
inline void
multiply_sc_and_thn(const int dst_idx,
                    const int sc_idx,
                    const int thn_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                    const bool extended_box = false)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> sc_data = patch->getPatchData(sc_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
            const SAMRAI::hier::Box<NDIM>& box = extended_box ? sc_data->getGhostBox() : patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAI::pdat::SideIterator<NDIM> si(box, axis); si; si++)
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& idx = si();
                    (*dst_data)(idx) =
                        0.5 * ((*thn_data)(idx.toCell(1)) + (*thn_data)(idx.toCell(0))) * (*sc_data)(idx);
                }
            }
        }
    }
}

// Multiply ths and side centered quantity and fill in the side centered dst_idx. Note that cell centered data should
// already have ghost cells filled. dst_idx can be the same as sc_idx.
inline void
multiply_sc_and_ths(const int dst_idx,
                    const int sc_idx,
                    const int thn_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                    const bool extended_box = false)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> sc_data = patch->getPatchData(sc_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
            const SAMRAI::hier::Box<NDIM>& box = extended_box ? sc_data->getGhostBox() : patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAI::pdat::SideIterator<NDIM> si(box, axis); si; si++)
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& idx = si();
                    (*dst_data)(idx) =
                        (1.0 - 0.5 * ((*thn_data)(idx.toCell(1)) + (*thn_data)(idx.toCell(0)))) * (*sc_data)(idx);
                }
            }
        }
    }
}

// Multiply ths and side centered quantity and fill in the side centered dst_idx. Note that cell centered data should
// already have ghost cells filled. dst_idx can be the same as sc_idx.
inline void
multiply_sc_and_thn_and_ths(const int dst_idx,
                            const int sc_idx,
                            const int thn_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                            const bool extended_box = false)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> sc_data = patch->getPatchData(sc_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
            const SAMRAI::hier::Box<NDIM>& box = extended_box ? sc_data->getGhostBox() : patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAI::pdat::SideIterator<NDIM> si(box, axis); si; si++)
                {
                    const SAMRAI::pdat::SideIndex<NDIM>& idx = si();
                    const double thn = 0.5 * ((*thn_data)(idx.toCell(1)) + (*thn_data)(idx.toCell(0)));
                    (*dst_data)(idx) = thn * convertToThs(thn) * (*sc_data)(idx);
                }
            }
        }
    }
}

/*!
 * \brief Allocate patch data on the specified levels and at the specified time.
 *
 * Note data should be deallocated manually to avoid memory leaks
 */
//\{
void allocate_patch_data(const SAMRAI::hier::ComponentSelector& comp,
                         const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                         double time,
                         int coarsest_ln,
                         int finest_ln);
void allocate_patch_data(const std::set<int>& idxs,
                         const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                         double time,
                         int coarsest_ln,
                         int finest_ln);
void allocate_patch_data(const int idx,
                         const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                         double time,
                         int coarsest_ln,
                         int finest_ln);
//\}

/*!
 * \brief Deallocate patch data on the specified levels.
 */
//\{
void deallocate_patch_data(const SAMRAI::hier::ComponentSelector& comp,
                           const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                           int coarsest_ln,
                           int finest_ln);
void deallocate_patch_data(const std::set<int>& idxs,
                           const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                           int coarsest_ln,
                           int finest_ln);
void deallocate_patch_data(const int idx,
                           const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                           int coarsest_ln,
                           int finest_ln);
//\}

template <typename T>
inline T
string_to_enum(const std::string& /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return -1;
} // string_to_enum

template <typename T>
inline std::string
enum_to_string(T /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return "UNKNOWN";
} // enum_to_string

enum class TimeSteppingType
{
    ADAMS_BASHFORTH,
    BACKWARD_EULER,
    FORWARD_EULER,
    MIDPOINT_RULE,
    TRAPEZOIDAL_RULE,
    BDF2,
    UNKNOWN_TIME_STEPPING_TYPE = -1
};

template <>
inline TimeSteppingType
string_to_enum<TimeSteppingType>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "ADAMS_BASHFORTH") == 0) return TimeSteppingType::ADAMS_BASHFORTH;
    if (strcasecmp(val.c_str(), "BACKWARD_EULER") == 0) return TimeSteppingType::BACKWARD_EULER;
    if (strcasecmp(val.c_str(), "FORWARD_EULER") == 0) return TimeSteppingType::FORWARD_EULER;
    if (strcasecmp(val.c_str(), "MIDPOINT_RULE") == 0) return TimeSteppingType::MIDPOINT_RULE;
    if (strcasecmp(val.c_str(), "TRAPEZOIDAL_RULE") == 0) return TimeSteppingType::TRAPEZOIDAL_RULE;
    if (strcasecmp(val.c_str(), "CRANK_NICOLSON") == 0) return TimeSteppingType::TRAPEZOIDAL_RULE;
    if (strcasecmp(val.c_str(), "BDF2") == 0) return TimeSteppingType::BDF2;
    return TimeSteppingType::UNKNOWN_TIME_STEPPING_TYPE;
} // string_to_enum

template <>
inline std::string
enum_to_string<TimeSteppingType>(TimeSteppingType val)
{
    if (val == TimeSteppingType::ADAMS_BASHFORTH) return "ADAMS_BASHFORTH";
    if (val == TimeSteppingType::BACKWARD_EULER) return "BACKWARD_EULER";
    if (val == TimeSteppingType::FORWARD_EULER) return "FORWARD_EULER";
    if (val == TimeSteppingType::MIDPOINT_RULE) return "MIDPOINT_RULE";
    if (val == TimeSteppingType::TRAPEZOIDAL_RULE) return "TRAPEZOIDAL_RULE";
    if (val == TimeSteppingType::BDF2) return "BDF2";
    return "UNKNOWN_TIME_STEPPING_TYPE";
} // enum_to_string

} // namespace multiphase

#include "multiphase/private/utility_functions_inc.h"
#endif
