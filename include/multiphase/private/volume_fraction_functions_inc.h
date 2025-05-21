#ifndef included_multiphase_volume_fraction_functions_inc
#define included_multiphase_volume_fraction_functions_inc
#include <multiphase/volume_fraction_functions.h>

namespace multiphase
{
inline void
fill_ghost_cells(const int idx,
                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                 const int coarsest_ln,
                 const int finest_ln,
                 double time,
                 std::string refine_op,
                 bool use_cf_bdry_interpolation,
                 std::string coarsen_op,
                 std::string phys_bdry_extrap_type,
                 bool consistent_type_2_bdry,
                 SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                 SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern,
                 const std::string& phys_bdry_type)
{
    HierarchyGhostCellInterpolation hier_ghost_fill;
    fill_ghost_cells(hier_ghost_fill,
                     idx,
                     hierarchy,
                     coarsest_ln,
                     finest_ln,
                     time,
                     refine_op,
                     use_cf_bdry_interpolation,
                     coarsen_op,
                     phys_bdry_extrap_type,
                     consistent_type_2_bdry,
                     bc_coef,
                     fill_pattern,
                     phys_bdry_type);
}

inline void
fill_ghost_cells(IBTK::HierarchyGhostCellInterpolation& hier_ghost_fill,
                 const int idx,
                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                 const int coarsest_ln,
                 const int finest_ln,
                 double time,
                 std::string refine_op,
                 bool use_cf_bdry_interpolation,
                 std::string coarsen_op,
                 std::string phys_bdry_extrap_type,
                 bool consistent_type_2_bdry,
                 SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                 SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern,
                 const std::string& phys_bdry_type)
{
    using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp{ ITC(idx,
                                          refine_op,
                                          use_cf_bdry_interpolation,
                                          coarsen_op,
                                          phys_bdry_extrap_type,
                                          consistent_type_2_bdry,
                                          bc_coef,
                                          fill_pattern,
                                          phys_bdry_type) };
    hier_ghost_fill.initializeOperatorState(ghost_cell_comp, hierarchy, coarsest_ln, finest_ln);
    hier_ghost_fill.fillData(time);
}

inline void
normalize_volume_fraction(const int thn_cc_idx,
                          SAMRAI::hier::PatchHierarchy<NDIM> hierarchy,
                          double regularize_thn,
                          int coarsest_ln,
                          int finest_ln)
{
    set_valid_level_numbers(hierarchy, coarsest_ln, finest_ln);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_cc_idx);
            for (SAMRAI::pdat::CellIterator<NDIM> ci(thn_data->getGhostBox()); ci; ci++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& idx = ci();
                (*thn_data)(idx) = std::max(std::min((*thn_data)(idx), 1.0 - regularize_thn), regularize_thn);
            }
        }
    }
}
} // namespace multiphase
#endif