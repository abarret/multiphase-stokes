#ifndef included_multiphase_volume_fraction_functions
#define included_multiphase_volume_fraction_functions
#include <ibtk/HierarchyGhostCellInterpolation.h>

namespace multiphase
{
void fill_ghost_cells(int idx,
                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                      int coarsest_ln,
                      int finest_ln,
                      double time,
                      std::string refine_op = "CONSERVATIVE_LINEAR_REFINE",

                      bool use_cf_bdry_interpolation = true,
                      std::string coarsen_op = "CONSERVATIVE_COARSEN",
                      std::string phys_bdry_extrap_type = "NONE",
                      bool consistent_type_2_bdry = false,
                      SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef = nullptr,
                      SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
                      const std::string& phys_bdry_type = "LINEAR");

void fill_ghost_cells(HierarchyGhostCellInterpolation& hier_ghost_cell,
                      int idx,
                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                      int coarsest_ln,
                      int finest_ln,
                      double time,
                      std::string refine_op = "CONSERVATIVE_LINEAR_REFINE",
                      bool use_cf_bdry_interpolation = true,
                      std::string coarsen_op = "CONSERVATIVE_COARSEN",
                      std::string phys_bdry_extrap_type = "NONE",
                      bool consistent_type_2_bdry = false,
                      SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef = nullptr,
                      SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
                      const std::string& phys_bdry_type = "LINEAR");

void normalize_volume_fraction(const int thn_cc_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                               double regularize_thn);
} // namespace multiphase

#include <multiphase/private/volume_fraction_functions_inc.h>
#endif