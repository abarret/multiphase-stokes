#ifndef included_multiphase_VolumeFractionDataManager
#define included_multiphase_VolumeFractionDataManager

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <NodeVariable.h>
#include <RobinBcCoefStrategy.h>
#include <SideVariable.h>

#include <set>

namespace multiphase
{
class VolumeFractionDataManager
{
public:
    VolumeFractionDataManager(std::string object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> thn_var,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                              SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr,
                              double regularize_thn = std::numeric_limits<double>::quiet_NaN());

    VolumeFractionDataManager(std::string object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                              SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr,
                              double regularize_thn = std::numeric_limits<double>::quiet_NaN());

    ~VolumeFractionDataManager();

    void updateVolumeFraction(int thn_cc_idx, SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, double time);
    void updateVolumeFraction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction>& thn_fcn,
                              SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                              double time,
                              IBTK::TimePoint time_pt,
                              bool initial_time = false);

    bool isAllocated(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy) const;

    void deallocateData(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy);
    void deallocateData();

    inline const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>& getCellVariable() const
    {
        return d_thn_cc_var;
    }
    inline int getCellIndex() const
    {
        return d_thn_cc_idx;
    }

    inline const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>>& getSideVariable() const
    {
        return d_thn_sc_var;
    }
    inline int getSideIndex() const
    {
        return d_thn_sc_idx;
    }

    inline const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>& getNodeVariable() const
    {
        return d_thn_nc_var;
    }
    inline int getNodeIndex() const
    {
        return d_thn_nc_idx;
    }

    inline const SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>& getContext() const
    {
        return d_ctx;
    }

private:
    void allocatePatchData(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, int coarsest_ln, int finest_ln, double time);
    void fillNodeAndSideData(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, double time);
    void updateVolumeFractionPrivate(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                     int coarsest_ln,
                                     int finest_ln,
                                     double time);

    std::string d_object_name;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_thn_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_thn_nc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_thn_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_cc_ndim_var;
    int d_thn_cc_idx = IBTK::invalid_index, d_thn_nc_idx = IBTK::invalid_index, d_thn_sc_idx = IBTK::invalid_index,
        d_cc_ndim_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;

    std::set<SAMRAI::hier::PatchHierarchy<NDIM>*> d_hierarchies;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_thn_bc_coef = nullptr;

    double d_regularize_thn = std::numeric_limits<double>::quiet_NaN();
};
} // namespace multiphase
#endif