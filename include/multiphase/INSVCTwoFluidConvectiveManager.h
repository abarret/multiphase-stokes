#ifndef included_multiphase_INSVCTwoFluidConvectiveManager
#define included_multiphase_INSVCTwoFluidConvectiveManager

#include <ibamr/ConvectiveOperator.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/ibtk_utilities.h>

#include <PatchHierarchy.h>
namespace multiphase
{
/*!
 * Class INSVCTwoFluidConvectiveManager is a lightweight container that can be used to manage the data needed to
 * evaluate the convective operator and store approximations across a time step.
 */
class INSVCTwoFluidConvectiveManager
{
public:
    INSVCTwoFluidConvectiveManager(std::string object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    ~INSVCTwoFluidConvectiveManager();

    /*!
     * Deallocate the patch data associated with this object.
     */
    void deallocateData();

    void approximateConvectiveOperator(int dst_un_idx,
                                       int dst_us_idx,
                                       IBAMR::TimeSteppingType ts_type,
                                       double current_time,
                                       double new_time,
                                       int un_cur_idx,
                                       int us_cur_idx,
                                       int thn_cur_idx,
                                       int un_new_idx,
                                       int us_new_idx,
                                       int thn_new_idx);
    void approximateConvectiveOperator(IBAMR::TimeSteppingType ts_type,
                                       double current_time,
                                       double new_time,
                                       int un_cur_idx,
                                       int us_cur_idx,
                                       int thn_cur_idx,
                                       int un_new_idx,
                                       int us_new_idx,
                                       int thn_new_idx);

    void fillWithConvectiveOperator(int dst_un_idx, int dst_us_idx);

    std::pair<int, int> getConvectiveIndices() const;

private:
    void allocateData(double time, int thn);

    bool getIsAllocated() const;

    void approximateForwardEuler(double current_time, double new_time, int un_cur_idx, int us_cur_idx, int thn_idx);
    void approximateMidpointRule(double current_time,
                                 double new_time,
                                 int un_cur_idx,
                                 int us_cur_idx,
                                 int thn_cur_idx,
                                 int un_new_idx,
                                 int us_new_idx,
                                 int thn_new_idx)
    {
    }
    void approximateTrapezoidalRule(double current_time,
                                    double new_time,
                                    int un_cur_idx,
                                    int us_cur_idx,
                                    int thn_cur_idx,
                                    int un_new_idx,
                                    int us_new_idx,
                                    int thn_new_idx)
    {
    }

    void approximateOperator(int dst_un_idx, int dst_us_idx, double eval_time, int un_idx, int us_idx, int thn_idx);

    void computeAdvectionVelocity(const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch,
                                  std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& u_adv_data,
                                  const SAMRAI::pdat::SideData<NDIM, double>& U_data);

    void findNetworkMomentum(int dst_idx, int thn_idx, int u_idx);

    void findSolventMomentum(int dst_idx, int thn_idx, int u_idx);

    void interpolateToSides(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch,
                            std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& interp_data,
                            const SAMRAI::pdat::SideData<NDIM, double>& mom_data,
                            const std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& u_data,
                            const IBAMR::LimiterType& limiter);

    void fluxDifference(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch,
                        SAMRAI::pdat::SideData<NDIM, double>& N_data,
                        const std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& mom_data,
                        const std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& U_adv_data);

    /*!
     * Object name
     */
    std::string d_object_name;

    /*!
     * Patch Hierarchy object
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    /*!
     * Is data currently allocated?
     */
    bool d_is_allocated = false;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_mom_var;
    int d_mom_un_idx = IBTK::invalid_index, d_mom_us_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_N_var, d_N0_var;
    int d_N_un_idx = IBTK::invalid_index, d_N_us_idx = IBTK::invalid_index;
    int d_N0_un_idx = IBTK::invalid_index, d_N0_us_idx = IBTK::invalid_index;

    SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> d_hier_sc_data_ops;

    IBAMR::LimiterType d_limiter = IBAMR::LimiterType::UPWIND;

    IBTK::HierarchyGhostCellInterpolation d_thn_ghost_fill, d_mom_ghost_fill;
};
} // namespace multiphase

#endif
