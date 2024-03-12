#ifndef included_multiphase_MultiphaseConvectiveManager
#define included_multiphase_MultiphaseConvectiveManager

#include <ibamr/ConvectiveOperator.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/ibtk_utilities.h>

#include <PatchHierarchy.h>
namespace multiphase
{
/*!
 * Class MultiphaseConvectiveManager is a class that encapsulates all the routines needed to approximate the momentum
 * convective operator for the multiphase equations.
 *
 * This class lazily allocates data but does not maintain state when the hierarchy changes. This class should be
 * destroyed (or at least, deallocated) before a regridding operation occurs.
 *
 * Because this class does not maintain state across a regridding operation, this class can not be used in isolation
 * with multistep algorithms.
 */
class MultiphaseConvectiveManager
{
public:
    /*!
     * Constructor that takes in a Database. The Database is searched for the string "limiter_type". The limiter must be
     * one of "UPWIND", "CUI", "FBICS", or "MGAMMA".
     */
    MultiphaseConvectiveManager(std::string object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& un_bc_coefs = {},
                                const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& us_bc_coefs = {},
                                SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr);

    /*!
     * Constructor that takes in the limiter. The limiter must be one of "UPWIND", "CUI", "FBICS", or "MGAMMA".
     */
    MultiphaseConvectiveManager(std::string object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                IBAMR::LimiterType limiter_type,
                                const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& un_bc_coefs = {},
                                const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& us_bc_coefs = {},
                                SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr);

    /*!
     * Destructor that deallocates patch data and removes patch indices from the variable database.
     */
    ~MultiphaseConvectiveManager();

    /*!
     * Set the boundary condition objects. Any can be null. This class does not assume ownership of any of the objects.
     */
    void setBoundaryConditions(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                               const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                               SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef);

    /*!
     * Deallocate the patch data associated with this object.
     */
    void deallocateData();

    /*!
     * Reset the manager without deallocating data, essentially leaving the operator in a blank state.
     *
     * This function is useful to reset data between time steps when the hierarchy does not change.
     */
    void resetData();

    /*!
     * Approximate the convective operator using the provided time stepping scheme. Either of dst_un_idx or dst_us_idx
     * may be invalid indices. If they are valid, the resulting approximation to the operator is copied into those
     * indices.
     *
     * Note that for ts_type = "FORWARD_EULER," the new indices are not used and may be invalid indices.
     *
     * Data may be copied later using fillWithConvectiveOperator() or getConvectiveIndices().
     */
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

    /*!
     * Approximate the convective operator using the provided time stepping scheme without copying resulting data to
     * indices.
     *
     * Note that for ts_type = "FORWARD_EULER," the new indices are not used and may be invalid indices.
     *
     * Data may be copied later using fillWithConvectiveOperator() or getConvectiveIndices().
     */
    void approximateConvectiveOperator(IBAMR::TimeSteppingType ts_type,
                                       double current_time,
                                       double new_time,
                                       int un_cur_idx,
                                       int us_cur_idx,
                                       int thn_cur_idx,
                                       int un_new_idx,
                                       int us_new_idx,
                                       int thn_new_idx);

    /*!
     * Fill the patch data indices with the most recent approximation of the convective operator. If either are invalid
     * indices, then that copy is skipped.
     */
    void fillWithConvectiveOperator(int dst_un_idx, int dst_us_idx);

    /*!
     * Return the pair of patch indices that contain the most recent approximation to the convective operator. The
     * network index is the first element of the pair.
     */
    std::pair<int, int> getConvectiveIndices() const;

private:
    /*!
     * Do any common construction. Sets up patch variables and indices.
     */
    void commonConstructor();

    /*!
     * Allocates patch data and initializes ghost filling algorithms.
     */
    void allocateData(double time);

    /*!
     * Return whether this object has allocated data.
     */
    bool getIsAllocated() const;

    /*!
     * Approximate the convective operator using forward Euler.
     */
    void approximateForwardEuler(double current_time, double new_time, int un_cur_idx, int us_cur_idx, int thn_idx);

    /*!
     * Approximate the convective operator using midpoint rule.
     */
    void approximateMidpointRule(double current_time,
                                 double new_time,
                                 int un_cur_idx,
                                 int us_cur_idx,
                                 int thn_cur_idx,
                                 int un_new_idx,
                                 int us_new_idx,
                                 int thn_new_idx);

    /*!
     * Approximate the convective operator using an explicit trapezoidal rule.
     */
    void approximateTrapezoidalRule(double current_time,
                                    double new_time,
                                    int un_cur_idx,
                                    int us_cur_idx,
                                    int thn_cur_idx,
                                    int un_new_idx,
                                    int us_new_idx,
                                    int thn_new_idx);

    /*!
     * Actually computes the action of the operator. Computes the fluid's momentum, computes a limited approximation at
     * control volume faces, then performs flux differencing.
     */
    void approximateOperator(int dst_un_idx, int dst_us_idx, double eval_time, int un_idx, int us_idx, int thn_idx);

    /*!
     * On a given patch, compute the advection velocity for a control volume. This uses linear interpolation to compute
     * velocities at faces of the control volumes for side centered data.
     *
     * Note that u_adv_data must be allocated on the correct boxes before hand.
     */
    void computeAdvectionVelocity(const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch,
                                  std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& u_adv_data,
                                  const SAMRAI::pdat::SideData<NDIM, double>& U_data);

    /*!
     * Compute the momentum for the network. Computes dst_idx = thn_idx * u_idx. Note that thn_idx should have at least
     * one layer of ghost cells, and they should be filled in prior to this call.
     */
    void findNetworkMomentum(int dst_idx, int thn_idx, int u_idx);

    /*!
     * Compute the momentum for the solvent. Computes dst_idx = (1 - thn_idx) * u_idx. Note that thn_idx should have at
     * least one layer of ghost cells, and they should be filled in prior to this call.
     */
    void findSolventMomentum(int dst_idx, int thn_idx, int u_idx);

    /*!
     * Interpolate the momentum to faces of the control volume of the side centered box using the specified limiter.
     * u_data should consist of the velocities at the faces of the control volumes. This can be computed via the
     * computeAdvectionVelocity() function.
     *
     * interp_data must be allocated on the correct boxes before hand.
     */
    void interpolateToSides(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch,
                            std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& interp_data,
                            const SAMRAI::pdat::SideData<NDIM, double>& mom_data,
                            const std::array<std::unique_ptr<SAMRAI::pdat::FaceData<NDIM, double>>, NDIM>& u_data,
                            const IBAMR::LimiterType& limiter);

    /*!
     * Compute flux differencing using the limited momentum data and the velocity field.
     */
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

    /*!
     * Is there an initial approximation that we can reuse?
     */
    bool d_is_initial_approximation_filled = false;

    /*!
     * Momentum variables
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_mom_var;
    int d_mom_un_idx = IBTK::invalid_index, d_mom_us_idx = IBTK::invalid_index;

    /*!
     * Convective operator variables and scratch velocity variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_N_var, d_N0_var, d_U_var;
    int d_N_un_idx = IBTK::invalid_index, d_N_us_idx = IBTK::invalid_index;
    int d_N0_un_idx = IBTK::invalid_index, d_N0_us_idx = IBTK::invalid_index;
    int d_un_scr_idx = IBTK::invalid_index, d_us_scr_idx = IBTK::invalid_index;

    /*!
     * Scratch volume fraction variables
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_thn_var;
    int d_thn_scr_idx = IBTK::invalid_index;

    /*!
     * Data operations
     */
    SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> d_hier_sc_data_ops;
    SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> d_hier_cc_data_ops;

    /*!
     * Limiter for momentum interpolation.
     */
    IBAMR::LimiterType d_limiter = IBAMR::LimiterType::UPWIND;

    /*!
     * Cached ghost filling routines.
     */
    IBTK::HierarchyGhostCellInterpolation d_thn_ghost_fill, d_u_ghost_fill;

    /*!
     * Boundary condition routines
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_un_bc_coefs, d_us_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_thn_bc_coef;

    /*!
     * Order of boundary interpolation
     */
    std::string d_bdry_interp_order = "LINEAR";
};
} // namespace multiphase

#endif
