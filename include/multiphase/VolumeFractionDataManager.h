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
/*!
 * \brief Class VolumeFractionDataManager is a simple container class that keeps a cell, node, and side centered
 * representation of the volume fraction on any number of provided patch hierarchies.
 *
 * One layer of ghost cells for the cell centered quantity is maintained by this class.
 *
 * Note that data should be deallocated before the patch hierarchy is destroyed, otherwise undefined behavior could
 * occur.
 */
class VolumeFractionDataManager
{
public:
    /*!
     * \brief Constructor that uses the provided context and CellVariable to manage representations.
     *
     * Note there are no bounds checking on regularize_thn.
     */
    VolumeFractionDataManager(std::string object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> thn_var,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                              SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr,
                              double regularize_thn = std::numeric_limits<double>::quiet_NaN());

    /*!
     * \brief Constructor that uses the provided context to manage representations.
     *
     * In this case, a CellVariable managed by the this class is created.
     *
     * Note there are no bounds checking on regularize_thn.
     */
    VolumeFractionDataManager(std::string object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx,
                              SAMRAI::solv::RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr,
                              double regularize_thn = std::numeric_limits<double>::quiet_NaN());

    /*!
     * Destructor that deallocates data managed by this class.
     */
    ~VolumeFractionDataManager();

    /*!
     * Update the volume fraction managed by this class. These functions will interpolate the cell centered quantity to
     * sides and nodes.
     *
     * Also fills one layer of ghost cells for the cell centered quantity. Note that this class is not capable of
     * detecting changes in ghost cells, so if ghost cells are tampered with outside this class, this class will not
     * update them until subsequant calls to updateVolumeFraction().
     */
    ///{
    /*!
     * Update the volume fraction by copying the provided cell centered quantity.
     */
    void updateVolumeFraction(int thn_cc_idx, SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, double time);

    /*!
     * Update the volume fraction using the provided function. This will first attempt to cast the function to a
     * ThnCartGridFunction before calling setDataOnPatchHierarchy
     */
    void updateVolumeFraction(IBTK::CartGridFunction& thn_fcn,
                              SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                              double time,
                              IBTK::TimePoint time_pt,
                              bool initial_time = false);
    ///\}

    /*!
     * Returns true if the provided hierarchy has data allocated by this class.
     */
    bool isAllocated(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy) const;

    /*!
     * Deallocates data owned by this class on the provided patch hierarchy.
     */
    void deallocateData(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy);

    /*!
     * Deallocate data owned by this class on all patch hierarchies.
     */
    void deallocateData();

    /*!
     * Return the cell variable managed by this class.
     */
    inline const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>& getCellVariable() const
    {
        return d_thn_cc_var;
    }

    /*!
     * Return the cell index corresponding to the patch data stored by this class.
     *
     * Note that ghost cells are maintained by this class.
     */
    inline int getCellIndex() const
    {
        return d_thn_cc_idx;
    }

    /*!
     * Return the side variable managed by this class.
     */
    inline const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>>& getSideVariable() const
    {
        return d_thn_sc_var;
    }

    /*!
     * Return the side index corresponding to the patch data stored by this class.
     */
    inline int getSideIndex() const
    {
        return d_thn_sc_idx;
    }

    /*!
     * Return the node variable managed by this class.
     */
    inline const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>>& getNodeVariable() const
    {
        return d_thn_nc_var;
    }

    /*!
     * Return the node index corresponding to the patch data stored by this class.
     */
    inline int getNodeIndex() const
    {
        return d_thn_nc_idx;
    }

    /*!
     * Return the VariableContext used by this class to manage data.
     */
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

    std::unordered_map<SAMRAI::hier::PatchHierarchy<NDIM>*, bool> d_hierarchy_allocated_map;

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_thn_bc_coef = nullptr;

    double d_regularize_thn = std::numeric_limits<double>::quiet_NaN();
};
} // namespace multiphase
#endif