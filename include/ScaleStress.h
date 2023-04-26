#ifndef included_ScaleStress
#define included_ScaleStress

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>

#include <CartesianGridGeometry.h>

// Local includes
/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
class ScaleStress : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    ScaleStress(std::string object_name,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> thn_var,
                SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> thn_integrator);

    /*!
     * \brief Destructor.
     */
    ~ScaleStress() = default;

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const
    {
        return true;
    }

    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = IBTK::invalid_level_number,
                                 const int finest_ln = IBTK::invalid_level_number) override;

    /*!
     * Set the data on the patch interior to the exact answer.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(NULL));

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    ScaleStress();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ScaleStress(const ScaleStress& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ScaleStress& operator=(const ScaleStress& that);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_thn_var;
    SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> d_thn_integrator;
};

/////////////////////////////// INLINE ///////////////////////////////////////

// #include "ScaleStress.I"

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_ScaleStress
