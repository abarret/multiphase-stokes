#ifndef included_CFMultiphaseOldroydB
#define included_CFMultiphaseOldroydB

#include <ibamr/CFStrategy.h>

#include <ibtk/HierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>

namespace multiphase
{
/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
class CFMultiphaseOldroydB : public IBAMR::CFStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    CFMultiphaseOldroydB(std::string object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> thn_var,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> z_var,
                         SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~CFMultiphaseOldroydB() = default;

    void computeRelaxation(int R_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> R_var,
                           int C_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> C_var,
                           IBAMR::TensorEvolutionType evolve_type,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                           double data_time) override;

    void computeStress(int sig_idx,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> sig_var,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                       double data_time) override;

private:
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_thn_var, d_z_var;
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_integrator;
    double d_relaxation_time = std::numeric_limits<double>::quiet_NaN();
    double d_alpha = std::numeric_limits<double>::quiet_NaN();
};
} // namespace multiphase
#endif // #ifndef included_CFMultiphaseOldroydB
