/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_multiphase_IBMultiphaseCrossLinks
#define included_multiphase_IBMultiphaseCrossLinks

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBMethod.h>

#include "ibtk/CartGridFunction.h"
#include <ibtk/LData.h>

#include "CoarsenPatchStrategy.h"
#include "IntVector.h"
#include "RefinePatchStrategy.h"
#include "StandardTagAndInitStrategy.h"
#include "VariableContext.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#include <multiphase/MultiphaseCrossLinksStrategy.h>

#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace multiphase
{
/*!
 * \brief Class IBMultiphaseCrossLinks computes the penalized cross-link force between the two representations of the
 * immersed structure.
 */
class IBMultiphaseCrossLinks : public MultiphaseCrossLinksStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBMultiphaseCrossLinks(
        SAMRAI::tbox::Pointer<IBAMR::IBMethod> ibn_ops,
        SAMRAI::tbox::Pointer<IBAMR::IBMethod> ibs_ops,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
        double kappa,
        double etan = 0.0,
        std::function<IBTK::VectorNd(double)> vel_fcn = [](double) -> IBTK::VectorNd
        { return IBTK::VectorNd::Zero(); });

    /*!
     * \brief Destructor.
     */
    virtual ~IBMultiphaseCrossLinks() = default;

private:
    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    virtual void doComputeLagrangianForce(double data_time, IBTK::TimePoint time_pt) override;

    void doComputeLagrangianForce(SAMRAI::tbox::Pointer<IBTK::LData>& Fn_data,
                                  SAMRAI::tbox::Pointer<IBTK::LData>& Fs_data,
                                  SAMRAI::tbox::Pointer<IBTK::LData>& Xn_data,
                                  SAMRAI::tbox::Pointer<IBTK::LData>& Un_data,
                                  SAMRAI::tbox::Pointer<IBTK::LData>& Xs_data,
                                  SAMRAI::tbox::Pointer<IBTK::LData>& Us_data,
                                  double time,
                                  int ln);

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    virtual void
    doSpreadForce(int f_data_idx,
                  IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                  const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>>& f_prolongation_scheds,
                  double data_time,
                  bool spread_network,
                  IBTK::TimePoint time_pt) override;

    void getNetworkPositionData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData>>& X_data, IBTK::TimePoint time_pt);
    void getNetworkVelocityData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData>>& U_data, IBTK::TimePoint time_pt);

    void getSolventPositionData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData>>& X_data, IBTK::TimePoint time_pt);
    void getSolventVelocityData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData>>& U_data, IBTK::TimePoint time_pt);

    std::vector<SAMRAI::tbox::Pointer<IBTK::LData>> d_Fn_data, d_Fs_data;

    SAMRAI::tbox::Pointer<IBAMR::IBMethod> d_ibn_ops, d_ibs_ops;

    double d_kappa = std::numeric_limits<double>::quiet_NaN();
    double d_eta = std::numeric_limits<double>::quiet_NaN();
    std::function<IBTK::VectorNd(double)> d_vel_fcn;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
};
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_multiphase_IBMultiphaseCrossLinks
