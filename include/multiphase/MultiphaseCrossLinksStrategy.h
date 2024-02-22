#ifndef included_multiphase_MultiphaseCrossLinksStrategy
#define included_multiphase_MultiphaseCrossLinksStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/CartGridFunction.h"
#include <ibtk/ibtk_enums.h>

#include "CoarsenPatchStrategy.h"
#include "IntVector.h"
#include "RefinePatchStrategy.h"
#include "StandardTagAndInitStrategy.h"
#include "VariableContext.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace multiphase
{
/*!
 * \brief Class MultiphaseCrossLinksStrategy computes the penalized cross-link force between the two representations of
 * the immersed structure.
 */
class MultiphaseCrossLinksStrategy : public SAMRAI::tbox::DescribedClass
{
public:

    /*!
     * \brief Constructor.
     */
    MultiphaseCrossLinksStrategy() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~MultiphaseCrossLinksStrategy() = default;

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    inline virtual void computeLagrangianForce(double data_time, IBTK::TimePoint time_pt)
    {
        doComputeLagrangianForce(data_time, time_pt);
    }

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    inline virtual void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>>& f_prolongation_scheds,
                double data_time,
                bool spread_network,
                IBTK::TimePoint time_pt)
    {
        doSpreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time, spread_network, time_pt);
    }

private:
    virtual void doComputeLagrangianForce(double data_time, IBTK::TimePoint time_pt) = 0;

    virtual void
    doSpreadForce(int f_data_idx,
                  IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                  const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>>& f_prolongation_scheds,
                  double data_time,
                  bool spread_network,
                  IBTK::TimePoint time_pt) = 0;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
};
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_multiphase_MultiphaseCrossLinksStrategy
