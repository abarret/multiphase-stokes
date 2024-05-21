#ifndef included_multiphase_MultiphaseParameters
#define included_multiphase_MultiphaseParameters

#include <ibtk/ibtk_utilities.h>

namespace multiphase
{
struct MultiphaseParameters
{
    double rho = std::numeric_limits<double>::quiet_NaN();

    double eta_n = std::numeric_limits<double>::quiet_NaN();
    double eta_s = std::numeric_limits<double>::quiet_NaN();
    double lambda_n = std::numeric_limits<double>::quiet_NaN();
    double lambda_s = std::numeric_limits<double>::quiet_NaN();

    double xi = std::numeric_limits<double>::quiet_NaN();
    double nu_n = std::numeric_limits<double>::quiet_NaN();
    double nu_s = std::numeric_limits<double>::quiet_NaN();

    int xi_idx = IBTK::invalid_index;

    bool isVariableDrag() const
    {
        return xi_idx != IBTK::invalid_index;
    }
};
} // namespace multiphase

#endif
