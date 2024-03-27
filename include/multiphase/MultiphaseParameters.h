#ifndef included_multiphase_MultiphaseParameters
#define included_multiphase_MultiphaseParameters

#include <ibtk/ibtk_utilities.h>

namespace multiphase
{
struct MultiphaseParameters
{
    double rho;

    double eta_n;
    double eta_s;

    double xi;
    double nu_n;
    double nu_s;

    int xi_idx = IBTK::invalid_index;

    bool isVariableDrag() const
    {
        return xi_idx != IBTK::invalid_index;
    }
};
} // namespace multiphase

#endif
