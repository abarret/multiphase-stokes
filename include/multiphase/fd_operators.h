#ifndef included_multiphase_fd_operators
#define included_multiphase_fd_operators

#include "multiphase/MultiphaseParameters.h"

#include "tbox/Pointer.h"

#include <Patch.h>

namespace multiphase
{
/*!
 * Accumulate the forces into respective patch indices for the network and solvent on a given patch. Assumes ghost cells
 * have been filled for the velocities and volume fraction.
 *
 * Specifically, computes
 *
 * [ C*thn + A_n + D_u*xi/nu_n*thn*ths   -D_u*xi/nu_n*thn*ths                D_p*thn*grad ][un]
 * [ -D_u*xi/nu_s*thn*ths                C*ths + A_s + D_u*xi/nu_s*thn*ths   D_p*ths*grad ][us]
 *                                                                                     [p ]
 * in which
 * A_i = D_u*eta_i*div(thn*((grad+grad^T)-div*I))
 *
 * These functions results in a runtime error if params.isVariableDrag() returns true.
 */
///\{
/*!
 * Accumulate the momentum forces for constant coefficient problems with the network volume fraction has been
 * interpolated to cell nodes and cell sides.
 */
void accumulateMomentumForcesOnPatchConstantCoefficient(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                        int A_un_idx,
                                                        int A_us_idx,
                                                        int p_idx,
                                                        int un_idx,
                                                        int us_idx,
                                                        int thn_idx,
                                                        int thn_nc_idx,
                                                        int thn_sc_idx,
                                                        const MultiphaseParameters& params,
                                                        double C,
                                                        double D_u,
                                                        double D_p);
/*!
 * Accumulate the momentum forces for constant coefficient problems with the network volume fraction is only provided at
 * cell centers. In this case, the volume fraction is linearly interpolated to respective sides and nodes when
 * necessary.
 *
 * Note that no synchronization is provided on the volume fraction when linear interpolation is done.
 */
void accumulateMomentumOnPatchConstantCoefficient(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                  int A_un_idx,
                                                  int A_us_idx,
                                                  int p_idx,
                                                  int un_idx,
                                                  int us_idx,
                                                  int thn_idx,
                                                  const MultiphaseParameters& params,
                                                  double C,
                                                  double D_u,
                                                  double D_p);
///\}

/*!
 * Accumulate the forces into respective patch indices for the network and solvent on a given patch. Assumes ghost cells
 * have been filled for the velocities and volume fraction.
 *
 * Specifically, computes
 *
 * [ C*thn + A_n + D_u*xi   -D_u*xi                D_p*thn*grad ][un]
 * [ -D_u*xi                C*ths + A_s + D_u*xi   D_p*ths*grad ][us]
 *                                                           [p ]
 * in which
 * A_i = D_u*eta_i*div(thn*((grad+grad^T)-div*I))
 *
 * Note that the drag coefficient in this case is not explicitly scaled by the volume fractions.
 *
 * This function results in a runtime error if params.isVariableDrag() returns false.
 *
 * This function linearly interpolates the volume fraction to cell sides and cell nodes as needed.
 */
void accumulateMomentumForcesOnPatchVariableDrag(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                 int A_un_idx,
                                                 int A_us_idx,
                                                 int p_idx,
                                                 int un_idx,
                                                 int us_idx,
                                                 int thn_idx,
                                                 const MultiphaseParameters& params,
                                                 double C,
                                                 double D_u,
                                                 double D_p);

/*!
 * Computes the divergence of the volume averaged velocity field on a given patch.
 *
 * Specifically, computes
 *
 * D * Div(thn*un + ths*us)
 *
 * This function linearly interpolates the volume fraction to cell sides as needed.
 */
void applyCoincompressibility(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                              int A_idx,
                              int un_idx,
                              int us_idx,
                              int thn_idx,
                              double D);
} // namespace multiphase

#include "multiphase/private/fd_operators_inc.h"
#endif
