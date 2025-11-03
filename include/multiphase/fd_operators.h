#ifndef included_multiphase_fd_operators
#define included_multiphase_fd_operators

#include "multiphase/MultiphaseParameters.h"
#include "multiphase/utility_functions.h"

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
 * These functions result in a runtime error if params.isVariableDrag() returns true.
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
void accumulateMomentumForcesConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
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
                                                 double D_p,
                                                 int coarsest_ln = IBTK::invalid_level_number,
                                                 int finest_ln = IBTK::invalid_level_number);
void accumulateMomentumForcesConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                                 int A_un_idx,
                                                 int A_us_idx,
                                                 int p_idx,
                                                 int un_idx,
                                                 int us_idx,
                                                 const VolumeFractionDataManager& thn_manager,
                                                 const MultiphaseParameters& params,
                                                 double C,
                                                 double D_u,
                                                 double D_p,
                                                 int coarsest_ln = IBTK::invalid_level_number,
                                                 int finest_ln = IBTK::invalid_level_number);

/*!
 * Accumulate the forces into respective patch indices for the network and solvent
 * on a given patch with variable drag.
 * Assumes ghost cells have been filled for the velocities and volume fraction.
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
                                                 int thn_nc_idx,
                                                 int thn_sc_idx,
                                                 const MultiphaseParameters& params,
                                                 double C,
                                                 double D_u,
                                                 double D_p);
void accumulateMomentumForcesVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
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
                                          double D_p,
                                          int coarsest_ln = IBTK::invalid_level_number,
                                          int finest_ln = IBTK::invalid_level_number);
void accumulateMomentumForcesVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                          int A_un_idx,
                                          int A_us_idx,
                                          int p_idx,
                                          int un_idx,
                                          int us_idx,
                                          const VolumeFractionDataManager& thn_manager,
                                          const MultiphaseParameters& params,
                                          double C,
                                          double D_u,
                                          double D_p,
                                          int coarsest_ln = IBTK::invalid_level_number,
                                          int finest_ln = IBTK::invalid_level_number);

/*!
 * Accumulate the forces into respective patch indices for the network and solvent on a given patch. Assumes ghost cells
 * have been filled for the velocities and volume fraction.
 *
 * Specifically, computes
 *
 * [ C*thn + A_n + D_u*xi   -D_u*xi                ][un]
 * [ -D_u*xi                 C*ths + A_s + D_u*xi  ][us]
 * in which
 * A_i = D_u*eta_i*div(thn*((grad+grad^T)-div*I))
 *
 * These functions result in a runtime error if params.isVariableDrag() returns true.
 */
///\{
/*!
 * Accumulate the momentum forces for constant coefficient problems with the network volume fraction interpolated
 * to cell nodes and cell sides.
 */
void accumulateMomentumWithoutPressureOnPatchConstantCoefficient(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                                 int A_un_idx,
                                                                 int A_us_idx,
                                                                 int un_idx,
                                                                 int us_idx,
                                                                 int thn_idx,
                                                                 int thn_nc_idx,
                                                                 int thn_sc_idx,
                                                                 const MultiphaseParameters& params,
                                                                 double C,
                                                                 double D_u);

void accumulateMomentumWithoutPressureConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                                          const int F_un_idx,
                                                          const int F_us_idx,
                                                          const int un_idx,
                                                          const int us_idx,
                                                          const int thn_idx,
                                                          const int thn_nc_idx,
                                                          const int thn_sc_idx,
                                                          const MultiphaseParameters& params,
                                                          const double C,
                                                          const double D_u,
                                                          int coarsest_ln = IBTK::invalid_index,
                                                          int finest_ln = IBTK::invalid_index);
void accumulateMomentumWithoutPressureConstantCoefficient(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                                          const int F_un_idx,
                                                          const int F_us_idx,
                                                          const int un_idx,
                                                          const int us_idx,
                                                          const VolumeFractionDataManager& thn_manager,
                                                          const MultiphaseParameters& params,
                                                          const double C,
                                                          const double D_u,
                                                          int coarsest_ln = IBTK::invalid_index,
                                                          int finest_ln = IBTK::invalid_index);
///\}

/*!
 * Accumulate the forces into respective patch indices for the network and solvent
 * on a given patch with variable drag.
 * Assumes ghost cells have been filled for the velocities and volume fraction.
 *
 * Specifically, computes
 *
 * [ C*thn + A_n + D_u*xi   -D_u*xi               ][un]
 * [ -D_u*xi                C*ths + A_s + D_u*xi  ][us]
 * in which
 * A_i = D_u*eta_i*div(thn*((grad+grad^T)-div*I))
 *
 * Note that the drag coefficient in this case is not explicitly scaled by the volume fractions.
 *
 * This function results in a runtime error if params.isVariableDrag() returns false.
 *
 * This function linearly interpolates the volume fraction to cell sides and cell nodes as needed.
 */
void accumulateMomentumWithoutPressureOnPatchVariableDrag(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                          const int F_un_idx,
                                                          const int F_us_idx,
                                                          const int un_idx,
                                                          const int us_idx,
                                                          const int thn_idx,
                                                          const int thn_nc_idx,
                                                          const int thn_sc_idx,
                                                          const MultiphaseParameters& params,
                                                          const double C,
                                                          const double D_u);
void accumulateMomentumWithoutPressureVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                                   const int F_un_idx,
                                                   const int F_us_idx,
                                                   const int un_idx,
                                                   const int us_idx,
                                                   const int thn_idx,
                                                   const int thn_nc_idx,
                                                   const int thn_sc_idx,
                                                   const MultiphaseParameters& params,
                                                   const double C,
                                                   const double D_u,
                                                   int coarsest_ln = IBTK::invalid_index,
                                                   int finest_ln = IBTK::invalid_index);
void accumulateMomentumWithoutPressureVariableDrag(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                                   const int F_un_idx,
                                                   const int F_us_idx,
                                                   const int un_idx,
                                                   const int us_idx,
                                                   const VolumeFractionDataManager& thn_manager,
                                                   const MultiphaseParameters& params,
                                                   const double C,
                                                   const double D_u,
                                                   int coarsest_ln = IBTK::invalid_index,
                                                   int finest_ln = IBTK::invalid_index);

/*!
 * Computes the divergence of the volume averaged velocity field on a given patch.
 *
 * Specifically, computes
 *
 * D * Div(thn*un + ths*us)
 *
 * This function linearly interpolates the volume fraction to cell sides as needed.
 */
///{
void applyCoincompressibility(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                              int A_idx,
                              int un_idx,
                              int us_idx,
                              int thn_sc_idx,
                              double D);
void applyCoincompressibility(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                              int A_idx,
                              int un_idx,
                              int us_idx,
                              int thn_sc_idx,
                              double D,
                              int coarsest_ln = IBTK::invalid_level_number,
                              int finest_ln = IBTK::invalid_level_number);
///}

/*!
 * Compute G*[u_n; u_s] with G = C*[thn*grad, ths*grad]^T * p.
 *
 * We assume that thn and p are cell centered and have ghost cells already filled.
 *
 * Gun_idx and Gus_idx should be a side centered quantities.
 *
 * If do_accumulate is true, the gradient is accumulated into Gun_idx and Gus_idx. Otherwise, the values in Gun_idx and
 * Gus_idx are overwritten.
 */
///{
void multiphase_grad_on_hierarchy(SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                                  int Gun_idx,
                                  int Gus_idx,
                                  int thn_sc_idx,
                                  int p_idx,
                                  double C,
                                  bool do_accumulate = true,
                                  int coarsest_ln = IBTK::invalid_level_number,
                                  int finest_ln = IBTK::invalid_level_number);

void multiphase_grad(const SAMRAI::hier::Patch<NDIM>& patch,
                     int Gun_idx,
                     int Gus_idx,
                     int thn_sc_idx,
                     int p_idx,
                     double C,
                     bool do_accumulate = true);
///}

/**
 * Apply the operator
 * [A_n + xi   -xi]
 * [A_n        A_s]
 * to the velocity variables in which
 * A_i = C * thn + D_u*div(eta_i*thn*((grad+grad^T)-lambda_i*div*I))
 */
///{
void computeVelocitySubBlockOnPatch(const SAMRAI::hier::Patch<NDIM>& patch,
                                    const int A_un_idx,
                                    const int A_us_idx,
                                    const int un_idx,
                                    const int us_idx,
                                    const int thn_idx,
                                    const int thn_nc_idx,
                                    const int thn_sc_idx,
                                    const MultiphaseParameters& params,
                                    const double C,
                                    const double D);
///}
} // namespace multiphase

#include "multiphase/private/fd_operators_inc.h"
#endif
