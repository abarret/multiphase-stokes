// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_StaggeredStokesOperator
#define included_IBAMR_StaggeredStokesOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/LinearOperator.h"
#include <ibtk/CartGridFunction.h>

#include "IntVector.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "tbox/Pointer.h"

#include <OutersideVariable.h>

#include <string>
#include <vector>

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesOperator is a concrete IBTK::LinearOperator which
 * implements a staggered-grid (MAC) discretization of the incompressible Stokes
 * operator.
 *
 * This class is intended to be used with an iterative (Krylov or Newton-Krylov)
 * incompressible flow solver.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class VCTwoFluidStaggeredStokesOperator : public IBTK::LinearOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    VCTwoFluidStaggeredStokesOperator(const std::string& object_name, bool homogeneous_bc, const double C);

    /*!
     * \brief Destructor.
     */
    ~VCTwoFluidStaggeredStokesOperator();

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    virtual void
    setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& un_problem_coefs,
                                     const SAMRAI::solv::PoissonSpecifications& us_problem_coefs); // might be modified

    /*!
     * \brief Get the PoissonSpecifications object used to specify the
     * coefficients for the network momentum equation in the incompressible Stokes
     * operator.
     */
    virtual const SAMRAI::solv::PoissonSpecifications& getNetworkVelocityPoissonSpecifications() const;

    /*!
     * \brief Get the PoissonSpecifications object used to specify the
     * coefficients for the solvent momentum equation in the incompressible Stokes
     * operator.
     */
    virtual const SAMRAI::solv::PoissonSpecifications& getSolventVelocityPoissonSpecifications() const;

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a un_bc_coefs and \a us_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param un_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the network velocity

     * \param us_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the network velocity

     * \param P_bc_coef   Pointer to object that can set the Robin boundary condition
     *coefficients
     *for the pressure
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                                    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                                    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    /*!
     * \brief Set the physical boundary condition helper object.
     */
    virtual void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper);

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * \note In general, the vectors x and y \em cannot be the same.
     *
     * Upon return from this function, the y vector will contain the result of
     * the application of A to x.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */

    void setThnIdx(int thn_idx);
    // brief This method sets up the patch data indices for Thn variable

    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Compute hierarchy dependent data required for computing y=Ax and
     * z=Ax+y.
     *
     * The vector arguments for apply(), applyAdjoint(), etc, need not match
     * those for initializeOperatorState().  However, there must be a certain
     * degree of similarity, including
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the input and output vectors
     *
     * \note It is generally necessary to reinitialize the operator state when
     * the hierarchy configuration changes.
     *
     * It is safe to call initializeOperatorState() when the state is already
     * initialized.  In this case, the operator state is first deallocated and
     * then reinitialized.
     *
     * Conditions on arguments:
     * - input and output vectors must have same hierarchy
     * - input and output vectors must have same structure, depth, etc.
     *
     * Call deallocateOperatorState() to remove any data allocated by this
     * method.
     *
     * \see deallocateOperatorState
     *
     * \param in input vector
     * \param out output vector
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * Remove all hierarchy dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() when the operator state is
     * already deallocated.
     *
     * \see initializeOperatorState
     */
    void deallocateOperatorState() override;

    /*!
     * \brief Modify the RHS vector to account for physical boundary conditions.
     */
    void modifyRhsForBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Modify the solution vector to account for physical boundary conditions.
     */
    void imposeSolBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u) override;

    //\}

protected:
    // Problem specification.
    SAMRAI::solv::PoissonSpecifications d_un_problem_coefs; // might create a different set for second fluid equation
    SAMRAI::solv::PoissonSpecifications d_us_problem_coefs; // might create a different set for second fluid equation
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_default_un_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_default_us_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_un_bc_coefs;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_us_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_default_P_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;
    double d_C = std::numeric_limits<double>::quiet_NaN(); // C*u

    // Reference patch hierarchy
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    // Boundary condition helper object.
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    // Cached communications operators.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> d_un_fill_pattern, d_us_fill_pattern,
        d_P_fill_pattern;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill, d_no_fill;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_x, d_b;
    int d_thn_idx;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCTwoFluidStaggeredStokesOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCTwoFluidStaggeredStokesOperator(const VCTwoFluidStaggeredStokesOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCTwoFluidStaggeredStokesOperator& operator=(const VCTwoFluidStaggeredStokesOperator& that) = delete;

    // Synchronization variable
    SAMRAI::tbox::Pointer<SAMRAI::pdat::OutersideVariable<NDIM, double>> d_os_var;
    int d_os_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM>> d_os_coarsen_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM>> d_os_coarsen_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM>>> d_os_coarsen_scheds;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_TwoFluidStaggeredStokesOperator
