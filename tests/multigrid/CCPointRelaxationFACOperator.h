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

#ifndef included_CCPointRelaxationFACOperator
#define included_CCPointRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/PoissonFACPreconditioner.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"

#include "IntVector.h"
#include "PoissonSpecifications.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <string>
#include <vector>

// Local includes
#include "CCPointFACPreconditionerStrategy.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
template <int DIM>
class BoxList;
template <int DIM>
class Patch;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class CellData;
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CCPointRelaxationFACOperator is a concrete
 * PoissonFACPreconditionerStrategy for solving elliptic equations of the form
 * \f$ \mbox{$L u$} = \mbox{$(C I + \nabla \cdot D \nabla) u$} = f \f$ using a
 * globally second-order accurate cell-centered finite-volume discretization,
 * with support for Robin and periodic boundary conditions.
 *
 * This class provides operators that are used by class FACPreconditioner to
 * solve scalar Poisson-type equations of the form \f[ (C I + \nabla \cdot D
 * \nabla) u = f \f] using a cell-centered, globally second-order accurate
 * finite-volume discretization, where
 *
 * - \f$ C \f$, \f$ D \f$ and \f$ f \f$ are independent of \f$ u \f$,
 * - \f$ C \f$ is a \em scalar,
 * - \f$ D \f$ is a \em scalar diffusion coefficient, and
 * - \f$ f \f$ is a cell-centered scalar function.
 *
 * Robin boundary conditions may be specified at physical boundaries; see class
 * SAMRAI::solv::RobinBcCoefStrategy.
 *
 * By default, the class is configured to solve the Poisson problem \f$
 * -\nabla^2 u = f \f$, subject to homogeneous Dirichlet boundary conditions.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_type = "PATCH_GAUSS_SEIDEL"         // see setSmootherType()
 prolongation_method = "LINEAR_REFINE"        // see setProlongationMethod()
 restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethod()
 coarse_solver_type = "HYPRE_LEVEL_SOLVER"    // see setCoarseSolverType()
 coarse_solver_rel_residual_tol = 1.0e-5      // see setCoarseSolverRelativeTolerance()
 coarse_solver_abs_residual_tol = 1.0e-50     // see setCoarseSolverAbsoluteTolerance()
 coarse_solver_max_iterations = 1             // see setCoarseSolverMaxIterations()
 coarse_solver_db {                           // SAMRAI::tbox::Database for initializing coarse
 level solver
    solver_type = "PFMG"
    num_pre_relax_steps = 0
    num_post_relax_steps = 2
 }
 \endverbatim
*/
class CCPointRelaxationFACOperator : public CCPointFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    CCPointRelaxationFACOperator(const std::string& object_name,
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                 std::string default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~CCPointRelaxationFACOperator();

    /*!
     * \name Functions for configuring the solver.
     */
    //\{

    /*!
     * \brief Specify the smoother type.
     *
     * Select from:
     * - \c "PATCH_GAUSS_SEIDEL"
     * - \c "PROCESSOR_GAUSS_SEIDEL"
     * - \c "RED_BLACK_GAUSS_SEIDEL"
     */
    void setSmootherType(const std::string& smoother_type) override;

    /*!
     * \brief Specify the coarse level solver.
     */
    void setCoarseSolverType(const std::string& coarse_solver_type) override;

    //\}

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are
     *being
     *performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are
     *being
     *performed
     */
    void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps) override;

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                            int coarsest_ln) override;

    /*!
     * \brief Compute composite grid residual on a range of levels.
     *
     * \param residual residual vector
     * \param solution solution vector
     * \param rhs source (right hand side) vector
     * \param coarsest_level_num coarsest level number
     * \param finest_level_num finest level number
     */
    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num) override;

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                                            int coarsest_reset_ln,
                                            int finest_reset_ln) override;

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void deallocateOperatorStateSpecialized(int coarsest_reset_ln, int finest_reset_ln) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CCPointRelaxationFACOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CCPointRelaxationFACOperator(const CCPointRelaxationFACOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCPointRelaxationFACOperator& operator=(const CCPointRelaxationFACOperator& that) = delete;

    /*
     * Coarse level solvers and solver parameters.
     */
    SAMRAI::tbox::Pointer<PoissonSolver> d_coarse_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_coarse_solver_db;

    /*
     * Patch overlap data.
     */
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM>>> d_patch_bc_box_overlap;
    std::vector<std::vector<std::map<int, SAMRAI::hier::Box<NDIM>>>> d_patch_neighbor_overlap;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CCPointRelaxationFACOperator
