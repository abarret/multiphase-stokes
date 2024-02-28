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

#ifndef included_multiphase_IBMultiphaseHierarchyIntegrator
#define included_multiphase_IBMultiphaseHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBHierarchyIntegrator.h>

#include <ibtk/MarkerPatchHierarchy.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <string>

// Local includes
#include "multiphase/INSVCTwoFluidStaggeredHierarchyIntegrator.h"
#include "multiphase/MultiphaseCrossLinksStrategy.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace multiphase
{
/*!
 * \brief Class IBMultiphaseHierarchyIntegrator is an implementation of a formally
 * second-order accurate, semi-implicit version of the immersed boundary method.
 *
 * <h2>Working with marker points</h2>
 * - Specify the IB kernel in the input database as IB_delta_fcn
 * - Specify the output directory for H5Part data in the input database as
     viz_dump_dirname
 * - Set marker point positions with setMarkerPoints()
 */
class IBMultiphaseHierarchyIntegrator : public IBAMR::IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBMultiphaseHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    IBMultiphaseHierarchyIntegrator(
        std::string object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<IBAMR::IBStrategy> ib_method_ops,
        SAMRAI::tbox::Pointer<IBAMR::IBStrategy> ibn_method_ops,
        SAMRAI::tbox::Pointer<INSVCTwoFluidStaggeredHierarchyIntegrator> ins_hier_integrator,
        bool register_for_restart = true);

    /*!
     * The destructor for class IBMultiphaseHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~IBMultiphaseHierarchyIntegrator() = default;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchySpecialized(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override;

    void initializePatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override;

    void registerCrossLinkStrategy(SAMRAI::tbox::Pointer<MultiphaseCrossLinksStrategy> cross_link_strategy);

protected:
    class IBMultiphaseEulerianForceFunction : public IBTK::CartGridFunction
    {
    public:
        IBMultiphaseEulerianForceFunction(const IBMultiphaseHierarchyIntegrator* ib_solver, int ib_idx);

        bool isTimeDependent() const override;

        void setDataOnPatchHierarchy(const int data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                     const double data_time,
                                     const bool initial_time = false,
                                     const int coarsest_ln = IBTK::invalid_level_number,
                                     const int finest_ln = IBTK::invalid_level_number) override;

        void setDataOnPatch(int data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                            double data_time,
                            bool initial_time = false,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level =
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(NULL)) override;

    private:
        const IBMultiphaseHierarchyIntegrator* const d_ib_solver;
        int d_ib_idx = IBTK::invalid_index;
    };

    friend class IBMultiphaseEulerianForceFunction;

    /*!
     * Perform necessary data movement, workload estimation, and logging prior
     * to regridding.
     */
    void regridHierarchyBeginSpecialized() override;

    /*!
     * Perform necessary data movement and logging after regridding.
     */
    void regridHierarchyEndSpecialized() override;

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    void initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM>> old_level,
                                        bool allocate_data) override;

    void resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                                int coarsest_level,
                                                int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too) override;

    void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                             const int workload_data_idx) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBMultiphaseHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMultiphaseHierarchyIntegrator(const IBMultiphaseHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMultiphaseHierarchyIntegrator& operator=(const IBMultiphaseHierarchyIntegrator& that) = delete;

    void getFromRestart();

    // Network velocity and force.
    // NOTE: Solvent force and velocity are handled by base class.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_un_var, d_fn_var, d_cross_links_un_var,
        d_cross_links_us_var;
    int d_un_idx = IBTK::invalid_index, d_fn_idx = IBTK::invalid_index;
    int d_fn_current_idx = IBTK::invalid_index;
    int d_cross_links_un_idx = IBTK::invalid_index, d_cross_links_current_un_idx = IBTK::invalid_index;
    int d_cross_links_us_idx = IBTK::invalid_index, d_cross_links_current_us_idx = IBTK::invalid_index;

    // Network IB method implementation object
    SAMRAI::tbox::Pointer<IBAMR::IBStrategy> d_ibn_method_ops;

    // IB forcing functions
    SAMRAI::tbox::Pointer<IBMultiphaseEulerianForceFunction> d_fn_fcn, d_fs_fcn;
    SAMRAI::tbox::Pointer<IBMultiphaseEulerianForceFunction> d_cross_links_un_fcn, d_cross_links_us_fcn;
    SAMRAI::tbox::Pointer<MultiphaseCrossLinksStrategy> d_cross_links_strategy;

    // Network communication algorithms.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> d_un_ghostfill_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM>> d_un_ghostfill_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM>> d_un_coarsen_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM>> d_un_coarsen_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> d_fn_prolong_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM>> d_fn_prolong_op;
};
} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_IBMultiphaseHierarchyIntegrator
