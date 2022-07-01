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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/FACPreconditioner.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "FullFACPreconditioner.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <ostream>
#include <string>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FullFACPreconditioner::FullFACPreconditioner(std::string object_name,
                                             Pointer<FACPreconditionerStrategy> fac_strategy,
                                             Pointer<Database> input_db,
                                             Pointer<GriddingAlgorithm<NDIM>> grid_alg,
                                             const std::string& default_options_prefix)
    : FACPreconditioner(std::move(object_name), fac_strategy, input_db, default_options_prefix), d_grid_alg(grid_alg)
{
    // intentionally blank
    return;
} // FACPreconditioner

FullFACPreconditioner::~FullFACPreconditioner()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~FACPreconditioner

bool
FullFACPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& u, SAMRAIVectorReal<NDIM, double>& f)
{
    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(u, f);

    // Set the initial guess to equal zero.
    u.setToScalar(0.0, /*interior_only*/ false);

    // We need to copy f into d_f
    transferToDense(f, *d_f);

    // Apply a single FAC cycle.
    if (d_cycle_type == V_CYCLE && d_num_pre_sweeps == 0)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(!d_f);
        TBOX_ASSERT(!d_r);
#endif
        // V-cycle MG without presmoothing keeps the residual equal to the
        // initial right-hand-side vector f, so we can simply use that vector
        // for the residual in the FAC algorithm.
        FACVCycleNoPreSmoothing(*d_u, *d_f, d_finest_ln);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_f);
        TBOX_ASSERT(d_r);
#endif
        d_r->copyVector(d_f, false);
        switch (d_cycle_type)
        {
        case F_CYCLE:
            FCycle(*d_u, *d_f, d_finest_ln);
            break;
        case FMG_CYCLE:
            FMGCycle(*d_u, *d_f, d_finest_ln, 1);
            break;
        case V_CYCLE:
            muCycle(*d_u, *d_f, d_finest_ln, 1);
            break;
        case W_CYCLE:
            muCycle(*d_u, *d_f, d_finest_ln, 2);
            break;
        default:
            TBOX_ERROR(d_object_name << "::solveSystem():\n"
                                     << "  unsupported FAC cycle type: " << enum_to_string<MGCycleType>(d_cycle_type)
                                     << "." << std::endl);
        }
    }

    transferToBase(u, *d_u);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();
    return true;
} // solveSystem

void
FullFACPreconditioner::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& solution,
                                             const SAMRAIVectorReal<NDIM, double>& rhs)
{
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        deallocateSolverState();
    }

    // Setup operator state.
    d_hierarchy = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln = solution.getFinestLevelNumber();
    int num_comps = solution.getNumberOfComponents();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == rhs.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == rhs.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == rhs.getFinestLevelNumber());
    TBOX_ASSERT(num_comps == rhs.getNumberOfComponents());
#endif

    generateDenseHierarchy(d_hierarchy);
    d_coarsest_ln = 0;
    d_finest_ln = d_dense_hierarchy->getFinestLevelNumber();

    // We need to allocate vectors on the dense hierarchy
    d_u = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::u", d_dense_hierarchy, 0, d_dense_hierarchy->getFinestLevelNumber());
    d_f = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::f", d_dense_hierarchy, 0, d_dense_hierarchy->getFinestLevelNumber());

    // We can use the same patch indicies in u and f as in solution and rhs.
    for (int comp_num = 0; comp_num < num_comps; ++comp_num)
    {
        d_u->addComponent(solution.getComponentVariable(comp_num),
                          solution.getComponentDescriptorIndex(comp_num),
                          solution.getControlVolumeIndex(comp_num));
        d_f->addComponent(rhs.getComponentVariable(comp_num),
                          rhs.getComponentDescriptorIndex(comp_num),
                          rhs.getControlVolumeIndex(comp_num));
    }

    d_u->allocateVectorData();
    d_f->allocateVectorData();

    d_fac_strategy->initializeOperatorState(*d_u, *d_f);

    // Create temporary vectors.
    if (!(d_cycle_type == V_CYCLE && d_num_pre_sweeps == 0))
    {
        d_r = d_f->cloneVector("");
        d_r->allocateVectorData();
    }

    // Allocate scratch data.
    d_fac_strategy->allocateScratchData();

    // Indicate the operator is initialized.
    d_is_initialized = true;
    return;
} // initializeSolverState

void
FullFACPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Deallocate scratch data.
    d_fac_strategy->deallocateScratchData();

    // Destroy temporary vectors.
    if (d_f)
    {
        d_f->deallocateVectorData();
        d_f.setNull();
    }

    if (d_r)
    {
        d_r->deallocateVectorData();
        d_r->freeVectorComponents();
        d_r.setNull();
    }

    if (d_u)
    {
        d_u->deallocateVectorData();
        d_u.setNull();
    }

    // Deallocate operator state.
    d_fac_strategy->deallocateOperatorState();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;
    return;
} // deallocateSolverState

void
FullFACPreconditioner::generateDenseHierarchy(Pointer<PatchHierarchy<NDIM>> base_hierarchy)
{
    // Build our new grid hierarchy.
    d_dense_hierarchy =
        new PatchHierarchy<NDIM>(d_object_name + "::DenseHierarchy", base_hierarchy->getGridGeometry(), false);
    int max_base_levels = base_hierarchy->getFinestLevelNumber();
    int max_dense_levels = d_grid_alg->getMaxLevels();
    int num_levels_not_in_base = max_dense_levels - max_base_levels;

    d_grid_alg->makeCoarsestLevel(d_dense_hierarchy, 0.0);
    // Now make levels up to the coarse level of the original hierarchy
    int ln = 0;
    for (; ln < num_levels_not_in_base; ++ln) d_grid_alg->makeFinerLevel(d_dense_hierarchy, 0.0, 0.0, /*tag_buffer*/ 1);
    // Now copy the levels from the base hierarchy to the dense hierarchy
    for (; ln <= max_dense_levels; ++ln)
    {
        Pointer<PatchLevel<NDIM>> base_level = base_hierarchy->getPatchLevel(ln - num_levels_not_in_base);
        d_dense_hierarchy->makeNewPatchLevel(ln,
                                             base_level->getRatio() * std::pow(2, num_levels_not_in_base),
                                             base_level->getBoxes(),
                                             base_level->getProcessorMapping());
    }
}

void
FullFACPreconditioner::transferToDense(const SAMRAIVectorReal<NDIM, double>& base_x,
                                       const SAMRAIVectorReal<NDIM, double>& dense_x)
{
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    // We need to create refinement schedules to go from the base to dense hierarchy
    int level_diff = d_dense_hierarchy->getFinestLevelNumber() - d_hierarchy->getFinestLevelNumber();
    Pointer<RefineAlgorithm<NDIM>> refine_alg = new RefineAlgorithm<NDIM>();
    std::vector<Pointer<RefineSchedule<NDIM>>> refine_scheds(d_hierarchy->getFinestLevelNumber());
    for (int ln = 0; ln < d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> base_level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM>> dense_level = d_dense_hierarchy->getPatchLevel(ln + level_diff);
        for (int comp = 0; comp < dense_x.getNumberOfComponents(); ++comp)
        {
            Pointer<RefineOperator<NDIM>> refine_op =
                grid_geom->lookupRefineOperator(dense_x.getComponentVariable(comp), "CONSERVATIVE_LINEAR_REFINE");
            int dense_idx = dense_x.getComponentDescriptorIndex(comp);
            int base_idx = base_x.getComponentDescriptorIndex(comp);
            refine_alg->registerRefine(dense_idx, base_idx, base_idx, refine_op);
        }
        refine_scheds[ln] = refine_alg->createSchedule(dense_level, base_level);
        refine_scheds[ln]->fillData(0.0, false);
    }

    // Now we need to coarsen the data from the dense hierarchy to the dense hierarchy.
    Pointer<CoarsenAlgorithm<NDIM>> coarsen_alg = new CoarsenAlgorithm<NDIM>();
    std::vector<Pointer<CoarsenSchedule<NDIM>>> coarsen_scheds(level_diff);
    for (int ln = level_diff; ln >= 0; --ln)
    {
        Pointer<PatchLevel<NDIM>> dst_level = d_dense_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM>> src_level = d_dense_hierarchy->getPatchLevel(ln + 1);
        for (int comp = 0; comp < dense_x.getNumberOfComponents(); ++comp)
        {
            Pointer<CoarsenOperator<NDIM>> coarsen_op =
                grid_geom->lookupCoarsenOperator(dense_x.getComponentVariable(comp), "CUBIC_COARSEN");
            int dense_idx = dense_x.getComponentDescriptorIndex(comp);
            coarsen_alg->registerCoarsen(dense_idx, dense_idx, coarsen_op);
        }
        coarsen_scheds[ln] = coarsen_alg->createSchedule(dst_level, src_level);
        coarsen_scheds[ln]->coarsenData();
    }
}

void
FullFACPreconditioner::transferToBase(const SAMRAIVectorReal<NDIM, double>& base_x,
                                      const SAMRAIVectorReal<NDIM, double>& dense_x)
{
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    // We need to create a refinement schedule from the dense hierarchy to the base hierarchy.
    Pointer<RefineAlgorithm<NDIM>> refine_alg = new RefineAlgorithm<NDIM>();
    std::vector<Pointer<RefineSchedule<NDIM>>> refine_scheds(d_hierarchy->getFinestLevelNumber());
    int level_diff = d_dense_hierarchy->getFinestLevelNumber() - d_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln < d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> base_level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM>> dense_level = d_dense_hierarchy->getPatchLevel(ln + level_diff);
        for (int comp = 0; comp < dense_x.getNumberOfComponents(); ++comp)
        {
            Pointer<RefineOperator<NDIM>> refine_op =
                grid_geom->lookupRefineOperator(dense_x.getComponentVariable(comp), "CONSERVATIVE_LINEAR_REFINE");
            int dense_idx = dense_x.getComponentDescriptorIndex(comp);
            int base_idx = base_x.getComponentDescriptorIndex(comp);
            refine_alg->registerRefine(base_idx, dense_idx, dense_idx, refine_op);
        }
        refine_scheds[ln] = refine_alg->createSchedule(base_level, dense_level);
        refine_scheds[ln]->fillData(0.0, false);
    }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
