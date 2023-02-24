// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/INSStaggeredConvectiveOperatorManager.h"
#include "ibamr/INSStaggeredPressureBcCoef.h"
#include "ibamr/INSStaggeredVelocityBcCoef.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellDoubleBoundsPreservingConservativeLinearRefine.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/CartSideDoubleSpecializedLinearRefine.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/ibtk_utilities.h"
#include <ibtk/PETScKrylovLinearSolver.h>

#include "ArrayData.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// Local includes
#include "VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.h"
#include "VCTwoFluidStaggeredStokesOperator.h"
#include "tests/multigrid/FullFACPreconditioner.h"
#include "time_stepping/INSVCTwoFluidStaggeredHierarchyIntegrator.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
INSVCTwoFluidStaggeredHierarchyIntegrator::INSVCTwoFluidStaggeredHierarchyIntegrator(
    std::string object_name,
    Pointer<Database> input_db,
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry,
    bool register_for_restart)
    : INSHierarchyIntegrator(std::move(object_name),
                             input_db,
                             new SideVariable<NDIM, double>(object_name + "::Us"),
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             new CellVariable<NDIM, double>(object_name + "::P"),
                             "CONSERVATIVE_COARSEN",
                             "LINEAR_REFINE",
                             new SideVariable<NDIM, double>(object_name + "::F"),
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             new CellVariable<NDIM, double>(object_name + "::Q"),
                             "CONSERVATIVE_COARSEN",
                             "CONSTANT_REFINE",
                             register_for_restart),
      d_thn_fcn("thn", input_db->getDatabase("thn"), grid_geometry),
      d_f_un_fcn("f_un", input_db->getDatabase("f_un"), grid_geometry),
      d_f_us_fcn("f_us", input_db->getDatabase("f_us"), grid_geometry),
      d_f_p_fcn("f_p", input_db->getDatabase("f_p"), grid_geometry)
{
    if (input_db->keyExists("rho")) d_rho = input_db->getDouble("rho");
    if (input_db->keyExists("solver")) d_solver_db = input_db->getDatabase("solver");
    if (input_db->keyExists("precond")) d_precond_db = input_db->getDatabase("precond");
    if (input_db->keyExists("w")) d_w = input_db->getDouble("w");
    if (input_db->keyExists("max_multigrid_levels"))
        d_max_multigrid_levels = input_db->getInteger("max_multigrid_levels");
    return;
} // INSVCTwoFluidStaggeredHierarchyIntegrator

INSVCTwoFluidStaggeredHierarchyIntegrator::~INSVCTwoFluidStaggeredHierarchyIntegrator()
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete d_U_bc_coefs[d];
        d_U_bc_coefs[d] = nullptr;
    }
    delete d_P_bc_coef;
    d_P_bc_coef = nullptr;
    d_velocity_solver.setNull();
    d_pressure_solver.setNull();
} // ~INSVCTwoFluidStaggeredHierarchyIntegrator

Pointer<ConvectiveOperator>
INSVCTwoFluidStaggeredHierarchyIntegrator::getConvectiveOperator()
{
    return nullptr;
} // getConvectiveOperator

Pointer<PoissonSolver>
INSVCTwoFluidStaggeredHierarchyIntegrator::getVelocitySubdomainSolver()
{
    return nullptr;
} // getVelocitySubdomainSolver

Pointer<PoissonSolver>
INSVCTwoFluidStaggeredHierarchyIntegrator::getPressureSubdomainSolver()
{
    return nullptr;
} // getPressureSubdomainSolver

void
INSVCTwoFluidStaggeredHierarchyIntegrator::setInitialData(Pointer<CartGridFunction> un_fcn,
                                                          Pointer<CartGridFunction> us_fcn,
                                                          Pointer<CartGridFunction> p_fcn)
{
    d_un_init_fcn = un_fcn;
    d_us_init_fcn = us_fcn;
    d_p_init_fcn = p_fcn;
    return;
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                                         Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Here we do all we need to ensure that calls to advanceHierarchy() or integrateHierarchy() are valid.
    // NOTE: This function is called before the patch hierarchy has valid patch levels.
    // To set initial data, we should do this in initializeLevelDataSpecialized().

    // First create the variables we need.
    d_un_sc_var = new SideVariable<NDIM, double>(d_object_name + "un_sc");
    d_us_sc_var = new SideVariable<NDIM, double>(d_object_name + "us_sc");
    d_p_cc_var = new CellVariable<NDIM, double>(d_object_name + "p_cc");
    d_thn_cc_var = new CellVariable<NDIM, double>(d_object_name + "thn_cc");
    d_f_un_sc_var = new SideVariable<NDIM, double>(d_object_name + "f_un_sc");
    d_f_us_sc_var = new SideVariable<NDIM, double>(d_object_name + "f_us_sc");
    d_f_cc_var = new CellVariable<NDIM, double>(d_object_name + "f_cc");

    // Now register patch indices to store data.
    // Note: We can retreive the data using the variable database along with the member variable from above and the
    // member context.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_ctx = var_db->getContext(d_object_name + "::context");
    var_db->registerVariableAndContext(d_un_sc_var, d_ctx, IntVector<NDIM>(1));
    var_db->registerVariableAndContext(d_us_sc_var, d_ctx, IntVector<NDIM>(1));
    var_db->registerVariableAndContext(d_p_cc_var, d_ctx, IntVector<NDIM>(1));
    var_db->registerVariableAndContext(d_thn_cc_var, d_ctx, IntVector<NDIM>(1));
    var_db->registerVariableAndContext(d_f_cc_var, d_ctx, IntVector<NDIM>(1));
    var_db->registerVariableAndContext(d_f_un_sc_var, d_ctx, IntVector<NDIM>(1));
    var_db->registerVariableAndContext(d_f_us_sc_var, d_ctx, IntVector<NDIM>(1));

    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
INSVCTwoFluidStaggeredHierarchyIntegrator::initializeLevelDataSpecialized(Pointer<BasePatchHierarchy<NDIM>> hierarchy,
                                                                          int level_number,
                                                                          double init_data_time,
                                                                          bool can_be_refined,
                                                                          bool initial_time,
                                                                          Pointer<BasePatchLevel<NDIM>> old_level,
                                                                          bool allocate_data)
{
    // Here is where we set initial data for our state variables.
    // These are solvent velocity, network velocity, and pressure.
    // NOTE: Pressure isn't really a "state" variable, but we let the data persist for initial guesses.
    if (initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int un_sc_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, d_ctx);
        const int us_sc_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, d_ctx);
        const int p_cc_idx = var_db->mapVariableAndContextToIndex(d_p_cc_var, d_ctx);
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
        if (!level->checkAllocated(un_sc_idx)) level->allocatePatchData(un_sc_idx, init_data_time);
        if (!level->checkAllocated(us_sc_idx)) level->allocatePatchData(us_sc_idx, init_data_time);
        if (!level->checkAllocated(p_cc_idx)) level->allocatePatchData(p_cc_idx, init_data_time);

        // Now let's set the initial data.
        d_un_init_fcn->setDataOnPatchLevel(un_sc_idx, d_un_sc_var, level, init_data_time, true);
        d_us_init_fcn->setDataOnPatchLevel(us_sc_idx, d_us_sc_var, level, init_data_time, true);
        d_p_init_fcn->setDataOnPatchLevel(p_cc_idx, d_p_cc_var, level, init_data_time, true);
    }

    // TODO: When a regrid operation occurs, we will have to copy data from the "old level" to the new hierarchy.
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                        const double new_time,
                                                                        const int num_cycles)
{
    INSHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Do anything that needs to be done before we call integrateHierarchy().

    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
INSVCTwoFluidStaggeredHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                              const double new_time,
                                                              const int cycle_num)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_sc_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, d_ctx);
    const int us_sc_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, d_ctx);
    const int p_cc_idx = var_db->mapVariableAndContextToIndex(d_p_cc_var, d_ctx);
    const int f_un_sc_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, d_ctx);
    const int f_us_sc_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, d_ctx);
    const int f_cc_idx = var_db->mapVariableAndContextToIndex(d_f_cc_var, d_ctx);
    const int thn_cc_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, d_ctx);

    // Allocate scratch data on each level of the patch hierarchy.
    // NOTE: State data (velocity and pressure) has already been allocated.
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(f_un_sc_idx, 0.0);
        level->allocatePatchData(f_us_sc_idx, 0.0);
        level->allocatePatchData(f_cc_idx, 0.0);
        level->allocatePatchData(thn_cc_idx, 0.0);
    }

    // Setup un(n), us(n), p(n) and right-hand-side vectors.
    HierarchyMathOps hier_math_ops("hier_math_ops", d_hierarchy);
    const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
    const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

    // not expensive to create these vectors each time
    SAMRAIVectorReal<NDIM, double> u_vec("u", d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
    SAMRAIVectorReal<NDIM, double> f_vec("f", d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());

    u_vec.addComponent(d_un_sc_var, un_sc_idx, h_sc_idx);
    u_vec.addComponent(d_us_sc_var, us_sc_idx, h_sc_idx);
    u_vec.addComponent(d_p_cc_var, p_cc_idx, h_cc_idx);
    f_vec.addComponent(d_f_un_sc_var, f_un_sc_idx, h_sc_idx);
    f_vec.addComponent(d_f_us_sc_var, f_us_sc_idx, h_sc_idx);
    f_vec.addComponent(d_f_cc_var, f_cc_idx, h_cc_idx);

    Pointer<SAMRAIVectorReal<NDIM, double>> u_new_vec;
    u_new_vec = u_vec.cloneVector("u_new");      // should delete the vector at the end
    u_new_vec->allocateVectorData();
    Pointer<SAMRAIVectorReal<NDIM, double>> u_vec_ptr(&u_vec, false);
    u_new_vec->copyVector(u_vec_ptr, true);

    // Only pressure nullspace
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> null_vecs(1, nullptr);
    null_vecs[0] = u_vec.cloneVector("PressureNull"); // should delete the vector at the end
    null_vecs[0]->allocateVectorData();
    null_vecs[0]->setToScalar(0.0);
    // Pull out pressure component and set to constant
    {
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> p_data = null_vecs[0]->getComponentPatchData(2, *patch);
                p_data->fillAll(1.0);
            }
        }
    }

    d_f_un_fcn.setDataOnPatchHierarchy(f_un_sc_idx, d_f_un_sc_var, d_hierarchy, 0.0);
    d_f_us_fcn.setDataOnPatchHierarchy(f_us_sc_idx, d_f_us_sc_var, d_hierarchy, 0.0);
    d_f_p_fcn.setDataOnPatchHierarchy(f_cc_idx, d_f_cc_var, d_hierarchy, 0.0);
    d_thn_fcn.setDataOnPatchHierarchy(thn_cc_idx, d_thn_cc_var, d_hierarchy, 0.0);

    const double dt = new_time - current_time;
    const double C = d_rho/dt;

    // set-up RHS for backward Euler scheme: f(n) + C*u_i(n) for  i = n, s
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_sc_idx);
            Pointer<SideData<NDIM, double>> us_data = patch->getPatchData(us_sc_idx);
            Pointer<SideData<NDIM, double>> f_un_data = patch->getPatchData(f_un_sc_idx);
            Pointer<SideData<NDIM, double>> f_us_data = patch->getPatchData(f_us_sc_idx);
            Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(thn_cc_idx);
            IntVector<NDIM> xp(1, 0), yp(0, 1);

            for (SideIterator<NDIM> si(patch->getBox(), 0); si; si++) // side-centers in x-dir
            {
                const SideIndex<NDIM>& idx = si();                                         // axis = 0, (i-1/2,j)
                CellIndex<NDIM> idx_c_low = idx.toCell(0);                                 // (i-1,j)
                CellIndex<NDIM> idx_c_up = idx.toCell(1);                                  // (i,j)
                double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i-1/2,j)
                (*f_un_data)(idx) = (*f_un_data)(idx) + C * thn_lower * (*un_data)(idx);
                (*f_us_data)(idx) = (*f_us_data)(idx) + C * (1.0 - thn_lower) * (*us_data)(idx);
            }

            for (SideIterator<NDIM> si(patch->getBox(), 1); si; si++) // side-centers in y-dir
            {
                const SideIndex<NDIM>& idx = si();                                         // axis = 1, (i,j-1/2)
                CellIndex<NDIM> idx_c_low = idx.toCell(0);                                 // (i,j-1)
                CellIndex<NDIM> idx_c_up = idx.toCell(1);                                  // (i,j)
                double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i,j-1/2)
                (*f_un_data)(idx) = (*f_un_data)(idx) + C * thn_lower * (*un_data)(idx);
                (*f_us_data)(idx) = (*f_us_data)(idx) + C * (1.0 - thn_lower) * (*us_data)(idx);
            }
        } // patches
    }     // levels

    // Setup the stokes operator
    Pointer<VCTwoFluidStaggeredStokesOperator> stokes_op = new VCTwoFluidStaggeredStokesOperator("stokes_op", true);
    const double D = -1.0; // This will depend on the time stepping scheme.
    stokes_op->setCandDCoefficients(C, D);
    stokes_op->setDragCoefficient(1.0, 1.0, 1.0);
    stokes_op->setViscosityCoefficient(1.0, 1.0);
    stokes_op->setThnIdx(thn_cc_idx);

    Pointer<PETScKrylovLinearSolver> krylov_solver = new PETScKrylovLinearSolver("solver", d_solver_db, "solver_");
    krylov_solver->setOperator(stokes_op);

    // Now create a preconditioner
    Pointer<VCTwoFluidStaggeredStokesBoxRelaxationFACOperator> fac_precondition_strategy =
        new VCTwoFluidStaggeredStokesBoxRelaxationFACOperator("KrylovPrecondStrategy", "Krylov_precond_", d_w, C);
    fac_precondition_strategy->setThnIdx(thn_cc_idx);
    Pointer<FullFACPreconditioner> Krylov_precond = new FullFACPreconditioner(
        "KrylovPrecond", fac_precondition_strategy, d_precond_db, d_max_multigrid_levels, "Krylov_precond_");
    Krylov_precond->setNullspace(false, null_vecs);
    krylov_solver->setPreconditioner(Krylov_precond);

    krylov_solver->setNullspace(false, null_vecs);
    krylov_solver->initializeSolverState(u_vec, f_vec);

    // Set thn_cc_idx on the dense hierarchy.
    Pointer<PatchHierarchy<NDIM>> dense_hierarchy = Krylov_precond->getDenseHierarchy();
    // Allocate data
    for (int ln = 0; ln <= dense_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = dense_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(thn_cc_idx)) level->allocatePatchData(thn_cc_idx, 0.0);
    }
    d_thn_fcn.setDataOnPatchHierarchy(
        thn_cc_idx, d_thn_cc_var, dense_hierarchy, 0.0, false, 0, dense_hierarchy->getFinestLevelNumber());
    // Also fill in theta ghost cells
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(1);
    ghost_cell_comp[0] = ITC(thn_cc_idx,
                             "CONSERVATIVE_LINEAR_REFINE",
                             false,
                             "NONE",
                             "LINEAR",
                             true,
                             nullptr); // defaults to fill corner
    HierarchyGhostCellInterpolation ghost_cell_fill;
    ghost_cell_fill.initializeOperatorState(
        ghost_cell_comp, dense_hierarchy, 0, dense_hierarchy->getFinestLevelNumber());
    ghost_cell_fill.fillData(0.0);

    // Solve for un(n+1), us(n+1), p(n+1).
    krylov_solver->solveSystem(*u_new_vec, f_vec);

    // Reset the solution vector: u^n = u^(n+1)
    u_vec.copyVector(u_new_vec, true);

    // Free storage from cloning vectors
    u_new_vec->freeVectorComponents();
    null_vecs[0]->freeVectorComponents();

    // Deallocate patch data
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(f_un_sc_idx);
        level->deallocatePatchData(f_us_sc_idx);
        level->deallocatePatchData(f_cc_idx);
        level->deallocatePatchData(thn_cc_idx);
    }

    // Deallocate patch data on the dense hierarchy
    for (int ln = 0; ln <= dense_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = dense_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(thn_cc_idx)) level->deallocatePatchData(thn_cc_idx);
    }

    return;
} // integrateHierarchy

void
INSVCTwoFluidStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                         const double new_time,
                                                                         const bool skip_synchronize_new_state_data,
                                                                         const int num_cycles)
{
    INSHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Do anything that needs to be done after integrateHierarchy().

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////