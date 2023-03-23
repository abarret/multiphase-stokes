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
#include "FullFACPreconditioner.h"
#include "INSVCTwoFluidStaggeredHierarchyIntegrator.h"
#include "VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.h"
#include "VCTwoFluidStaggeredStokesOperator.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
static double
convertToThs(const double thn)
{
    return 1.0 - thn;
}

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
                             register_for_restart)
{
    if (input_db->keyExists("viscous_time_stepping_type"))
        d_viscous_time_stepping_type =
            string_to_enum<TimeSteppingType>(input_db->getString("viscous_time_stepping_type"));
    if (input_db->keyExists("rho")) d_rho = input_db->getDouble("rho");
    if (input_db->keyExists("solver_db")) d_solver_db = input_db->getDatabase("solver_db");
    if (input_db->keyExists("precond_db")) d_precond_db = input_db->getDatabase("precond_db");
    if (input_db->keyExists("w")) d_w = input_db->getDouble("w");
    if (input_db->keyExists("use_preconditioner")) d_use_preconditioner = input_db->getBool("use_preconditioner");
    return;
} // INSVCTwoFluidStaggeredHierarchyIntegrator

INSVCTwoFluidStaggeredHierarchyIntegrator::~INSVCTwoFluidStaggeredHierarchyIntegrator()
{
    // intentionally blank
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

Pointer<SideVariable<NDIM, double>>
INSVCTwoFluidStaggeredHierarchyIntegrator::getSolventVariable() const
{
    return d_us_sc_var;
}

Pointer<SideVariable<NDIM, double>>
INSVCTwoFluidStaggeredHierarchyIntegrator::getNetworkVariable() const
{
    return d_un_sc_var;
}

Pointer<CellVariable<NDIM, double>>
INSVCTwoFluidStaggeredHierarchyIntegrator::getPressureVariable() const
{
    return d_P_var;
}

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
INSVCTwoFluidStaggeredHierarchyIntegrator::setForcingFunctions(Pointer<CartGridFunction> fn_fcn,
                                                               Pointer<CartGridFunction> fs_fcn,
                                                               Pointer<CartGridFunction> fp_fcn)
{
    d_f_un_fcn = fn_fcn;
    d_f_us_fcn = fs_fcn;
    d_f_p_fcn = fp_fcn;
    return;
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::setNetworkVolumeFractionFunction(Pointer<CartGridFunction> thn_fcn)
{
    // Make sure we have a valid pointer.
    TBOX_ASSERT(thn_fcn);
    d_thn_fcn = thn_fcn;
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
    // NOTE: d_P_var is a member variable of the base class.
    d_un_sc_var = new SideVariable<NDIM, double>(d_object_name + "un_sc");
    d_us_sc_var = new SideVariable<NDIM, double>(d_object_name + "us_sc");
    d_thn_cc_var = new CellVariable<NDIM, double>(d_object_name + "thn_cc");
    d_f_un_sc_var = new SideVariable<NDIM, double>(d_object_name + "f_un_sc");
    d_f_us_sc_var = new SideVariable<NDIM, double>(d_object_name + "f_us_sc");
    d_f_cc_var = new CellVariable<NDIM, double>(d_object_name + "f_cc");

    // Register variables with the integrator. Those with states (Velocities) have associated current, new, and scratch
    // indices. NOTE: Pressure is NOT a state variable, but we keep track of current and new values for initial guesses
    // for the solver.
    int un_cur_idx, un_scr_idx, un_new_idx;
    int us_cur_idx, us_scr_idx, us_new_idx;
    int p_cur_idx, p_scr_idx, p_new_idx;
    registerVariable(un_cur_idx,
                     un_new_idx,
                     un_scr_idx,
                     d_un_sc_var,
                     IntVector<NDIM>(1),
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_un_init_fcn);
    registerVariable(us_cur_idx,
                     us_new_idx,
                     us_scr_idx,
                     d_us_sc_var,
                     IntVector<NDIM>(1),
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_us_init_fcn);
    registerVariable(p_cur_idx,
                     p_new_idx,
                     p_scr_idx,
                     d_P_var,
                     IntVector<NDIM>(1),
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_p_init_fcn);

    // Everything else only gets a scratch context, which is deallocated at the end of each time step.
    // Note the forces need ghost cells for modifying the RHS to account for non-homogenous boundary conditions.
    int thn_idx, f_p_idx, f_un_idx, f_us_idx;
    registerVariable(thn_idx, d_thn_cc_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(f_p_idx, d_f_cc_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(f_un_idx, d_f_un_sc_var, IntVector<NDIM>(1), getScratchContext());
    registerVariable(f_us_idx, d_f_us_sc_var, IntVector<NDIM>(1), getScratchContext());

    // Drawing variables. Velocities are node centered that get allocated with the current context.
    // Pressure is cell centered, so we can just use the current context for that.
    // Note: Using the current context for the velocities means drawing variables are allocated with the state
    // variables. This is important because the interface has no way of separately allocating and deallocating "drawing"
    // data.
    if (d_visit_writer)
    {
        d_un_draw_var = new NodeVariable<NDIM, double>(d_object_name + "::Un_draw", NDIM);
        d_us_draw_var = new NodeVariable<NDIM, double>(d_object_name + "::Us_draw", NDIM);
        int un_draw_idx, us_draw_idx;
        registerVariable(un_draw_idx, d_un_draw_var, IntVector<NDIM>(0), getCurrentContext());
        registerVariable(us_draw_idx, d_us_draw_var, IntVector<NDIM>(0), getCurrentContext());

        d_visit_writer->registerPlotQuantity("Un", "VECTOR", un_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Un_" + std::to_string(d), "SCALAR", un_draw_idx, d, 1.0, "NODE");

        d_visit_writer->registerPlotQuantity("Us", "VECTOR", us_draw_idx, 0, 1.0, "NODE");
        for (int d = 0; d < NDIM; ++d)
            d_visit_writer->registerPlotQuantity("Us_" + std::to_string(d), "SCALAR", us_draw_idx, d, 1.0, "NODE");

        d_visit_writer->registerPlotQuantity("P", "SCALAR", p_cur_idx, 0, 1.0, "NODE");
    }

    // Create the hierarchy data operations
    auto hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops =
        hier_ops_manager->getOperationsDouble(new CellVariable<NDIM, double>("cc_var"), hierarchy, true /*get_unique*/);
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), hierarchy, true /*get_unique*/);

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
    // Do any kind of Hierarchy specific initialization on a given patch level.
    // Note: All initialization and regridding of state variables is managed by HierarchyIntegrator.
    // Example uses of this include mantaining a vorticity value to see where things should be refined or checking norms
    // of quantities before and after regridding. intentionally blank.
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM>> hierarchy,
    const int level_num,
    const double error_data_time,
    const int tag_idx,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    // Fill in the tag_idx with 1 if the cell index on a given level should be refined.
    // tag_idx corresponds to patch data CellData<NDIM, int>.
    // TODO: Implement reasonable refinement criteria.
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
    int coarsest_level,
    int finest_level)
{
    // Do any kind of Hierarchy specific configuration on a new patch hierarchy.
    // If there are synchronization objects that need to be created, do that here.

    // Reset the hierarchy operations objects.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(coarsest_level, finest_level);

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(coarsest_level, finest_level);
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                        const double new_time,
                                                                        const int num_cycles)
{
    // Do anything that needs to be done before we call integrateHierarchy().
    INSHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Pull out current solution components
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_cur_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getCurrentContext());
    const int us_cur_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getCurrentContext());
    const int p_cur_idx = var_db->mapVariableAndContextToIndex(d_P_var, getCurrentContext());
    const int un_new_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getNewContext());
    const int us_new_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getNewContext());
    const int p_new_idx = var_db->mapVariableAndContextToIndex(d_P_var, getNewContext());
    const int un_scr_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getScratchContext());
    const int us_scr_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getScratchContext());
    const int p_scr_idx = var_db->mapVariableAndContextToIndex(d_P_var, getScratchContext());

    // Allocate scratch and new data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Set initial guess
    d_hier_sc_data_ops->copyData(un_new_idx, un_cur_idx);
    d_hier_sc_data_ops->copyData(us_new_idx, us_cur_idx);
    d_hier_cc_data_ops->copyData(p_new_idx, p_cur_idx);
    d_hier_sc_data_ops->copyData(un_scr_idx, un_cur_idx);
    d_hier_sc_data_ops->copyData(us_scr_idx, us_cur_idx);
    d_hier_cc_data_ops->copyData(p_scr_idx, p_cur_idx);

    // Set up our solution vector
    d_sol_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::sol_vec", d_hierarchy, coarsest_ln, finest_ln);

    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    d_sol_vec->addComponent(d_un_sc_var, un_scr_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_sol_vec->addComponent(d_us_sc_var, us_scr_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_sol_vec->addComponent(d_P_var, p_scr_idx, wgt_cc_idx, d_hier_cc_data_ops);

    // Set up the RHS vector
    const int f_un_idx = var_db->mapVariableAndContextToIndex(d_f_un_sc_var, getScratchContext());
    const int f_us_idx = var_db->mapVariableAndContextToIndex(d_f_us_sc_var, getScratchContext());
    const int f_p_idx = var_db->mapVariableAndContextToIndex(d_f_cc_var, getScratchContext());
    d_rhs_vec = new SAMRAIVectorReal<NDIM, double>(d_object_name + "::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    d_rhs_vec->addComponent(d_f_un_sc_var, f_un_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_rhs_vec->addComponent(d_f_us_sc_var, f_us_idx, wgt_sc_idx, d_hier_sc_data_ops);
    d_rhs_vec->addComponent(d_f_cc_var, f_p_idx, wgt_cc_idx, d_hier_cc_data_ops);

    // Grab the theta index
    const int thn_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getScratchContext());
    // Set the correct data if applicable
    if (d_thn_fcn)
        d_thn_fcn->setDataOnPatchHierarchy(
            thn_idx, d_thn_cc_var, d_hierarchy, current_time, false, coarsest_ln, finest_ln);

    // Set up null vectors (if applicable)
    d_nul_vecs.resize(1);
    d_nul_vecs[0] = d_sol_vec->cloneVector("PressureNull"); // should delete the vector at the end
    d_nul_vecs[0]->allocateVectorData();
    d_nul_vecs[0]->setToScalar(0.0);
    // Pull out pressure component and set to constant
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> p_data = d_nul_vecs[0]->getComponentPatchData(2, *patch);
            p_data->fillAll(1.0);
        }
    }

    // Set the forcing data if applicable.
    d_rhs_vec->setToScalar(0.0);
    if (d_f_un_fcn) d_f_un_fcn->setDataOnPatchHierarchy(f_un_idx, d_f_un_sc_var, d_hierarchy, current_time);
    if (d_f_us_fcn) d_f_us_fcn->setDataOnPatchHierarchy(f_us_idx, d_f_us_sc_var, d_hierarchy, current_time);
    if (d_f_p_fcn) d_f_p_fcn->setDataOnPatchHierarchy(f_p_idx, d_f_cc_var, d_hierarchy, current_time);

    // Need to fill in ghost cells for theta to do interpolation
    {
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comp(1);
        ghost_cell_comp[0] = ITC(thn_idx,
                                 "CONSERVATIVE_LINEAR_REFINE",
                                 false,
                                 "NONE",
                                 "LINEAR",
                                 true,
                                 nullptr); // defaults to fill corner
        HierarchyGhostCellInterpolation ghost_cell_fill;
        ghost_cell_fill.initializeOperatorState(ghost_cell_comp, d_hierarchy, 0, d_hierarchy->getFinestLevelNumber());
        ghost_cell_fill.fillData(0.0);
    }

    // set-up RHS to treat viscosity and drag with backward Euler or Implicit Trapezoidal Rule:
    // RHS = f(n) + C*theta_i(n)*u_i(n) + D1*(pressure + viscous + drag) for  i = n, s
    double D1 = std::numeric_limits<double>::signaling_NaN();
    double D2 = std::numeric_limits<double>::signaling_NaN();
    const double C = d_rho / dt;

    switch (d_viscous_time_stepping_type)
    {
    case TRAPEZOIDAL_RULE:
        D1 = 0.5;
        D2 = -0.5;
        break;

    case BACKWARD_EULER:
        D1 = 0.0;
        D2 = -1.0;
        break;

    default:
        TBOX_ERROR("Unknown time stepping type " + enum_to_string<TimeSteppingType>(d_viscous_time_stepping_type) +
                   ". Valid options are BACKWARD_EULER and TRAPEZOIDAL_RULE.");
    }

    VCTwoFluidStaggeredStokesOperator RHS_op("RHS_op", true);
    RHS_op.setCandDCoefficients(C, D1);
    RHS_op.setDragCoefficient(1.0, 1.0, 1.0);
    RHS_op.setViscosityCoefficient(1.0, 1.0);
    RHS_op.setThnIdx(thn_idx);

    // Store results of applying stokes operator in f2_vec
    Pointer<SAMRAIVectorReal<NDIM, double>> f2_vec = d_rhs_vec->cloneVector("f2_vec");
    f2_vec->allocateVectorData();
    RHS_op.initializeOperatorState(*d_sol_vec, *f2_vec);
    RHS_op.apply(*d_sol_vec, *f2_vec);

    // Sum f_vec and f2_vec and store result in f_vec
    d_rhs_vec->add(d_rhs_vec, f2_vec, true);
    f2_vec->deallocateVectorData();

    // Set up the operators and solvers needed to solve the linear system.
    d_stokes_op = new VCTwoFluidStaggeredStokesOperator("stokes_op", true);
    d_stokes_op->setCandDCoefficients(C, D2);
    d_stokes_op->setDragCoefficient(1.0, 1.0, 1.0);
    d_stokes_op->setViscosityCoefficient(1.0, 1.0);
    d_stokes_op->setThnIdx(thn_idx);

    d_stokes_solver = new PETScKrylovLinearSolver("solver", d_solver_db, "solver_");
    d_stokes_solver->setOperator(d_stokes_op);

    // Now create a preconditioner
    if (d_use_preconditioner)
    {
        d_precond_op = new VCTwoFluidStaggeredStokesBoxRelaxationFACOperator(
            "KrylovPrecondStrategy", "Krylov_precond_", d_w, C, D2);
        d_precond_op->setThnIdx(thn_idx);
        d_stokes_precond = new FullFACPreconditioner("KrylovPrecond", d_precond_op, d_precond_db, "Krylov_precond_");
        d_stokes_precond->setNullspace(false, d_nul_vecs);
        d_stokes_solver->setPreconditioner(d_stokes_precond);
    }

    d_stokes_solver->setNullspace(false, d_nul_vecs);
    d_stokes_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec);

    // Set thn_cc_idx on the dense hierarchy.
    if (d_use_preconditioner)
    {
        Pointer<PatchHierarchy<NDIM>> dense_hierarchy = d_stokes_precond->getDenseHierarchy();
        // Allocate data
        for (int ln = 0; ln <= dense_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = dense_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(thn_idx)) level->allocatePatchData(thn_idx, 0.0);
        }
        d_thn_fcn->setDataOnPatchHierarchy(
            thn_idx, d_thn_cc_var, dense_hierarchy, 0.0, false, 0, dense_hierarchy->getFinestLevelNumber());
        {
            // Also fill in theta ghost cells
            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comp(1);
            ghost_cell_comp[0] = ITC(thn_idx,
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
        }
    }

    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
INSVCTwoFluidStaggeredHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                              const double new_time,
                                                              const int cycle_num)
{
    INSHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_new_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getNewContext());
    const int us_new_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getNewContext());
    const int p_new_idx = var_db->mapVariableAndContextToIndex(d_P_var, getNewContext());

    const double dt = new_time - current_time;

    // Set the initial guess for the system to be the most recent approximation to t^{n+1}
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(0), un_new_idx);
    d_hier_sc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(1), un_new_idx);
    d_hier_cc_data_ops->copyData(d_sol_vec->getComponentDescriptorIndex(2), p_new_idx);

    // Solve for un(n+1), us(n+1), p(n+1).
    d_stokes_solver->solveSystem(*d_sol_vec, *d_rhs_vec);

    // Reset the solve vector to copy the "scratch" data into the "new" data
    d_hier_sc_data_ops->copyData(un_new_idx, d_sol_vec->getComponentDescriptorIndex(0));
    d_hier_sc_data_ops->copyData(us_new_idx, d_sol_vec->getComponentDescriptorIndex(1));
    d_hier_cc_data_ops->copyData(p_new_idx, d_sol_vec->getComponentDescriptorIndex(2));
    return;
} // integrateHierarchy

void
INSVCTwoFluidStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                         const double new_time,
                                                                         const bool skip_synchronize_new_state_data,
                                                                         const int num_cycles)
{
    // Do anything that needs to be done after integrateHierarchy().
    INSHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Deallocate patch data on the dense hierarchy
    if (d_use_preconditioner)
    {
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int thn_cc_idx = var_db->mapVariableAndContextToIndex(d_thn_cc_var, getScratchContext());
        Pointer<PatchHierarchy<NDIM>> dense_hierarchy = d_stokes_precond->getDenseHierarchy();
        for (int ln = 0; ln <= dense_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = dense_hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(thn_cc_idx)) level->deallocatePatchData(thn_cc_idx);
        }
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
INSVCTwoFluidStaggeredHierarchyIntegrator::regridProjection()
{
    // During regridding, the coarsening and refining operations can introduce errors. Here, we project the velocity
    // onto a divergence free field.
    // TODO: Do we need to implement this? How do we project something that has a co-incompressibility condition.
    // TODO: Check divergence norms to ensure our regridding process does not introduce too large of errors.
}

double
INSVCTwoFluidStaggeredHierarchyIntegrator::getStableTimestep(Pointer<Patch<NDIM>> /*patch*/) const
{
    return d_dt_init;
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::synchronizeHierarchyDataSpecialized(VariableContextType ctx_type)
{
    // intentionally blank.
}

void
INSVCTwoFluidStaggeredHierarchyIntegrator::setupPlotDataSpecialized()
{
    // Interpolate velocities to node centered.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int un_draw_idx = var_db->mapVariableAndContextToIndex(d_un_draw_var, getCurrentContext());
    const int us_draw_idx = var_db->mapVariableAndContextToIndex(d_us_draw_var, getCurrentContext());
    const int un_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getCurrentContext());
    const int us_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getCurrentContext());
    const int un_scr_idx = var_db->mapVariableAndContextToIndex(d_un_sc_var, getScratchContext());
    const int us_scr_idx = var_db->mapVariableAndContextToIndex(d_us_sc_var, getScratchContext());

    static const bool synch_cf_interface = true;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate scratch data if necessary
    allocatePatchData(un_scr_idx, 0.0, coarsest_ln, finest_ln);
    allocatePatchData(us_scr_idx, 0.0, coarsest_ln, finest_ln);

    // We need ghost cells to interpolate to nodes.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_comp = {
        ITC(un_scr_idx, un_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr),
        ITC(us_scr_idx, us_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr)
    };
    HierarchyGhostCellInterpolation hier_bdry_fill;
    hier_bdry_fill.initializeOperatorState(ghost_comp, d_hierarchy);
    hier_bdry_fill.fillData(0.0);

    // interpolate
    HierarchyMathOps hier_math_ops("hier_math_ops", d_hierarchy);
    hier_math_ops.setPatchHierarchy(d_hierarchy);
    hier_math_ops.resetLevels(coarsest_ln, finest_ln);
    hier_math_ops.interp(un_draw_idx,
                         d_un_draw_var,
                         synch_cf_interface,
                         un_scr_idx,
                         d_un_sc_var,
                         nullptr,
                         d_integrator_time,
                         synch_cf_interface);
    hier_math_ops.interp(us_draw_idx,
                         d_us_draw_var,
                         synch_cf_interface,
                         us_scr_idx,
                         d_us_sc_var,
                         nullptr,
                         d_integrator_time,
                         synch_cf_interface);

    // Deallocate scratch data
    deallocatePatchData(un_scr_idx, coarsest_ln, finest_ln);
    deallocatePatchData(us_scr_idx, coarsest_ln, finest_ln);
}

int
INSVCTwoFluidStaggeredHierarchyIntegrator::getNumberOfCycles() const
{
    return 1;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
