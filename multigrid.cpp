// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/StaggeredStokesSolverManager.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/app_namespaces.h>

#include "ibtk/CartCellDoubleQuadraticRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LinearOperator.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

// Local includes
#include "VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.h"
#include "VCTwoFluidStaggeredStokesOperator.h"
#include "tests/multigrid/FullFACPreconditioner.h"

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "stokes.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        grid_geometry->addSpatialRefineOperator(new CartCellDoubleQuadraticRefine()); // refine op for cell-centered
                                                                                      // variables
        grid_geometry->addSpatialRefineOperator(new CartSideDoubleRT0Refine()); // refine op for side-centered variables
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // State variables: Velocity and pressure.
        Pointer<SideVariable<NDIM, double>> un_sc_var = new SideVariable<NDIM, double>("un_sc");
        Pointer<SideVariable<NDIM, double>> us_sc_var = new SideVariable<NDIM, double>("us_sc");
        Pointer<CellVariable<NDIM, double>> p_cc_var = new CellVariable<NDIM, double>("p_cc");

        // variable coefficient: theta_n
        Pointer<CellVariable<NDIM, double>> thn_cc_var = new CellVariable<NDIM, double>("thn_cc");

        // Results of operator "forces" and "divergence"
        Pointer<SideVariable<NDIM, double>> f_un_sc_var = new SideVariable<NDIM, double>("f_un_sc");
        Pointer<SideVariable<NDIM, double>> f_us_sc_var = new SideVariable<NDIM, double>("f_us_sc");
        Pointer<CellVariable<NDIM, double>> f_cc_var = new CellVariable<NDIM, double>("f_cc");

        // Error terms.
        Pointer<SideVariable<NDIM, double>> e_un_sc_var = new SideVariable<NDIM, double>("e_un_sc");
        Pointer<SideVariable<NDIM, double>> e_us_sc_var = new SideVariable<NDIM, double>("e_us_sc");
        Pointer<CellVariable<NDIM, double>> e_cc_var = new CellVariable<NDIM, double>("e_cc");

        // Register patch data indices...
        const int un_sc_idx = var_db->registerVariableAndContext(un_sc_var, ctx, IntVector<NDIM>(1));
        const int us_sc_idx = var_db->registerVariableAndContext(us_sc_var, ctx, IntVector<NDIM>(1));
        const int p_cc_idx = var_db->registerVariableAndContext(p_cc_var, ctx, IntVector<NDIM>(1));
        const int thn_cc_idx =
            var_db->registerVariableAndContext(thn_cc_var, ctx, IntVector<NDIM>(1)); // 1 layer of ghost cells
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(1));
        const int f_un_sc_idx = var_db->registerVariableAndContext(f_un_sc_var, ctx, IntVector<NDIM>(1));
        const int f_us_sc_idx = var_db->registerVariableAndContext(f_us_sc_var, ctx, IntVector<NDIM>(1));
        const int e_un_sc_idx = var_db->registerVariableAndContext(e_un_sc_var, ctx, IntVector<NDIM>(1));
        const int e_us_sc_idx = var_db->registerVariableAndContext(e_us_sc_var, ctx, IntVector<NDIM>(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));

        // Drawing variables
        Pointer<CellVariable<NDIM, double>> draw_un_var = new CellVariable<NDIM, double>("draw_un", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_fn_var = new CellVariable<NDIM, double>("draw_fn", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_en_var = new CellVariable<NDIM, double>("draw_en", NDIM);
        const int draw_un_idx = var_db->registerVariableAndContext(draw_un_var, ctx);
        const int draw_fn_idx = var_db->registerVariableAndContext(draw_fn_var, ctx);
        const int draw_en_idx = var_db->registerVariableAndContext(draw_en_var, ctx);

        Pointer<CellVariable<NDIM, double>> draw_us_var = new CellVariable<NDIM, double>("draw_us", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_fs_var = new CellVariable<NDIM, double>("draw_fs", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_es_var = new CellVariable<NDIM, double>("draw_es", NDIM);
        const int draw_us_idx = var_db->registerVariableAndContext(draw_us_var, ctx);
        const int draw_fs_idx = var_db->registerVariableAndContext(draw_fs_var, ctx);
        const int draw_es_idx = var_db->registerVariableAndContext(draw_es_var, ctx);

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity("Pressure", "SCALAR", p_cc_idx);
        visit_data_writer->registerPlotQuantity("Thn", "SCALAR", thn_cc_idx);
        visit_data_writer->registerPlotQuantity("RHS_P", "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity("error_p", "SCALAR", e_cc_idx);

        visit_data_writer->registerPlotQuantity("Un", "VECTOR", draw_un_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("Un_" + std::to_string(d), "SCALAR", draw_un_idx, d);
        }

        visit_data_writer->registerPlotQuantity("RHS_Un", "VECTOR", draw_fn_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("RHS_Un_" + std::to_string(d), "SCALAR", draw_fn_idx, d);
        }

        visit_data_writer->registerPlotQuantity("error_un", "VECTOR", draw_en_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("error_un_" + std::to_string(d), "SCALAR", draw_en_idx, d);
        }

        visit_data_writer->registerPlotQuantity("Us", "VECTOR", draw_us_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("Us_" + std::to_string(d), "SCALAR", draw_us_idx, d);
        }

        visit_data_writer->registerPlotQuantity("RHS_Us", "VECTOR", draw_fs_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("RHS_Us_" + std::to_string(d), "SCALAR", draw_fs_idx, d);
        }

        visit_data_writer->registerPlotQuantity("error_us", "VECTOR", draw_es_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("error_us_" + std::to_string(d), "SCALAR", draw_es_idx, d);
        }

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(un_sc_idx, 0.0);
            level->allocatePatchData(us_sc_idx, 0.0);
            level->allocatePatchData(f_un_sc_idx, 0.0);
            level->allocatePatchData(f_us_sc_idx, 0.0);
            level->allocatePatchData(e_un_sc_idx, 0.0);
            level->allocatePatchData(e_us_sc_idx, 0.0);
            level->allocatePatchData(p_cc_idx, 0.0);
            level->allocatePatchData(thn_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(draw_un_idx, 0.0);
            level->allocatePatchData(draw_fn_idx, 0.0);
            level->allocatePatchData(draw_en_idx, 0.0);
            level->allocatePatchData(draw_us_idx, 0.0);
            level->allocatePatchData(draw_fs_idx, 0.0);
            level->allocatePatchData(draw_es_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(un_sc_var, un_sc_idx, h_sc_idx);
        u_vec.addComponent(us_sc_var, us_sc_idx, h_sc_idx);
        u_vec.addComponent(p_cc_var, p_cc_idx, h_cc_idx);
        f_vec.addComponent(f_un_sc_var, f_un_sc_idx, h_sc_idx);
        f_vec.addComponent(f_us_sc_var, f_us_sc_idx, h_sc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        e_vec.addComponent(e_un_sc_var, e_un_sc_idx, h_sc_idx);
        e_vec.addComponent(e_us_sc_var, e_us_sc_idx, h_sc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(0.0);

        // There is one pressure and two velocity nullspaces (one for each component)
        std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> null_vecs(3);
        null_vecs[0] = u_vec.cloneVector("PressureNull");
        null_vecs[1] = u_vec.cloneVector("VelNull0");
        null_vecs[2] = u_vec.cloneVector("VelNull1");
        null_vecs[0]->allocateVectorData();
        null_vecs[0]->setToScalar(0.0);
        null_vecs[1]->allocateVectorData();
        null_vecs[1]->setToScalar(0.0);
        null_vecs[2]->allocateVectorData();
        null_vecs[2]->setToScalar(0.0);
        // Pull out pressure component and set to constant
        {
            for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double>> p_data = null_vecs[0]->getComponentPatchData(2, *patch);
                    p_data->fillAll(1.0);

                    Pointer<SideData<NDIM, double>> un_data = null_vecs[1]->getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double>> us_data = null_vecs[1]->getComponentPatchData(1, *patch);
                    for (SideIterator<NDIM> si(patch->getBox(), 0); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        (*un_data)(idx) = 1.0;
                        (*us_data)(idx) = 1.0;
                    }

                    un_data = null_vecs[2]->getComponentPatchData(0, *patch);
                    us_data = null_vecs[2]->getComponentPatchData(1, *patch);
                    for (SideIterator<NDIM> si(patch->getBox(), 1); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        (*un_data)(idx) = 1.0;
                        (*us_data)(idx) = 1.0;
                    }
                }
            }
        }

        // Setup velocity and pressures functions.
        muParserCartGridFunction un_fcn("un", app_initializer->getComponentDatabase("un"), grid_geometry);
        muParserCartGridFunction us_fcn("us", app_initializer->getComponentDatabase("us"), grid_geometry);
        muParserCartGridFunction p_fcn("p", app_initializer->getComponentDatabase("p"), grid_geometry);

        // Set up Thn functions
        muParserCartGridFunction thn_fcn("thn", app_initializer->getComponentDatabase("thn"), grid_geometry);

        // Setup exact solution functions
        muParserCartGridFunction f_un_fcn("f_un", app_initializer->getComponentDatabase("f_un"), grid_geometry);
        muParserCartGridFunction f_us_fcn("f_us", app_initializer->getComponentDatabase("f_us"), grid_geometry);
        muParserCartGridFunction f_p_fcn("f_p", app_initializer->getComponentDatabase("f_p"), grid_geometry);

        f_un_fcn.setDataOnPatchHierarchy(f_un_sc_idx, f_un_sc_var, patch_hierarchy, 0.0);
        f_us_fcn.setDataOnPatchHierarchy(f_us_sc_idx, f_us_sc_var, patch_hierarchy, 0.0);
        f_p_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);
        thn_fcn.setDataOnPatchHierarchy(thn_cc_idx, thn_cc_var, patch_hierarchy, 0.0);

        un_fcn.setDataOnPatchHierarchy(e_un_sc_idx, e_un_sc_var, patch_hierarchy, 0.0);
        us_fcn.setDataOnPatchHierarchy(e_us_sc_idx, e_us_sc_var, patch_hierarchy, 0.0);
        p_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);

        // Setup the stokes operator
        Pointer<VCTwoFluidStaggeredStokesOperator> stokes_op = new VCTwoFluidStaggeredStokesOperator("stokes_op", true);

        stokes_op->setThnIdx(thn_cc_idx);

        Pointer<PETScKrylovLinearSolver> krylov_solver =
            new PETScKrylovLinearSolver("solver", app_initializer->getComponentDatabase("KrylovSolver"), "solver_");
        krylov_solver->setOperator(stokes_op);

        // Now create a preconditioner
        Pointer<VCTwoFluidStaggeredStokesBoxRelaxationFACOperator> fac_precondition_strategy =
            new VCTwoFluidStaggeredStokesBoxRelaxationFACOperator(
                "KrylovPrecondStrategy",
                // app_initializer->getComponentDatabase("KrylovPrecondStrategy"),
                "Krylov_precond_", input_db->getDouble("w"));
        fac_precondition_strategy->setThnIdx(thn_cc_idx);
        Pointer<FullFACPreconditioner> Krylov_precond =
            new FullFACPreconditioner("KrylovPrecond",
                                      fac_precondition_strategy,
                                      app_initializer->getComponentDatabase("KrylovPrecond"),
                                      input_db->getInteger("multigrid_max_levels"),
                                      "Krylov_precond_");
        bool use_precond = input_db->getBool("USE_PRECOND");
        Krylov_precond->setNullspace(false, null_vecs);
        if (use_precond) krylov_solver->setPreconditioner(Krylov_precond);

        krylov_solver->setNullspace(false, null_vecs);
        krylov_solver->initializeSolverState(u_vec, f_vec);
        // We need to set thn_cc_idx on the dense hierarchy.
        // TODO: find a better way to do this
        if (use_precond)
        {
            Pointer<PatchHierarchy<NDIM>> dense_hierarchy = Krylov_precond->getDenseHierarchy();
            // Allocate data
            for (int ln = 0; ln <= dense_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = dense_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(thn_cc_idx)) level->allocatePatchData(thn_cc_idx, 0.0);
            }
            thn_fcn.setDataOnPatchHierarchy(
                thn_cc_idx, thn_cc_var, dense_hierarchy, 0.0, false, 0, dense_hierarchy->getFinestLevelNumber());
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
        }

        krylov_solver->solveSystem(u_vec, f_vec);

        // Deallocate data
        if (use_precond)
        {
            Pointer<PatchHierarchy<NDIM>> dense_hierarchy = Krylov_precond->getDenseHierarchy();
            for (int ln = 0; ln <= dense_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = dense_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(thn_cc_idx)) level->deallocatePatchData(thn_cc_idx);
            }
        }

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false),  // numerical
                       Pointer<SAMRAIVectorReal<NDIM, double>>(&e_vec, false)); // analytical
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(0, patch_hierarchy->getFinestLevelNumber());
        // just computes quadrature weights
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        pout << "Error in u_n :\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_un_sc_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_un_sc_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(e_un_sc_idx, wgt_sc_idx) << "\n";

        pout << "Error in u_s :\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_us_sc_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_us_sc_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(e_us_sc_idx, wgt_sc_idx) << "\n";

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        pout << "Error in p :\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(e_cc_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(e_cc_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(e_cc_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(draw_un_idx, draw_un_var, un_sc_idx, un_sc_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_fn_idx, draw_fn_var, f_un_sc_idx, f_un_sc_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_en_idx, draw_en_var, e_un_sc_idx, e_un_sc_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_us_idx, draw_us_var, us_sc_idx, us_sc_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_fs_idx, draw_fs_var, f_us_sc_idx, f_us_sc_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_es_idx, draw_es_var, e_us_sc_idx, e_us_sc_var, nullptr, 0.0, synch_cf_interface);

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
        // Deallocate level data
        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(un_sc_idx);
            level->deallocatePatchData(us_sc_idx);
            level->deallocatePatchData(f_un_sc_idx);
            level->deallocatePatchData(f_us_sc_idx);
            level->deallocatePatchData(e_un_sc_idx);
            level->deallocatePatchData(e_us_sc_idx);
            level->deallocatePatchData(p_cc_idx);
            level->deallocatePatchData(thn_cc_idx);
            level->deallocatePatchData(f_cc_idx);
            level->deallocatePatchData(e_cc_idx);
            level->deallocatePatchData(draw_un_idx);
            level->deallocatePatchData(draw_fn_idx);
            level->deallocatePatchData(draw_en_idx);
            level->deallocatePatchData(draw_us_idx);
            level->deallocatePatchData(draw_fs_idx);
            level->deallocatePatchData(draw_es_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
