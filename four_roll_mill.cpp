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

#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
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
#include "FRMThnForcing.h"
#include "INSVCTwoFluidStaggeredHierarchyIntegrator.h"

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
        Pointer<INSVCTwoFluidStaggeredHierarchyIntegrator> ins_integrator =
            new INSVCTwoFluidStaggeredHierarchyIntegrator(
                "FluidSolver",
                app_initializer->getComponentDatabase("INSVCTwoFluidStaggeredHierarchyIntegrator"),
                grid_geometry,
                true /*register_for_restart*/);

        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator("AdvDiffIntegrator",
                                                       app_initializer->getComponentDatabase("AdvDiffIntegrator"),
                                                       true /*register_for_restart*/);
        grid_geometry->addSpatialRefineOperator(new CartCellDoubleQuadraticRefine()); // refine op for cell-centered
                                                                                      // variables
        grid_geometry->addSpatialRefineOperator(new CartSideDoubleRT0Refine()); // refine op for side-centered variables
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ins_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        const double eta_n = input_db->getDouble("ETA_N");
        const double eta_s = input_db->getDouble("ETA_S");
        const double xi = input_db->getDouble("XI");
        const double nu_n = input_db->getDouble("NU_N");
        const double nu_s = input_db->getDouble("NU_S");
        ins_integrator->setViscosityCoefficient(eta_n, eta_s);
        ins_integrator->setDragCoefficient(xi, nu_n, nu_s);

        // Setup velocity and pressures functions.
        Pointer<CartGridFunction> un_init =
            new muParserCartGridFunction("un", app_initializer->getComponentDatabase("un"), grid_geometry);
        Pointer<CartGridFunction> us_init =
            new muParserCartGridFunction("us", app_initializer->getComponentDatabase("us"), grid_geometry);
        Pointer<CartGridFunction> p_init =
            new muParserCartGridFunction("p", app_initializer->getComponentDatabase("p"), grid_geometry);
        ins_integrator->setInitialData(un_init, us_init, p_init);

        // Set up Thn functions
        Pointer<CartGridFunction> thn_init_fcn =
            new muParserCartGridFunction("thn_init", app_initializer->getComponentDatabase("thn"), grid_geometry);
        ins_integrator->setInitialNetworkVolumeFraction(thn_init_fcn);
        ins_integrator->advectNetworkVolumeFraction(adv_diff_integrator);

        // Set up forcing function
        Pointer<CartGridFunction> fn_fcn =
            new FRMThnForcing(ins_integrator->getNetworkVolumeFractionVariable(), adv_diff_integrator, true);
        Pointer<CartGridFunction> fs_fcn =
            new FRMThnForcing(ins_integrator->getNetworkVolumeFractionVariable(), adv_diff_integrator, false);
        ins_integrator->setForcingFunctions(fn_fcn, fs_fcn, nullptr);

        // Set up visualizations
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        ins_integrator->registerVisItDataWriter(visit_data_writer);

        // Initialize the INS integrator
        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Get some time stepping information.
        unsigned int iteration_num = ins_integrator->getIntegratorStep();
        double loop_time = ins_integrator->getIntegratorTime();
        double time_end = ins_integrator->getEndTime();
        double dt = 0.0;

        // Visualization files info.
        double viz_dump_time_interval = input_db->getDouble("VIZ_DUMP_TIME_INTERVAL");
        double next_viz_dump_time = 0.0;
        // At specified intervals, write visualization files
        if (IBTK::abs_equal_eps(loop_time, next_viz_dump_time, 0.1 * dt) || loop_time >= next_viz_dump_time)
        {
            pout << "\nWriting visualization files...\n\n";
            ins_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            next_viz_dump_time += viz_dump_time_interval;
        }
        // Main time step loop
        while (!IBTK::rel_equal_eps(loop_time, time_end) && ins_integrator->stepsRemaining())
        {
            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = ins_integrator->getMaximumTimeStepSize();
            ins_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            iteration_num += 1;
            // At specified intervals, write visualization files
            if (IBTK::abs_equal_eps(loop_time, next_viz_dump_time, 0.1 * dt) || loop_time >= next_viz_dump_time)
            {
                pout << "\nWriting visualization files...\n\n";
                ins_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                next_viz_dump_time += viz_dump_time_interval;
            }
        }

        // Compute errors.
        Pointer<CartGridFunction> un_exact_fcn =
            new muParserCartGridFunction("UN_exact", app_initializer->getComponentDatabase("un_exact"), grid_geometry);
        Pointer<CartGridFunction> us_exact_fcn =
            new muParserCartGridFunction("US_exact", app_initializer->getComponentDatabase("us_exact"), grid_geometry);
        Pointer<CartGridFunction> p_exact_fcn =
            new muParserCartGridFunction("P_exact", app_initializer->getComponentDatabase("p_exact"), grid_geometry);

        // Get the current data.
        Pointer<SideVariable<NDIM, double>> un_var = ins_integrator->getNetworkVariable();
        Pointer<SideVariable<NDIM, double>> us_var = ins_integrator->getSolventVariable();
        Pointer<CellVariable<NDIM, double>> p_var = ins_integrator->getPressureVariable();

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = ins_integrator->getCurrentContext();
        const int un_idx = var_db->mapVariableAndContextToIndex(un_var, ctx);
        const int us_idx = var_db->mapVariableAndContextToIndex(us_var, ctx);
        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, ctx);

        // Create scratch data.
        Pointer<VariableContext> exa_ctx = var_db->getContext("Exact");
        const int un_exa_idx = var_db->registerVariableAndContext(un_var, exa_ctx);
        const int us_exa_idx = var_db->registerVariableAndContext(us_var, exa_ctx);
        const int p_exa_idx = var_db->registerVariableAndContext(p_var, exa_ctx);

        // Allocate scratch data
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(un_exa_idx, loop_time);
            level->allocatePatchData(us_exa_idx, loop_time);
            level->allocatePatchData(p_exa_idx, loop_time);
        }

        // Fill in exact values
        un_exact_fcn->setDataOnPatchHierarchy(
            un_exa_idx, un_var, patch_hierarchy, loop_time, false, coarsest_ln, finest_ln);
        us_exact_fcn->setDataOnPatchHierarchy(
            us_exa_idx, us_var, patch_hierarchy, loop_time, false, coarsest_ln, finest_ln);
        p_exact_fcn->setDataOnPatchHierarchy(
            p_exa_idx, p_var, patch_hierarchy, loop_time, false, coarsest_ln, finest_ln);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy, coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        // Compute error and print out norms.
        hier_sc_data_ops.subtract(un_idx, un_idx, un_exa_idx);
        hier_sc_data_ops.subtract(us_idx, us_idx, us_exa_idx);
        hier_cc_data_ops.subtract(p_idx, p_idx, p_exa_idx);
        pout << "Printing error norms\n\n";
        pout << "Newtork velocity\n";
        pout << "Un L1-norm:  " << hier_sc_data_ops.L1Norm(un_idx, wgt_sc_idx) << "\n";
        pout << "Un L1-norm:  " << hier_sc_data_ops.L2Norm(un_idx, wgt_sc_idx) << "\n";
        pout << "Un max-norm: " << hier_sc_data_ops.maxNorm(un_idx, wgt_sc_idx) << "\n\n";

        pout << "Solvent velocity\n";
        pout << "Us L1-norm:  " << hier_sc_data_ops.L1Norm(us_idx, wgt_sc_idx) << "\n";
        pout << "Us L1-norm:  " << hier_sc_data_ops.L2Norm(us_idx, wgt_sc_idx) << "\n";
        pout << "Us max-norm: " << hier_sc_data_ops.maxNorm(us_idx, wgt_sc_idx) << "\n\n";

        pout << "Pressure\n";
        pout << "P L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n";
        pout << "P L1-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n";
        pout << "P max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n";

        // Print extra viz files for the error
        pout << "\nWriting visualization files...\n\n";
        ins_integrator->setupPlotData();
        visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(un_exa_idx);
            level->deallocatePatchData(us_exa_idx);
            level->deallocatePatchData(p_exa_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main