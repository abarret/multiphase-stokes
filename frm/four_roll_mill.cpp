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

        // Set up visualizations
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();

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
            new muParserCartGridFunction("FN_FCN", app_initializer->getComponentDatabase("FN_FCN"), grid_geometry);
        Pointer<CartGridFunction> fs_fcn =
            new muParserCartGridFunction("FS_FCN", app_initializer->getComponentDatabase("FS_FCN"), grid_geometry);
        ins_integrator->setForcingFunctionsScaled(fn_fcn, fs_fcn);

        // Initialize the INS integrator
        ins_integrator->registerVisItDataWriter(visit_data_writer);
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
    } // cleanup dynamically allocated objects prior to shutdown
} // main
