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
#include <ibamr/CFINSForcing.h>
#include <ibamr/StaggeredStokesSolverManager.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/app_namespaces.h>

#include "ibtk/CartCellDoubleQuadraticRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
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
#include "CFMultiphaseOldroydB.h"
#include "INSVCTwoFluidStaggeredHierarchyIntegrator.h"

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSVCTwoFluidStaggeredHierarchyIntegrator> ins_integrator,
                 Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                 Pointer<CFINSForcing> conformation_tensor_handler,  
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

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

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellVariable<NDIM, double>> EE_var = new CellVariable<NDIM, double>("EE", 3);
        const int EE_idx = var_db->registerVariableAndContext(EE_var, var_db->getContext("CTX"));
        visit_data_writer->registerPlotQuantity("EE_0", "SCALAR", EE_idx, 0);
        visit_data_writer->registerPlotQuantity("EE_1", "SCALAR", EE_idx, 1);
        visit_data_writer->registerPlotQuantity("EE_2", "SCALAR", EE_idx, 2);

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

        bool use_cf = input_db->getBool("USE_CF");
        Pointer<CFINSForcing> cf_un_forcing;
        Pointer<CFStrategy> cf_strategy;
        if (use_cf)
        {
            Pointer<INSHierarchyIntegrator> ins_cf_integrator = ins_integrator;
            cf_un_forcing = new CFINSForcing("CFINSForcing",
                                             app_initializer->getComponentDatabase("CFINSForcing"),
                                             ins_cf_integrator,
                                             grid_geometry,
                                             adv_diff_integrator,
                                             visit_data_writer);
            cf_strategy = new CFMultiphaseOldroydB("CFOldroydB",
                                                   ins_integrator->getNetworkVolumeFractionVariable(),
                                                   adv_diff_integrator,
                                                   input_db->getDatabase("CFINSForcing"));
            cf_un_forcing->registerCFStrategy(cf_strategy);
            ins_integrator->setForcingFunctions(cf_un_forcing, nullptr);
        }

        // Initialize the INS integrator
        ins_integrator->registerVisItDataWriter(visit_data_writer);
        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Get some time stepping information.
        unsigned int iteration_num = ins_integrator->getIntegratorStep();
        double loop_time = ins_integrator->getIntegratorTime();
        double time_end = ins_integrator->getEndTime();
        double dt = 0.0;

        // Get simulation restart information
        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        // Get hierarchy data dump information
        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        // Visualization files info.
        double viz_dump_time_interval = input_db->getDouble("VIZ_DUMP_TIME_INTERVAL");
        double next_viz_dump_time = 0.0;

        // At specified intervals, write visualization files
        if (IBTK::abs_equal_eps(loop_time, next_viz_dump_time, 0.1 * dt) || loop_time >= next_viz_dump_time)
        {
            pout << "\nWriting visualization files...\n\n";
            ins_integrator->setupPlotData();
            int coarsest_ln = 0;
            int finest_ln = patch_hierarchy->getFinestLevelNumber();
            const int u_cur_idx = var_db->mapVariableAndContextToIndex(ins_integrator->getNetworkVariable(),
                                                                       ins_integrator->getCurrentContext());
            const int u_scr_idx = var_db->mapVariableAndContextToIndex(ins_integrator->getNetworkVariable(),
                                                                       ins_integrator->getScratchContext());
            ins_integrator->allocatePatchData(u_scr_idx, loop_time, coarsest_ln, finest_ln);
            ins_integrator->allocatePatchData(EE_idx, loop_time, coarsest_ln, finest_ln);
            HierarchyMathOps hier_math_ops("HierMathOps", patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<ITC> ghost_cell_comps = { ITC(
                u_scr_idx, u_cur_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE", "LINEAR", true, nullptr) };
            Pointer<HierarchyGhostCellInterpolation> hier_ghost_fill = new HierarchyGhostCellInterpolation();
            hier_ghost_fill->initializeOperatorState(ghost_cell_comps, patch_hierarchy, coarsest_ln, finest_ln);
            hier_math_ops.strain_rate(
                EE_idx, EE_var, u_scr_idx, ins_integrator->getNetworkVariable(), hier_ghost_fill, loop_time);
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            ins_integrator->deallocatePatchData(u_scr_idx, coarsest_ln, finest_ln);
            ins_integrator->deallocatePatchData(EE_idx, coarsest_ln, finest_ln);
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
                int coarsest_ln = 0;
                int finest_ln = patch_hierarchy->getFinestLevelNumber();
                const int u_cur_idx = var_db->mapVariableAndContextToIndex(ins_integrator->getNetworkVariable(),
                                                                           ins_integrator->getCurrentContext());
                const int u_scr_idx = var_db->mapVariableAndContextToIndex(ins_integrator->getNetworkVariable(),
                                                                           ins_integrator->getScratchContext());
                ins_integrator->allocatePatchData(u_scr_idx, loop_time, coarsest_ln, finest_ln);
                ins_integrator->allocatePatchData(EE_idx, loop_time, coarsest_ln, finest_ln);
                HierarchyMathOps hier_math_ops("HierMathOps", patch_hierarchy);
                hier_math_ops.resetLevels(coarsest_ln, finest_ln);
                using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
                std::vector<ITC> ghost_cell_comps = { ITC(
                    u_scr_idx, u_cur_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE", "LINEAR", true, nullptr) };
                Pointer<HierarchyGhostCellInterpolation> hier_ghost_fill = new HierarchyGhostCellInterpolation();
                hier_ghost_fill->initializeOperatorState(ghost_cell_comps, patch_hierarchy, coarsest_ln, finest_ln);
                hier_math_ops.strain_rate(
                    EE_idx, EE_var, u_scr_idx, ins_integrator->getNetworkVariable(), hier_ghost_fill, loop_time);
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                ins_integrator->deallocatePatchData(u_scr_idx, coarsest_ln, finest_ln);
                ins_integrator->deallocatePatchData(EE_idx, coarsest_ln, finest_ln);
                next_viz_dump_time += viz_dump_time_interval;
            }

            // At specified intervals, write restart files.
            const bool last_step = !ins_integrator->stepsRemaining();
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }

            // At specified intervals, store hierarchy data for post processing.
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                pout << "\nWriting hierarchy data files...\n\n";
                output_data(patch_hierarchy, ins_integrator, adv_diff_integrator, cf_un_forcing, iteration_num, loop_time, postproc_data_dump_dirname);
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSVCTwoFluidStaggeredHierarchyIntegrator> ins_integrator,
                 Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                 Pointer<CFINSForcing> conformation_tensor_handler,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;
    string file_name = data_dump_dirname + "/" + "hier_data_cf";
    char temp_buf[128];
    sprintf(temp_buf, ".%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getNetworkVariable(),    // Network velocity
                                                        ins_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getSolventVariable(),    // Solvent velocity
                                                        ins_integrator->getCurrentContext()));                                                    
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getPressureVariable(),   // Pressure
                                                        ins_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getNetworkVolumeFractionVariable(), // Network volume fraction 
                                                        adv_diff_integrator->getCurrentContext())); 
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(conformation_tensor_handler->getVariable(),      // Conformation tensor
                                                            adv_diff_integrator->getCurrentContext()));    
    pout << "network variable name: " << ins_integrator->getNetworkVariable()->getName() << "\n";
    pout << "solvent variable name: " << ins_integrator->getSolventVariable()->getName() << "\n";
    pout << "pressure variable name: " << ins_integrator->getPressureVariable()->getName() << "\n";
    pout << "Thn variable name: " << ins_integrator->getNetworkVolumeFractionVariable()->getName() << "\n";
    pout << "Conformation tensor variable name: " << conformation_tensor_handler->getVariable()->getName() << "\n";
    pout << "ins context: " << ins_integrator->getCurrentContext()->getName() << "\n";
    pout << "ins context: " << adv_diff_integrator->getCurrentContext()->getName() << "\n";
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data
