#include "multiphase/MultiphaseStandardHierarchyIntegrator.h"

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

using namespace multiphase;

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
    SAMRAIManager::setMaxNumberPatchDataEntries(8192);

    PetscOptionsSetValue(nullptr, "-solver_ksp_rtol", "1.0e-12");

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "stokes.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<MultiphaseStandardHierarchyIntegrator> ins_integrator = new MultiphaseStandardHierarchyIntegrator(
            "FluidSolver", app_initializer->getComponentDatabase("INSVCTwoFluidStaggeredHierarchyIntegrator"), false);
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

        bool using_var_xi = input_db->getBool("USING_VAR_XI");
        const double eta_n = input_db->getDouble("ETA_N");
        const double eta_s = input_db->getDouble("ETA_S");
        ins_integrator->setViscosityCoefficient(eta_n, eta_s);
        if (using_var_xi)
        {
            Pointer<CartGridFunction> xi_fcn =
                new muParserCartGridFunction("xi", app_initializer->getComponentDatabase("xi"), grid_geometry);
            ins_integrator->setDragCoefficientFunction(xi_fcn);
        }
        else
        {
            const double xi = input_db->getDouble("XI");
            const double nu = input_db->getDouble("NU");
            ins_integrator->setDragCoefficient(xi, nu, nu);
        }

        // Setup boundary conditions if necessary
        std::vector<RobinBcCoefStrategy<NDIM>*> un_bc_coefs(NDIM, nullptr), us_bc_coefs(NDIM, nullptr);
        std::unique_ptr<RobinBcCoefStrategy<NDIM>> p_bc_coef, thn_bc_coef;
        const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
        if (!periodic_domain)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                std::string un_bc_coef_name = "un_bcs_" + std::to_string(d);
                un_bc_coefs[d] =
                    new muParserRobinBcCoefs(un_bc_coef_name, input_db->getDatabase(un_bc_coef_name), grid_geometry);
                std::string us_bc_coef_name = "us_bcs_" + std::to_string(d);
                us_bc_coefs[d] =
                    new muParserRobinBcCoefs(us_bc_coef_name, input_db->getDatabase(us_bc_coef_name), grid_geometry);
            }

            // Create boundary condition objects for volume fraction
            thn_bc_coef =
                std::make_unique<muParserRobinBcCoefs>("thn_bc", input_db->getDatabase("thn_bc"), grid_geometry);
            p_bc_coef = std::make_unique<muParserRobinBcCoefs>("p_bc", input_db->getDatabase("p_bc"), grid_geometry);
        }

        ins_integrator->registerPhysicalBoundaryConditions(un_bc_coefs, us_bc_coefs, p_bc_coef.get());
        ins_integrator->registerVolumeFractionBoundaryConditions(thn_bc_coef.get());

        // Setup velocity and pressures functions.
        Pointer<CartGridFunction> un_fcn =
            new muParserCartGridFunction("un", app_initializer->getComponentDatabase("un"), grid_geometry);
        Pointer<CartGridFunction> us_fcn =
            new muParserCartGridFunction("us", app_initializer->getComponentDatabase("us"), grid_geometry);
        Pointer<CartGridFunction> p_fcn =
            new muParserCartGridFunction("p", app_initializer->getComponentDatabase("p"), grid_geometry);
        ins_integrator->setInitialData(un_fcn, us_fcn, p_fcn);

        // Set up Thn functions
        Pointer<CartGridFunction> thn_fcn =
            new muParserCartGridFunction("thn", app_initializer->getComponentDatabase("thn"), grid_geometry);
        ins_integrator->setNetworkVolumeFractionFunction(thn_fcn, true);

        // Set up forcing function
        Pointer<CartGridFunction> fn_fcn =
            new muParserCartGridFunction("f_un", app_initializer->getComponentDatabase("f_un"), grid_geometry);
        Pointer<CartGridFunction> fs_fcn =
            new muParserCartGridFunction("f_us", app_initializer->getComponentDatabase("f_us"), grid_geometry);
        Pointer<CartGridFunction> fp_fcn =
            new muParserCartGridFunction("f_p", app_initializer->getComponentDatabase("f_p"), grid_geometry);
        ins_integrator->setForcingFunctions(fn_fcn, fs_fcn, fp_fcn);

        // Set up visualizations
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        ins_integrator->registerVisItDataWriter(visit_data_writer);

        // Initialize the INS integrator
        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Create exact indices
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<NodeVariable<NDIM, double>> un_exact_var = new NodeVariable<NDIM, double>("un_exact", NDIM);
        Pointer<NodeVariable<NDIM, double>> us_exact_var = new NodeVariable<NDIM, double>("us_exact", NDIM);
        Pointer<CellVariable<NDIM, double>> p_exact_var = new CellVariable<NDIM, double>("p_exact");
        const int un_exact_idx = var_db->registerVariableAndContext(un_exact_var, var_db->getContext("Exact"));
        const int us_exact_idx = var_db->registerVariableAndContext(us_exact_var, var_db->getContext("Exact"));
        const int p_exact_idx = var_db->registerVariableAndContext(p_exact_var, var_db->getContext("Exact"));

        visit_data_writer->registerPlotQuantity("UN_exact", "VECTOR", un_exact_idx, 0, 1.0, "NODE");
        for (unsigned int d = 0; d < NDIM; ++d)
            visit_data_writer->registerPlotQuantity(
                "UN_exact_" + std::to_string(d), "SCALAR", un_exact_idx, d, 1.0, "NODE");
        visit_data_writer->registerPlotQuantity("US_exact", "VECTOR", us_exact_idx, 0, 1.0, "NODE");
        for (unsigned int d = 0; d < NDIM; ++d)
            visit_data_writer->registerPlotQuantity(
                "US_exact_" + std::to_string(d), "SCALAR", us_exact_idx, d, 1.0, "NODE");
        visit_data_writer->registerPlotQuantity("p_exact", "SCALAR", p_exact_idx, 0, 1.0, "CELL");

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
            ins_integrator->allocatePatchData(un_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
            ins_integrator->allocatePatchData(us_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
            ins_integrator->allocatePatchData(p_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
            un_fcn->setDataOnPatchHierarchy(un_exact_idx, un_exact_var, patch_hierarchy, loop_time);
            us_fcn->setDataOnPatchHierarchy(us_exact_idx, us_exact_var, patch_hierarchy, loop_time);
            p_fcn->setDataOnPatchHierarchy(p_exact_idx, p_exact_var, patch_hierarchy, loop_time);
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            ins_integrator->deallocatePatchData(un_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());
            ins_integrator->deallocatePatchData(us_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());
            ins_integrator->deallocatePatchData(p_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());
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
                ins_integrator->allocatePatchData(un_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
                ins_integrator->allocatePatchData(us_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
                ins_integrator->allocatePatchData(p_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
                un_fcn->setDataOnPatchHierarchy(un_exact_idx, un_exact_var, patch_hierarchy, loop_time);
                us_fcn->setDataOnPatchHierarchy(us_exact_idx, us_exact_var, patch_hierarchy, loop_time);
                p_fcn->setDataOnPatchHierarchy(p_exact_idx, p_exact_var, patch_hierarchy, loop_time);
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                ins_integrator->deallocatePatchData(un_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());
                ins_integrator->deallocatePatchData(us_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());
                ins_integrator->deallocatePatchData(p_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());

                next_viz_dump_time += viz_dump_time_interval;
            }
        }

        // Compute errors.
        // Get the current data.
        Pointer<SideVariable<NDIM, double>> un_var = ins_integrator->getNetworkVariable();
        Pointer<SideVariable<NDIM, double>> us_var = ins_integrator->getSolventVariable();
        Pointer<CellVariable<NDIM, double>> p_var = ins_integrator->getPressureVariable();

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
        un_fcn->setDataOnPatchHierarchy(un_exa_idx, un_var, patch_hierarchy, loop_time, false, coarsest_ln, finest_ln);
        us_fcn->setDataOnPatchHierarchy(us_exa_idx, us_var, patch_hierarchy, loop_time, false, coarsest_ln, finest_ln);
        p_fcn->setDataOnPatchHierarchy(p_exa_idx, p_var, patch_hierarchy, loop_time, false, coarsest_ln, finest_ln);

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
        pout << "Network velocity\n";
        pout << "Un L1-norm:  " << hier_sc_data_ops.L1Norm(un_idx, wgt_sc_idx) << "\n";
        pout << "Un L2-norm:  " << hier_sc_data_ops.L2Norm(un_idx, wgt_sc_idx) << "\n";
        pout << "Un max-norm: " << hier_sc_data_ops.maxNorm(un_idx, wgt_sc_idx) << "\n\n";

        pout << "Solvent velocity\n";
        pout << "Us L1-norm:  " << hier_sc_data_ops.L1Norm(us_idx, wgt_sc_idx) << "\n";
        pout << "Us L2-norm:  " << hier_sc_data_ops.L2Norm(us_idx, wgt_sc_idx) << "\n";
        pout << "Us max-norm: " << hier_sc_data_ops.maxNorm(us_idx, wgt_sc_idx) << "\n\n";

        pout << "Pressure\n";
        pout << "P L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n";
        pout << "P L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n";
        pout << "P max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n";

        std::ofstream out("output");
        out << "Printing error norms\n\n";
        out << "Network velocity\n";
        out << "Un L1-norm:  " << hier_sc_data_ops.L1Norm(un_idx, wgt_sc_idx) << "\n";
        out << "Un L2-norm:  " << hier_sc_data_ops.L2Norm(un_idx, wgt_sc_idx) << "\n";
        out << "Un max-norm: " << hier_sc_data_ops.maxNorm(un_idx, wgt_sc_idx) << "\n\n";

        out << "Solvent velocity\n";
        out << "Us L1-norm:  " << hier_sc_data_ops.L1Norm(us_idx, wgt_sc_idx) << "\n";
        out << "Us L2-norm:  " << hier_sc_data_ops.L2Norm(us_idx, wgt_sc_idx) << "\n";
        out << "Us max-norm: " << hier_sc_data_ops.maxNorm(us_idx, wgt_sc_idx) << "\n\n";

        out << "Pressure\n";
        out << "P L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n";
        out << "P L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n";
        out << "P max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n";

        // Print extra viz files for the error
        ins_integrator->setupPlotData();
        ins_integrator->allocatePatchData(un_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
        ins_integrator->allocatePatchData(us_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
        ins_integrator->allocatePatchData(p_exact_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
        un_fcn->setDataOnPatchHierarchy(un_exact_idx, un_exact_var, patch_hierarchy, loop_time);
        us_fcn->setDataOnPatchHierarchy(us_exact_idx, us_exact_var, patch_hierarchy, loop_time);
        p_fcn->setDataOnPatchHierarchy(p_exact_idx, p_exact_var, patch_hierarchy, loop_time);
        visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        ins_integrator->deallocatePatchData(un_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());
        ins_integrator->deallocatePatchData(us_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());
        ins_integrator->deallocatePatchData(p_exact_idx, 0, patch_hierarchy->getFinestLevelNumber());

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(un_exa_idx);
            level->deallocatePatchData(us_exa_idx);
            level->deallocatePatchData(p_exa_idx);
        }

        if (periodic_domain)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                delete un_bc_coefs[d];
                delete us_bc_coefs[d];
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
