// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBLagrangianForceStrategySet.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <array>

// Local includes
#include "IBMultiphaseCrossLinks.h"
#include "IBMultiphaseHierarchyIntegrator.h"
#include "INSVCTwoFluidStaggeredHierarchyIntegrator.h"

int finest_ln;
std::array<int, NDIM> N;
SAMRAI::tbox::Array<int> num_node;
SAMRAI::tbox::Array<double> ds;
int num_node_circum, num_node_radial;
double dr;
void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* /*ctx*/)
{
    if (ln != finest_ln)
    {
        num_vertices = 0;
        vertex_posn.resize(num_vertices);
    }
    double beta = 0.35;
    double alpha = 0.25 * 0.25 / beta;
    double A = M_PI * alpha * beta; // Area of ellipse
    double R = sqrt(A / M_PI);      // Radius of disc with equivalent area as ellipse
    double perim = 2 * M_PI * R;    // Perimenter of equivalent disc
    double dx = 1.0 / static_cast<double>(N[0]);
    //         double sdf = perim/(1.0/(3.0*64.0))/4.0;
    //         double aa = ceil(sdf);
    //         double aaa = aa*4.0;
    //         double aaaa = dx/64.0*aaa;
    //         int aaaaa = static_cast<int>(aaaa);
    num_node[strct_num] = static_cast<int>((1.0 / (dx * 64.0)) * ceil(perim / (1.0 / (3.0 * 64.0)) / 4.0) * 4.0);
    ds[strct_num] = 2.0 * M_PI * R / num_node[strct_num];
    num_vertices = num_node[strct_num];
    vertex_posn.resize(num_vertices);
    for (std::vector<IBTK::Point>::iterator it = vertex_posn.begin(); it != vertex_posn.end(); ++it)
    {
        Point& X = *it;
        int num = std::distance(vertex_posn.begin(), it);
        double theta = 2.0 * M_PI * num / num_vertices;
        X(0) = 0.5 + alpha * std::cos(theta);
        X(1) = 0.5 + beta * std::sin(theta);
    }
    return;
}

void
generate_springs(
    const unsigned int& strct_num,
    const int& ln,
    std::multimap<int, IBRedundantInitializer::Edge>& spring_map,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>&
        spring_spec,
    void* /*ctx*/)
{
    if (ln != finest_ln) return;
    double K = 1.0;
    for (int k = 0; k < num_node[strct_num]; ++k)
    {
        IBRedundantInitializer::Edge e;
        std::vector<double> parameters(2);
        int force_fcn_idx = 0;
        e.first = k;
        e.second = (k + 1) % num_node[strct_num];
        if (e.first > e.second) std::swap(e.first, e.second);
        spring_map.insert(std::make_pair(e.first, e));
        IBRedundantInitializer::SpringSpec spec_data;
        parameters[0] = K / ds[strct_num]; // spring constant
        parameters[1] = 0.0;               // resting length
        spec_data.parameters = parameters;
        spec_data.force_fcn_idx = force_fcn_idx;
        spring_spec.insert(std::make_pair(e, spec_data));
    }
    return;
}

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    IBTKInit init(argc, argv);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create dummy boundary condition objects
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (int d = 0; d < NDIM; ++d)
            u_bc_coefs[d] =
                new muParserRobinBcCoefs("U_dummy", app_initializer->getComponentDatabase("U_dummy"), grid_geometry);
        RobinBcCoefStrategy<NDIM>* p_bc_coef =
            new muParserRobinBcCoefs("P_dummy", app_initializer->getComponentDatabase("P_dummy"), grid_geometry);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSVCTwoFluidStaggeredHierarchyIntegrator> ins_integrator =
            new INSVCTwoFluidStaggeredHierarchyIntegrator(
                "INSIntegrator", app_initializer->getComponentDatabase("INSIntegrator"), u_bc_coefs, p_bc_coef, false);
        Pointer<IBMethod> ib_solvent_ops =
            new IBMethod("IBSolventMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<IBMethod> ib_network_ops =
            new IBMethod("IBNetworkMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<IBMultiphaseHierarchyIntegrator> time_integrator =
            new IBMultiphaseHierarchyIntegrator("IBHierarchyIntegrator",
                                                app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                ib_solvent_ops,
                                                ib_network_ops,
                                                ins_integrator);
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator("AdvDiffIntegrator",
                                                       app_initializer->getComponentDatabase("AdvDiffIntegrator"),
                                                       true /*register_for_restart*/);

        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
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

        // Configure the IB solver.
        Pointer<IBRedundantInitializer> ibn_initializer = new IBRedundantInitializer(
            "IBNRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        Pointer<IBRedundantInitializer> ibs_initializer = new IBRedundantInitializer(
            "IBSRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        std::vector<std::string> struct_list_network = { "network" };
        std::vector<std::string> struct_list_solvent = { "solvent" };
        ds.resizeArray(1);
        num_node.resizeArray(1);
        N[0] = N[1] = input_db->getInteger("N");
        finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
        double beta = 0.35;
        double alpha = 0.25 * 0.25 / beta;
        double A = M_PI * alpha * beta; // Area of ellipse
        double R = sqrt(A / M_PI);      // Radius of disc with equivalent area as ellipse
        ibn_initializer->setStructureNamesOnLevel(finest_ln, struct_list_network);
        ibn_initializer->registerInitStructureFunction(generate_structure);
        ibn_initializer->registerInitSpringDataFunction(generate_springs);
        ib_network_ops->registerLInitStrategy(ibn_initializer);
        Pointer<IBStandardForceGen> ibn_spring_forces = new IBStandardForceGen();
        Pointer<IBMultiphaseCrossLinks> ibn_cross_forces = new IBMultiphaseCrossLinks(
            ib_solvent_ops->getLDataManager(), input_db->getDouble("KAPPA") / (2.0 * M_PI * R / num_node[0]));
        std::vector<Pointer<IBLagrangianForceStrategy>> net_forces = { ibn_spring_forces, ibn_cross_forces };
        Pointer<IBLagrangianForceStrategySet> ibn_force_fcn =
            new IBLagrangianForceStrategySet(net_forces.begin(), net_forces.end());
        ib_network_ops->registerIBLagrangianForceFunction(ibn_force_fcn);

        ibs_initializer->setStructureNamesOnLevel(finest_ln, struct_list_solvent);
        ibs_initializer->registerInitStructureFunction(generate_structure);
        ibs_initializer->registerInitSpringDataFunction(generate_springs);
        ib_solvent_ops->registerLInitStrategy(ibs_initializer);
        Pointer<IBStandardForceGen> ibs_spring_forces = new IBStandardForceGen();
        Pointer<IBMultiphaseCrossLinks> ibs_cross_forces = new IBMultiphaseCrossLinks(
            ib_network_ops->getLDataManager(), input_db->getDouble("KAPPA") / (2.0 * M_PI * R / num_node[0]));
        std::vector<Pointer<IBLagrangianForceStrategy>> sol_forces = { ibs_spring_forces, ibs_cross_forces };
        Pointer<IBLagrangianForceStrategySet> ibs_force_fcn =
            new IBLagrangianForceStrategySet(sol_forces.begin(), sol_forces.end());
        ib_solvent_ops->registerIBLagrangianForceFunction(ibs_force_fcn);

        // Set up coefficients
        const double eta_n = input_db->getDouble("ETA_N");
        const double eta_s = input_db->getDouble("ETA_S");
        const double xi = input_db->getDouble("XI");
        const double nu_n = input_db->getDouble("NU_N");
        const double nu_s = input_db->getDouble("NU_S");
        ins_integrator->setViscosityCoefficient(eta_n, eta_s);
        ins_integrator->setDragCoefficient(xi, nu_n, nu_s);

        // Set up volume fractions
        Pointer<CartGridFunction> thn_init_fcn =
            new muParserCartGridFunction("thn_init", app_initializer->getComponentDatabase("thn"), grid_geometry);
        ins_integrator->setInitialNetworkVolumeFraction(thn_init_fcn);
        ins_integrator->advectNetworkVolumeFraction(adv_diff_integrator);

        // Boundary conditions are periodic.
        // TODO: Remove this assumption
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, nullptr);
        ins_integrator->registerPhysicalBoundaryConditions(bc_coefs);

        // Create Eulerian initial condition specification objects.  These
        // objects also are used to specify exact solution values for error
        // analysis.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        ins_integrator->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        ins_integrator->registerPressureInitialConditions(p_init);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        Pointer<LSiloDataWriter> network_data_writer = new LSiloDataWriter("NetworkDataWriter", "viz_n_IB2d", false);
        if (uses_visit)
        {
            ibn_initializer->registerLSiloDataWriter(network_data_writer);
            ib_network_ops->registerLSiloDataWriter(network_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ibs_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_solvent_ops->registerLSiloDataWriter(silo_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_solvent_ops->freeLInitStrategy();
        ib_network_ops->freeLInitStrategy();
        ibs_initializer.setNull();
        ibn_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
            network_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
                network_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        for (int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete p_bc_coef;
    } // cleanup dynamically allocated objects prior to shutdown
} // main
