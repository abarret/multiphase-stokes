#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/CFINSForcing.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBLagrangianForceStrategySet.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/IBTargetPointForceSpec.h>
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
#include "multiphase/IBMultiphaseCrossLinks.h"
#include "multiphase/IBMultiphaseHierarchyIntegrator.h"
#include "multiphase/INSVCTwoFluidStaggeredHierarchyIntegrator.h"
#include "multiphase/ScaleStress.h"
#include "multiphase/StressRelaxation.h"
using namespace multiphase;

int finest_ln;
std::array<int, NDIM> N;
double alpha = 0.0;
double g = 0.0;
double upper_perim = 0.0;
double lower_perim = 0.0;
double MFAC = 0.0;
double dx = 0.0;
double L = 0.0;
double K = 0.0;
std::vector<int> num_nodes_per_struct;
std::vector<double> ds;

VectorNd
upper_channel(const double s, const double t)
{
    VectorNd x;
    x(0) = s;
    x(1) = alpha / (2.0 * M_PI) * (1.0 + g * std::sin(2.0 * M_PI * (s - t)));
    return x;
}

VectorNd
lower_channel(const double s, const double t)
{
    VectorNd x;
    x(0) = s;
    x(1) = -alpha / (2.0 * M_PI) * (1.0 + g * std::sin(2.0 * M_PI * (s - t)));
    return x;
}
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
    if (strct_num == 0)
    {
        // Generating upper level of channel.
        // Determine lag grid spacing.
        ds[strct_num] = MFAC * upper_perim * dx;
        num_vertices = std::floor(L / ds[strct_num]);
        vertex_posn.resize(num_vertices);
        for (int i = 0; i < num_vertices; ++i)
        {
            VectorNd x = upper_channel(ds[strct_num] * (static_cast<double>(i) + 0.5), 0.0);
            vertex_posn[i] = x;
        }
    }
    else if (strct_num == 1)
    {
        // Generating lower level of channel.
        ds[strct_num] = MFAC * lower_perim * dx;
        num_vertices = std::floor(L / ds[strct_num]);
        vertex_posn.resize(num_vertices);
        for (int i = 0; i < num_vertices; ++i)
        {
            VectorNd x = lower_channel(ds[strct_num] * (static_cast<double>(i) + 0.5), 0.0);
            vertex_posn[i] = x;
        }
    }
    num_nodes_per_struct[strct_num] = num_vertices;
    return;
}

void
generate_tethers(const unsigned int& strct_num,
                 const int& ln,
                 std::multimap<int, IBRedundantInitializer::TargetSpec>& tg_pt_spec,
                 void* /*ctx*/)
{
    if (ln != finest_ln) return;
    for (int k = 0; k < num_nodes_per_struct[strct_num]; ++k)
    {
        IBRedundantInitializer::TargetSpec e;
        e.stiffness = K / ds[strct_num];
        e.damping = 0.0;
        tg_pt_spec.insert(std::make_pair(k, e));
    }
}

void
move_tethers(LDataManager* data_manager, const double time)
{
    const std::pair<int, int>& upper_lag_idxs = data_manager->getLagrangianStructureIndexRange(0, finest_ln);

    // Update both local and ghost nodes.
    Pointer<LMesh> mesh = data_manager->getLMesh(finest_ln);
    std::vector<LNode*> nodes = mesh->getLocalNodes();
    const std::vector<LNode*>& ghost_nodes = mesh->getGhostNodes();
    nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());

    // Also need reference information.
    Pointer<LData> X_ref_data = data_manager->getLData(data_manager->INIT_POSN_DATA_NAME, finest_ln);
    double* X_ref_vals = X_ref_data->getVecArray()->data();

    for (const auto node : nodes)
    {
        const auto& force_spec = node->getNodeDataItem<IBTargetPointForceSpec>();
        if (!force_spec) continue;
        const int lag_idx = node->getLagrangianIndex();
        const int petsc_idx = node->getLocalPETScIndex();
        IBTK::Point& X_target = force_spec->getTargetPointPosition();
        // Detect with side of channel we are on.
        if (upper_lag_idxs.first <= lag_idx && lag_idx < upper_lag_idxs.second)
            X_target = upper_channel(X_ref_vals[petsc_idx * NDIM], time);
        else
            X_target = lower_channel(X_ref_vals[petsc_idx * NDIM], time);
    }

    X_ref_data->restoreArrays();
}

void
print_max_dist_from_target(LDataManager* data_manager, const double time)
{
    double max_disp = 0.0;
    Pointer<LMesh> mesh = data_manager->getLMesh(finest_ln);
    const std::pair<int, int>& upper_lag_idxs = data_manager->getLagrangianStructureIndexRange(0, finest_ln);
    std::vector<LNode*> nodes = mesh->getLocalNodes();

    Pointer<LData> X_data = data_manager->getLData(data_manager->POSN_DATA_NAME, finest_ln);
    double* X_vals = X_data->getVecArray()->data();
    Pointer<LData> X_ref_data = data_manager->getLData(data_manager->INIT_POSN_DATA_NAME, finest_ln);
    double* X_ref_vals = X_ref_data->getVecArray()->data();

    for (const auto& node : nodes)
    {
        const int petsc_idx = node->getLocalPETScIndex();
        const int lag_idx = node->getLagrangianIndex();
        Eigen::Map<VectorNd> x(&X_vals[petsc_idx * NDIM]), x_ref(&X_ref_vals[petsc_idx * NDIM]);

        VectorNd x_target;
        if (upper_lag_idxs.first <= lag_idx && lag_idx < upper_lag_idxs.second)
            x_target = upper_channel(x_ref(0), time);
        else
            x_target = lower_channel(x_ref(0), time);
        max_disp = std::max(max_disp, (x - x_target).norm());
    }

    pout << "Max distance from target point: " << max_disp << "\n";
    pout << "Fraction of grid spacing:       " << max_disp / dx << "\n";
    if ((max_disp / dx) > 0.1)
    {
        // Error out if the distance from target point is more than 10% of a grid cell.
        TBOX_ERROR("Too far away from target point!");
    }
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
        std::vector<std::string> struct_list_network = { "network_upper", "newtork_lower" };
        std::vector<std::string> struct_list_solvent = { "solvent_upper", "solvent_lower" };
        ds.resize(2);
        num_nodes_per_struct.resize(2);
        N[0] = N[1] = input_db->getInteger("N");
        finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
        alpha = input_db->getDouble("ALPHA");
        K = input_db->getDouble("K_TETHER");
        g = input_db->getDouble("GAMMA");
        L = input_db->getDouble("L");
        MFAC = input_db->getDouble("MFAC");
        dx = L / N[0];
        upper_perim = 1.02423522856; // alpha * (L*M_PI + g * std::sin(L*M_PI)*std::sin(L*M_PI)) / (2.0 * M_PI * M_PI);
        lower_perim = upper_perim;

        ibs_initializer->setStructureNamesOnLevel(finest_ln, struct_list_solvent);
        ibs_initializer->registerInitStructureFunction(generate_structure);
        ibs_initializer->registerInitTargetPtFunction(generate_tethers);
        ib_solvent_ops->registerLInitStrategy(ibs_initializer);
        Pointer<IBStandardForceGen> ibs_spring_forces = new IBStandardForceGen();
        ib_solvent_ops->registerIBLagrangianForceFunction(ibs_spring_forces);

        Pointer<IBMultiphaseCrossLinks> ib_cross_forces = new IBMultiphaseCrossLinks(
            ib_network_ops, ib_solvent_ops, patch_hierarchy, input_db->getDouble("KAPPA") / 0.0160037, 0.0);
        time_integrator->registerCrossLinkStrategy(ib_cross_forces);

        ibn_initializer->setStructureNamesOnLevel(finest_ln, struct_list_network);
        ibn_initializer->registerInitStructureFunction(generate_structure);
        ib_network_ops->registerLInitStrategy(ibn_initializer);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> solvent_data_writer = new LSiloDataWriter(
            "SolventDataWriter", input_db->getDatabase("Main")->getString("viz_dump_dirname") + "/solvent", false);
        Pointer<LSiloDataWriter> network_data_writer = new LSiloDataWriter(
            "NetworkDataWriter", input_db->getDatabase("Main")->getString("viz_dump_dirname") + "/network", false);
        if (uses_visit)
        {
            ibn_initializer->registerLSiloDataWriter(network_data_writer);
            ib_network_ops->registerLSiloDataWriter(network_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ibs_initializer->registerLSiloDataWriter(solvent_data_writer);
            ib_solvent_ops->registerLSiloDataWriter(solvent_data_writer);
            adv_diff_integrator->registerVisItDataWriter(visit_data_writer);
        }

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

        Pointer<CFINSForcing> cf_forcing;
        Pointer<ScaleStress> stress_scale =
            new ScaleStress("ScaleStress", ins_integrator->getNetworkVolumeFractionVariable(), adv_diff_integrator);
        if (input_db->getBool("USE_CF"))
        {
            Pointer<INSHierarchyIntegrator> ins_cf_integrator = ins_integrator;
            cf_forcing = new CFINSForcing("CFINSForcing",
                                          app_initializer->getComponentDatabase("CFINSForcing"),
                                          ins_cf_integrator,
                                          grid_geometry,
                                          adv_diff_integrator,
                                          visit_data_writer);
            cf_forcing->setSigmaScaleFcn(stress_scale);
            Pointer<StressRelaxation> stress_relax =
                new StressRelaxation("StressRelax",
                                     app_initializer->getComponentDatabase("CFINSForcing"),
                                     ins_integrator->getNetworkVariable(),
                                     ins_integrator,
                                     ins_integrator->getNetworkVolumeFractionVariable(),
                                     adv_diff_integrator);
            cf_forcing->registerRelaxationOperator(stress_relax);
            ins_integrator->setForcingFunctions(cf_forcing, nullptr);
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
            solvent_data_writer->writePlotData(iteration_num, loop_time);
            network_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();
            print_max_dist_from_target(ib_network_ops->getLDataManager(), loop_time);
            print_max_dist_from_target(ib_solvent_ops->getLDataManager(), loop_time);
            move_tethers(ib_network_ops->getLDataManager(), loop_time);
            move_tethers(ib_solvent_ops->getLDataManager(), loop_time);

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
                solvent_data_writer->writePlotData(iteration_num, loop_time);
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
