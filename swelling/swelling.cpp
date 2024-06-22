#include "multiphase/CFMultiphaseOldroydB.h"
#include "multiphase/MultiphaseStandardHierarchyIntegrator.h"

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

using namespace multiphase;

class OsmoticPressure : public CartGridFunction
{
public:
    OsmoticPressure(std::string object_name, Pointer<MultiphaseStandardHierarchyIntegrator> integrator)
        : CartGridFunction(std::move(object_name)), d_integrator(integrator)
    {
    }

    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatch(int data_idx,
                        Pointer<hier::Variable<NDIM>> var,
                        Pointer<Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        Pointer<PatchLevel<NDIM>> level = Pointer<PatchLevel<NDIM>>(NULL)) override
    {
        Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(data_idx);
        Pointer<CellData<NDIM, double>> thn_data =
            patch->getPatchData(d_integrator->getNetworkVolumeFractionVariable(), d_integrator->getCurrentContext());

        F_data->fillAll(0.0);
        if (initial_time) return;
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
            {
                const SideIndex<NDIM>& idx = si();
                const CellIndex<NDIM>& idx_up = idx.toCell(1);
                const CellIndex<NDIM>& idx_low = idx.toCell(0);
                double thn_up = (*thn_data)(idx_up);
                double thn_low = (*thn_data)(idx_low);
                double val_up = thn_up * psi(thn_up);
                double val_low = thn_low * psi(thn_low);
                (*F_data)(idx) = -(val_up - val_low) / dx[axis];
            }
        }
    }

private:
    Pointer<MultiphaseStandardHierarchyIntegrator> d_integrator;

    double psi(const double thn)
    {
        return thn - 0.1;
    }
};

struct ExactParams
{
    ExactParams(Pointer<Database> input_db)
    {
        thn_init = input_db->getDouble("thn_init");
        R0 = input_db->getDouble("R0");
        a = input_db->getDouble("a");
        psi0 = input_db->getDouble("psi0");
        mu = input_db->getDouble("mu");
        lambda = input_db->getDouble("lambda");
    }
    double thn_init = std::numeric_limits<double>::quiet_NaN();
    double R0 = std::numeric_limits<double>::quiet_NaN();
    double a = std::numeric_limits<double>::quiet_NaN();
    double psi0 = std::numeric_limits<double>::quiet_NaN();
    double mu = std::numeric_limits<double>::quiet_NaN();
    double lambda = std::numeric_limits<double>::quiet_NaN();
};

class ThnExactFcn : public CartGridFunction
{
public:
    ThnExactFcn(std::string object_name, const ExactParams& params)
        : CartGridFunction(std::move(object_name)), d_params(params)
    {
    }

    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatch(int data_idx,
                        Pointer<hier::Variable<NDIM>> var,
                        Pointer<Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        Pointer<PatchLevel<NDIM>> level = Pointer<PatchLevel<NDIM>>(NULL)) override
    {
        Pointer<CellData<NDIM, double>> thn_data = patch->getPatchData(data_idx);
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        const double* const xlow = pgeom->getXLower();
        const hier::Index<NDIM>& idx_low = patch->getBox().lower();

        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            VectorNd x;
            for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
            const double R_idx = x.norm();

            const double R_val = R(data_time);

            if (R_val > R_idx)
            {
                (*thn_data)(idx) = d_params.thn_init * d_params.R0 * d_params.R0 / (R_val * R_val);
            }
            else
            {
                (*thn_data)(idx) = 0.0;
            }
        }
    }

private:
    const ExactParams& d_params;

    double R(const double t)
    {
        double R0_sq = d_params.R0 * d_params.R0;
        double NT = d_params.thn_init * R0_sq;
        double Req_sq = NT / d_params.a;

        double R_sq =
            Req_sq + (R0_sq - Req_sq) * std::exp(-d_params.psi0 * d_params.a / (d_params.mu + d_params.lambda) * t);

        return std::sqrt(R_sq);
    }
};

class VelExactFcn : public CartGridFunction
{
public:
    VelExactFcn(std::string object_name, const ExactParams& params)
        : CartGridFunction(std::move(object_name)), d_params(params)
    {
    }

    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatch(int data_idx,
                        Pointer<hier::Variable<NDIM>> var,
                        Pointer<Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        Pointer<PatchLevel<NDIM>> level = Pointer<PatchLevel<NDIM>>(NULL)) override
    {
        Pointer<SideData<NDIM, double>> u_sc_data = patch->getPatchData(data_idx);
        Pointer<NodeData<NDIM, double>> u_nc_data = patch->getPatchData(data_idx);
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        const double* const xlow = pgeom->getXLower();
        const hier::Index<NDIM>& idx_low = patch->getBox().lower();

        if (u_sc_data)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SideIndex<NDIM>& idx = si();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                    const double th_idx = std::atan2(x[1], x[0]);
                    const double u_r = u(x, data_time);

                    if (axis == 0)
                        (*u_sc_data)(idx) = u_r * std::cos(th_idx);
                    else
                        (*u_sc_data)(idx) = u_r * std::sin(th_idx);
                }
            }
        }
        else if (u_nc_data)
        {
            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = ni();
                VectorNd x;
                for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
                const double th_idx = std::atan2(x[1], x[0]);
                const double u_r = u(x, data_time);
                (*u_nc_data)(idx, 0) = u_r * std::cos(th_idx);
                (*u_nc_data)(idx, 1) = u_r * std::sin(th_idx);
            }
        }
    }

private:
    const ExactParams& d_params;
    double u(const VectorNd& x, const double t)
    {
        const double R_idx = x.norm();
        const double R_val = R(t);

        return Rp(t) * R_idx / R_val;
    }

    double R(const double t)
    {
        double R0_sq = d_params.R0 * d_params.R0;
        double NT = d_params.thn_init * R0_sq;
        double Req_sq = NT / d_params.a;

        double R_sq =
            Req_sq + (R0_sq - Req_sq) * std::exp(-d_params.psi0 * d_params.a / (d_params.mu + d_params.lambda) * t);

        return std::sqrt(R_sq);
    }

    double Rp(const double t)
    {
        double R0_sq = d_params.R0 * d_params.R0;
        double NT = d_params.thn_init * R0_sq;
        double Req_sq = NT / d_params.a;
        double R_inv = 1.0 / R(t);

        double Rp_val = (R0_sq - Req_sq) * (-d_params.psi0 * d_params.a / (d_params.mu + d_params.lambda)) *
                        std::exp(-d_params.psi0 * d_params.a / (d_params.mu + d_params.lambda) * t);

        return 0.5 * Rp_val * R_inv;
    }
};

void compute_errors(Pointer<ThnExactFcn> thn_fcn,
                    Pointer<VelExactFcn> un_fcn,
                    const int thn_exact_idx,
                    Pointer<CellVariable<NDIM, double>> thn_exact_var,
                    const int thn_error_idx,
                    Pointer<CellVariable<NDIM, double>> thn_error_var,
                    const int un_exact_idx,
                    Pointer<SideVariable<NDIM, double>> un_exact_var,
                    const int un_error_idx,
                    Pointer<SideVariable<NDIM, double>> un_error_var,
                    const int un_draw_idx,
                    Pointer<NodeVariable<NDIM, double>> un_draw_var,
                    Pointer<PatchHierarchy<NDIM>> hierarchy,
                    Pointer<MultiphaseStandardHierarchyIntegrator> integrator,
                    const double t);

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
        Pointer<MultiphaseStandardHierarchyIntegrator> ins_integrator = new MultiphaseStandardHierarchyIntegrator(
            "FluidSolver",
            app_initializer->getComponentDatabase("INSVCTwoFluidStaggeredHierarchyIntegrator"),
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
        ins_integrator->setViscosityCoefficient(eta_n, eta_s, 0.0, 0.0);
        ins_integrator->setDragCoefficient(xi, nu_n, nu_s);

        // Create boundary condition objects for velocity
        std::vector<RobinBcCoefStrategy<NDIM>*> un_bc_coefs(NDIM, nullptr), us_bc_coefs(NDIM, nullptr);
        for (int d = 0; d < NDIM; ++d)
        {
            std::string un_bc_coef_name = "bc_coefs_" + std::to_string(d);
            un_bc_coefs[d] =
                new muParserRobinBcCoefs(un_bc_coef_name, input_db->getDatabase(un_bc_coef_name), grid_geometry);
            std::string us_bc_coef_name = "bc_coefs_" + std::to_string(d);
            us_bc_coefs[d] =
                new muParserRobinBcCoefs(us_bc_coef_name, input_db->getDatabase(us_bc_coef_name), grid_geometry);
        }
        auto thn_bc_coef =
            std::make_unique<muParserRobinBcCoefs>("thn_bc", input_db->getDatabase("thn_bc_coefs"), grid_geometry);
        ins_integrator->registerPhysicalBoundaryConditions(un_bc_coefs, us_bc_coefs);

        ExactParams params(input_db->getDatabase("ExactParams"));
        Pointer<VelExactFcn> un_exact_fcn = new VelExactFcn("VelExactFcn", params);
        Pointer<ThnExactFcn> thn_exact_fcn = new ThnExactFcn("ThnExactFcn", params);

        // Set up visualizations
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();

        // Set up exact variables
        Pointer<CellVariable<NDIM, double>> thn_exact_var = new CellVariable<NDIM, double>("thn_exact");
        Pointer<CellVariable<NDIM, double>> thn_error_var = new CellVariable<NDIM, double>("thn_error");
        Pointer<SideVariable<NDIM, double>> un_exact_var = new SideVariable<NDIM, double>("un_exact");
        Pointer<SideVariable<NDIM, double>> un_error_var = new SideVariable<NDIM, double>("us_error");
        Pointer<NodeVariable<NDIM, double>> un_draw_var = new NodeVariable<NDIM, double>("un_draw", NDIM);

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int thn_exact_idx = var_db->registerVariableAndContext(thn_exact_var, var_db->getContext("CTX"));
        const int thn_error_idx = var_db->registerVariableAndContext(thn_error_var, var_db->getContext("CTX"));
        const int un_exact_idx = var_db->registerVariableAndContext(un_exact_var, var_db->getContext("CTX"));
        const int un_error_idx = var_db->registerVariableAndContext(un_error_var, var_db->getContext("CTX"));
        const int un_draw_idx = var_db->registerVariableAndContext(un_draw_var, var_db->getContext("CTX"));
        std::set<int> idxs = { thn_exact_idx, thn_error_idx, un_exact_idx, un_error_idx, un_draw_idx };
        visit_data_writer->registerPlotQuantity("thn_exact", "SCALAR", thn_exact_idx);
        visit_data_writer->registerPlotQuantity("thn_error", "SCALAR", thn_error_idx);
        visit_data_writer->registerPlotQuantity("un_exact", "VECTOR", un_draw_idx);

        // Setup velocity and pressures functions.
        Pointer<CartGridFunction> un_init =
            new muParserCartGridFunction("un", app_initializer->getComponentDatabase("un"), grid_geometry);
        Pointer<CartGridFunction> us_init =
            new muParserCartGridFunction("us", app_initializer->getComponentDatabase("us"), grid_geometry);
        Pointer<CartGridFunction> p_init =
            new muParserCartGridFunction("p", app_initializer->getComponentDatabase("p"), grid_geometry);
        ins_integrator->setInitialData(un_init, us_init, p_init);

        // Set up Thn functions
        ins_integrator->setInitialNetworkVolumeFraction(thn_exact_fcn);
        ins_integrator->registerVolumeFractionBoundaryConditions(thn_bc_coef.get());
        ins_integrator->advectNetworkVolumeFraction(adv_diff_integrator);

        Pointer<CartGridFunction> op_fcn = new OsmoticPressure("OPFcn", ins_integrator);
        ins_integrator->setForcingFunctions(op_fcn, nullptr);

        // Initialize the INS integrator
        ins_integrator->registerVisItDataWriter(visit_data_writer);
        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        input_db->printClassData(plog);

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
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            allocate_patch_data(idxs, patch_hierarchy, loop_time, coarsest_ln, finest_ln);
            compute_errors(thn_exact_fcn,
                           un_exact_fcn,
                           thn_exact_idx,
                           thn_exact_var,
                           thn_error_idx,
                           thn_error_var,
                           un_exact_idx,
                           un_exact_var,
                           un_error_idx,
                           un_error_var,
                           un_draw_idx,
                           un_draw_var,
                           patch_hierarchy,
                           ins_integrator,
                           loop_time);
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            deallocate_patch_data(idxs, patch_hierarchy, coarsest_ln, finest_ln);
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
                allocate_patch_data(idxs, patch_hierarchy, loop_time, coarsest_ln, finest_ln);
                compute_errors(thn_exact_fcn,
                               un_exact_fcn,
                               thn_exact_idx,
                               thn_exact_var,
                               thn_error_idx,
                               thn_error_var,
                               un_exact_idx,
                               un_exact_var,
                               un_error_idx,
                               un_error_var,
                               un_draw_idx,
                               un_draw_var,
                               patch_hierarchy,
                               ins_integrator,
                               loop_time);
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                deallocate_patch_data(idxs, patch_hierarchy, coarsest_ln, finest_ln);
                next_viz_dump_time += viz_dump_time_interval;
            }

            // At specified intervals, write restart files.
            const bool last_step = !ins_integrator->stepsRemaining();
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
compute_errors(Pointer<ThnExactFcn> thn_fcn,
               Pointer<VelExactFcn> un_fcn,
               const int thn_exact_idx,
               Pointer<CellVariable<NDIM, double>> thn_exact_var,
               const int thn_error_idx,
               Pointer<CellVariable<NDIM, double>> thn_error_var,
               const int un_exact_idx,
               Pointer<SideVariable<NDIM, double>> un_exact_var,
               const int un_error_idx,
               Pointer<SideVariable<NDIM, double>> un_error_var,
               const int un_draw_idx,
               Pointer<NodeVariable<NDIM, double>> un_draw_var,
               Pointer<PatchHierarchy<NDIM>> hierarchy,
               Pointer<MultiphaseStandardHierarchyIntegrator> integrator,
               const double t)
{
    thn_fcn->setDataOnPatchHierarchy(thn_exact_idx, thn_exact_var, hierarchy, t);
    un_fcn->setDataOnPatchHierarchy(un_exact_idx, un_exact_var, hierarchy, t);
    un_fcn->setDataOnPatchHierarchy(un_draw_idx, un_draw_var, hierarchy, t);
    // Compute errors
    HierarchyMathOps hier_math_ops("hier_math_ops", hierarchy);
    hier_math_ops.setPatchHierarchy(hierarchy);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy);
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(hierarchy);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int thn_idx = var_db->mapVariableAndContextToIndex(integrator->getNetworkVolumeFractionVariable(),
                                                             integrator->getCurrentContext());
    const int un_idx =
        var_db->mapVariableAndContextToIndex(integrator->getNetworkVariable(), integrator->getCurrentContext());
    hier_cc_data_ops.subtract(thn_error_idx, thn_idx, thn_exact_idx);
    hier_sc_data_ops.subtract(un_error_idx, un_idx, un_exact_idx);
    // Print error norms
    pout << "Computing error norms at time " << t << "\n";
    pout << "Errors in volume fraction:\n";
    pout << "  L1-norm:  " << hier_cc_data_ops.L1Norm(thn_error_idx, wgt_cc_idx) << "\n";
    pout << "  L2-norm:  " << hier_cc_data_ops.L2Norm(thn_error_idx, wgt_cc_idx) << "\n";
    pout << "  max-norm: " << hier_cc_data_ops.maxNorm(thn_error_idx, wgt_cc_idx) << "\n";
    pout << "Errors in network velocity:\n";
    pout << "  L1-norm:  " << hier_sc_data_ops.L1Norm(un_error_idx, wgt_sc_idx) << "\n";
    pout << "  L2-norm:  " << hier_sc_data_ops.L2Norm(un_error_idx, wgt_sc_idx) << "\n";
    pout << "  max-norm: " << hier_sc_data_ops.maxNorm(un_error_idx, wgt_sc_idx) << "\n";
}
