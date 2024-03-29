#include <ibamr/PETScKrylovStaggeredStokesSolver.h>
#include <ibamr/StaggeredStokesSolverManager.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
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
#include "multiphase/MultiphaseConvectiveManager.h"
#include "multiphase/utility_functions.h"
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

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "sc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
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

        // Create boundary condition objects for velocity
        std::vector<RobinBcCoefStrategy<NDIM>*> un_bc_coefs(NDIM, nullptr), us_bc_coefs(NDIM, nullptr);
        std::unique_ptr<RobinBcCoefStrategy<NDIM>> thn_bc_coef = nullptr;
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
        }

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // State variables: Velocity and pressure.
        Pointer<SideVariable<NDIM, double>> un_var = new SideVariable<NDIM, double>("un");
        Pointer<SideVariable<NDIM, double>> us_var = new SideVariable<NDIM, double>("us");
        Pointer<CellVariable<NDIM, double>> thn_var = new CellVariable<NDIM, double>("thn");

        // Results of operator
        Pointer<SideVariable<NDIM, double>> N_un_var = new SideVariable<NDIM, double>("N_un");
        Pointer<SideVariable<NDIM, double>> N_us_var = new SideVariable<NDIM, double>("N_us");

        // Error terms.
        Pointer<SideVariable<NDIM, double>> e_un_var = new SideVariable<NDIM, double>("e_un");
        Pointer<SideVariable<NDIM, double>> e_us_var = new SideVariable<NDIM, double>("e_us");

        // Register patch data indices...
        const int un_idx = var_db->registerVariableAndContext(un_var, ctx, IntVector<NDIM>(3));
        const int us_idx = var_db->registerVariableAndContext(us_var, ctx, IntVector<NDIM>(3));
        const int thn_idx = var_db->registerVariableAndContext(thn_var, ctx, IntVector<NDIM>(3));
        const int N_un_idx = var_db->registerVariableAndContext(N_un_var, ctx, IntVector<NDIM>(1));
        const int N_us_idx = var_db->registerVariableAndContext(N_us_var, ctx, IntVector<NDIM>(1));
        const int e_un_idx = var_db->registerVariableAndContext(e_un_var, ctx, IntVector<NDIM>(1));
        const int e_us_idx = var_db->registerVariableAndContext(e_us_var, ctx, IntVector<NDIM>(1));

        // Drawing variables
        Pointer<CellVariable<NDIM, double>> draw_un_var = new CellVariable<NDIM, double>("draw_un", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_Nn_var = new CellVariable<NDIM, double>("draw_Nn", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_en_var = new CellVariable<NDIM, double>("draw_en", NDIM);
        const int draw_un_idx = var_db->registerVariableAndContext(draw_un_var, ctx);
        const int draw_Nn_idx = var_db->registerVariableAndContext(draw_Nn_var, ctx);
        const int draw_en_idx = var_db->registerVariableAndContext(draw_en_var, ctx);

        Pointer<CellVariable<NDIM, double>> draw_us_var = new CellVariable<NDIM, double>("draw_us", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_Ns_var = new CellVariable<NDIM, double>("draw_Ns", NDIM);
        Pointer<CellVariable<NDIM, double>> draw_es_var = new CellVariable<NDIM, double>("draw_es", NDIM);
        const int draw_us_idx = var_db->registerVariableAndContext(draw_us_var, ctx);
        const int draw_Ns_idx = var_db->registerVariableAndContext(draw_Ns_var, ctx);
        const int draw_es_idx = var_db->registerVariableAndContext(draw_es_var, ctx);

        Pointer<CellVariable<NDIM, double>> exact_N_un_var = new CellVariable<NDIM, double>("exact_N_un", NDIM);
        Pointer<CellVariable<NDIM, double>> exact_N_us_var = new CellVariable<NDIM, double>("exact_N_us", NDIM);
        const int exact_N_un_idx = var_db->registerVariableAndContext(exact_N_un_var, ctx);
        const int exact_N_us_idx = var_db->registerVariableAndContext(exact_N_us_var, ctx);

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity("Thn", "SCALAR", thn_idx);
        visit_data_writer->registerPlotQuantity("Un", "VECTOR", draw_un_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("Un_" + std::to_string(d), "SCALAR", draw_un_idx, d);
        }

        visit_data_writer->registerPlotQuantity("N_un", "VECTOR", draw_Nn_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("N_un_" + std::to_string(d), "SCALAR", draw_Nn_idx, d);
        }

        visit_data_writer->registerPlotQuantity("error_N_un", "VECTOR", draw_en_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("error_N_un_" + std::to_string(d), "SCALAR", draw_en_idx, d);
        }

        visit_data_writer->registerPlotQuantity("Us", "VECTOR", draw_us_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("Us_" + std::to_string(d), "SCALAR", draw_us_idx, d);
        }

        visit_data_writer->registerPlotQuantity("N_us", "VECTOR", draw_Ns_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("N_us_" + std::to_string(d), "SCALAR", draw_Ns_idx, d);
        }

        visit_data_writer->registerPlotQuantity("error_N_us", "VECTOR", draw_es_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("error_N_us_" + std::to_string(d), "SCALAR", draw_es_idx, d);
        }

        visit_data_writer->registerPlotQuantity("exact_N_un", "VECTOR", exact_N_un_idx);
        visit_data_writer->registerPlotQuantity("exact_N_us", "VECTOR", exact_N_us_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("exact_N_un_" + std::to_string(d), "SCALAR", exact_N_un_idx, d);
            visit_data_writer->registerPlotQuantity("exact_N_us_" + std::to_string(d), "SCALAR", exact_N_us_idx, d);
        }

        // Apply convective operator. This should be created before the patch hierarchy is created.
        MultiphaseConvectiveManager convec_op("convec_op",
                                              patch_hierarchy,
                                              input_db->getDatabase("ConvecOp"),
                                              un_bc_coefs,
                                              us_bc_coefs,
                                              thn_bc_coef.get());

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
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        allocate_patch_data({ un_idx,
                              us_idx,
                              N_un_idx,
                              N_us_idx,
                              e_un_idx,
                              e_us_idx,
                              thn_idx,
                              draw_un_idx,
                              draw_us_idx,
                              draw_en_idx,
                              draw_es_idx,
                              draw_Nn_idx,
                              draw_Ns_idx,
                              exact_N_un_idx,
                              exact_N_us_idx },
                            patch_hierarchy,
                            0.0,
                            coarsest_ln,
                            finest_ln);

        // Setup velocity and pressures functions.
        muParserCartGridFunction un_fcn("un", app_initializer->getComponentDatabase("un"), grid_geometry);
        muParserCartGridFunction us_fcn("us", app_initializer->getComponentDatabase("us"), grid_geometry);

        // Set up Thn functions
        muParserCartGridFunction thn_fcn("thn", app_initializer->getComponentDatabase("thn"), grid_geometry);

        // Setup exact solution functions
        muParserCartGridFunction N_un_fcn("N_un", app_initializer->getComponentDatabase("N_un"), grid_geometry);
        muParserCartGridFunction N_us_fcn("N_us", app_initializer->getComponentDatabase("N_us"), grid_geometry);

        un_fcn.setDataOnPatchHierarchy(un_idx, un_var, patch_hierarchy, 0.0);
        us_fcn.setDataOnPatchHierarchy(us_idx, us_var, patch_hierarchy, 0.0);
        thn_fcn.setDataOnPatchHierarchy(thn_idx, thn_var, patch_hierarchy, 0.0);

        N_un_fcn.setDataOnPatchHierarchy(e_un_idx, e_un_var, patch_hierarchy, 0.0);
        N_us_fcn.setDataOnPatchHierarchy(e_us_idx, e_us_var, patch_hierarchy, 0.0);

        convec_op.approximateConvectiveOperator(N_un_idx,
                                                N_us_idx,
                                                TimeSteppingType::FORWARD_EULER,
                                                0.0,
                                                1.0,
                                                un_idx,
                                                us_idx,
                                                thn_idx,
                                                un_idx,
                                                us_idx,
                                                thn_idx);

        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        // Compute error and print error norms.
        hier_sc_data_ops.subtract(e_un_idx, e_un_idx, N_un_idx);
        hier_sc_data_ops.subtract(e_us_idx, e_us_idx, N_us_idx);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(0, patch_hierarchy->getFinestLevelNumber());
        // just computes quadrature weights
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        pout << "Error in N_un :\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_un_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_un_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(e_un_idx, wgt_sc_idx) << "\n";

        pout << "Error in N_us :\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_us_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_us_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(e_us_idx, wgt_sc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(draw_un_idx, draw_un_var, un_idx, un_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_Nn_idx, draw_Nn_var, N_un_idx, N_un_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_en_idx, draw_en_var, e_un_idx, e_un_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_us_idx, draw_us_var, us_idx, us_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_Ns_idx, draw_Ns_var, N_us_idx, N_us_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(draw_es_idx, draw_es_var, e_us_idx, e_us_var, nullptr, 0.0, synch_cf_interface);

        N_un_fcn.setDataOnPatchHierarchy(exact_N_un_idx, exact_N_un_var, patch_hierarchy, 0.0, false);
        N_us_fcn.setDataOnPatchHierarchy(exact_N_us_idx, exact_N_us_var, patch_hierarchy, 0.0, false);

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        deallocate_patch_data({ un_idx,
                                us_idx,
                                N_un_idx,
                                N_us_idx,
                                e_un_idx,
                                e_us_idx,
                                thn_idx,
                                draw_un_idx,
                                draw_us_idx,
                                draw_en_idx,
                                draw_es_idx,
                                draw_Nn_idx,
                                draw_Ns_idx,
                                exact_N_un_idx,
                                exact_N_us_idx },
                              patch_hierarchy,
                              coarsest_ln,
                              finest_ln);
    } // cleanup dynamically allocated objects prior to shutdown
} // main
