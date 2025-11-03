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

#include "multiphase/MultiphaseStaggeredStokesOperator.h"

#include <ibamr/PETScKrylovStaggeredStokesSolver.h>
#include <ibamr/StaggeredStokesSolverManager.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
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

        // Grab the boundary condition objects
        std::vector<RobinBcCoefStrategy<NDIM>*> un_bc_coefs(NDIM, nullptr), us_bc_coefs(NDIM, nullptr);
        RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr;
        bool is_periodic = grid_geometry->getPeriodicShift().max() == 1;
        if (!is_periodic)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                std::string un_bc_coef_name = "un_bc_" + std::to_string(d);
                un_bc_coefs[d] =
                    new muParserRobinBcCoefs(un_bc_coef_name, input_db->getDatabase(un_bc_coef_name), grid_geometry);
                std::string us_bc_coef_name = "us_bc_" + std::to_string(d);
                us_bc_coefs[d] =
                    new muParserRobinBcCoefs(us_bc_coef_name, input_db->getDatabase(us_bc_coef_name), grid_geometry);
            }
            std::string thn_bc_name = "thn_bc";
            thn_bc_coef = new muParserRobinBcCoefs(thn_bc_name, input_db->getDatabase(thn_bc_name), grid_geometry);
        }

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

        auto thn_manager =
            std::make_unique<VolumeFractionDataManager>("ThnManager", ctx, thn_cc_var, thn_bc_coef, 1.0e-5);

        // Register patch data indices...
        const int un_sc_idx = var_db->registerVariableAndContext(un_sc_var, ctx, IntVector<NDIM>(1));
        const int us_sc_idx = var_db->registerVariableAndContext(us_sc_var, ctx, IntVector<NDIM>(1));
        const int p_cc_idx = var_db->registerVariableAndContext(p_cc_var, ctx, IntVector<NDIM>(1));
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

        Pointer<CellVariable<NDIM, double>> exact_cc_un_var = new CellVariable<NDIM, double>("exact_un", NDIM);
        Pointer<CellVariable<NDIM, double>> exact_cc_us_var = new CellVariable<NDIM, double>("exact_us", NDIM);
        Pointer<CellVariable<NDIM, double>> exact_cc_p_var = new CellVariable<NDIM, double>("exact_p");
        const int exact_cc_un_idx = var_db->registerVariableAndContext(exact_cc_un_var, ctx);
        const int exact_cc_us_idx = var_db->registerVariableAndContext(exact_cc_us_var, ctx);
        const int exact_cc_p_idx = var_db->registerVariableAndContext(exact_cc_p_var, ctx);

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity("Pressure", "SCALAR", p_cc_idx);
        visit_data_writer->registerPlotQuantity("Thn", "SCALAR", thn_manager->getCellIndex());
        visit_data_writer->registerPlotQuantity("RHS_P", "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity("error_RHS_p", "SCALAR", e_cc_idx);

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

        visit_data_writer->registerPlotQuantity("error_RHS_un", "VECTOR", draw_en_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("error_RHS_un_" + std::to_string(d), "SCALAR", draw_en_idx, d);
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

        visit_data_writer->registerPlotQuantity("error_RHS_us", "VECTOR", draw_es_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("error_RHS_us_" + std::to_string(d), "SCALAR", draw_es_idx, d);
        }

        visit_data_writer->registerPlotQuantity("exact_RHS_un", "VECTOR", exact_cc_un_idx);
        visit_data_writer->registerPlotQuantity("exact_RHS_us", "VECTOR", exact_cc_us_idx);
        visit_data_writer->registerPlotQuantity("exact_RHS_p", "SCALAR", exact_cc_p_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("exact_RHS_un_" + std::to_string(d), "SCALAR", exact_cc_un_idx, d);
            visit_data_writer->registerPlotQuantity("exact_RHS_us_" + std::to_string(d), "SCALAR", exact_cc_us_idx, d);
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
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(draw_un_idx, 0.0);
            level->allocatePatchData(draw_fn_idx, 0.0);
            level->allocatePatchData(draw_en_idx, 0.0);
            level->allocatePatchData(draw_us_idx, 0.0);
            level->allocatePatchData(draw_fs_idx, 0.0);
            level->allocatePatchData(draw_es_idx, 0.0);
            level->allocatePatchData(exact_cc_un_idx, 0.0);
            level->allocatePatchData(exact_cc_us_idx, 0.0);
            level->allocatePatchData(exact_cc_p_idx, 0.0);
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

        un_fcn.setDataOnPatchHierarchy(un_sc_idx, un_sc_var, patch_hierarchy, 0.0);
        us_fcn.setDataOnPatchHierarchy(us_sc_idx, us_sc_var, patch_hierarchy, 0.0);
        p_fcn.setDataOnPatchHierarchy(p_cc_idx, p_cc_var, patch_hierarchy, 0.0);
        thn_manager->updateVolumeFraction(thn_fcn, patch_hierarchy, 0.0, TimePoint::CURRENT_TIME);

        f_un_fcn.setDataOnPatchHierarchy(e_un_sc_idx, e_un_sc_var, patch_hierarchy, 0.0);
        f_us_fcn.setDataOnPatchHierarchy(e_us_sc_idx, e_us_sc_var, patch_hierarchy, 0.0);
        f_p_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);

        // Setup the stokes operator
        MultiphaseParameters params;
        const double C = input_db->getDouble("C");
        const double D = input_db->getDouble("D");
        params.xi = input_db->getDouble("XI");
        params.eta_n = input_db->getDouble("ETAN");
        params.eta_s = input_db->getDouble("ETAS");
        params.nu_n = params.nu_s = input_db->getDouble("NU");
        MultiphaseStaggeredStokesOperator stokes_op("stokes_op", false, params, thn_manager);
        stokes_op.setCandDCoefficients(C, D);
        stokes_op.setPhysicalBcCoefs(un_bc_coefs, us_bc_coefs, nullptr);

        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_un_helper = new StaggeredStokesPhysicalBoundaryHelper();
        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_us_helper = new StaggeredStokesPhysicalBoundaryHelper();
        stokes_op.setPhysicalBoundaryHelper(bc_un_helper, bc_us_helper);

        stokes_op.initializeOperatorState(u_vec, f_vec);

        // Apply the operator
        stokes_op.apply(u_vec, f_vec);

        // Print out f_un_sc_idx and f_us_sc_idx.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                if (pgeom->getTouchesRegularBoundary())
                {
                    Pointer<SideData<NDIM, double>> fun_data = patch->getPatchData(f_un_sc_idx);
                    Pointer<SideData<NDIM, double>> fus_data = patch->getPatchData(f_us_sc_idx);
                    Pointer<SideData<NDIM, double>> exact_us_data = patch->getPatchData(e_us_sc_idx);
                    Pointer<SideData<NDIM, double>> exact_un_data = patch->getPatchData(e_un_sc_idx);
                    const tbox::Array<BoundaryBox<NDIM>>& physical_codim1_boxes =
                        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                    const int num_bdry_boxes = physical_codim1_boxes.size();
                    for (int n = 0; n < num_bdry_boxes; ++n)
                    {
                        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                        const int location_index = bdry_box.getLocationIndex();
                        const int bdry_normal_axis = location_index / 2;
                        const int upper_lower = ((location_index % 2) + 1) % 2;
                        for (CellIterator<NDIM> ci(bdry_box.getBox()); ci; ci++)
                        {
                            const CellIndex<NDIM>& cell_idx = ci();
                            SideIndex<NDIM> idx(cell_idx, bdry_normal_axis, upper_lower);
                            pout << "On idx " << cell_idx << " on bdry " << bdry_normal_axis << " and "
                                 << (upper_lower == 0 ? "lower" : "upper") << " side\n";
                            pout << "Un: " << (*fun_data)(idx) << "\n";
                            pout << "Exact un: " << (*exact_un_data)(idx) << "\n";
                            pout << "Us: " << (*fus_data)(idx) << "\n";
                            pout << "Exact us: " << (*exact_us_data)(idx) << "\n";
                        }
                    }
                }
            }
        }

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false),  // numerical
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
        pout << "Error in RHS_un :\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_un_sc_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_un_sc_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(e_un_sc_idx, wgt_sc_idx) << "\n";

        pout << "Error in RHS_us :\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_us_sc_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_us_sc_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(e_us_sc_idx, wgt_sc_idx) << "\n";

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        pout << "Error in RHS_p :\n"
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

        f_un_fcn.setDataOnPatchHierarchy(exact_cc_un_idx, exact_cc_un_var, patch_hierarchy, 0.0, false);
        f_us_fcn.setDataOnPatchHierarchy(exact_cc_us_idx, exact_cc_us_var, patch_hierarchy, 0.0, false);
        f_p_fcn.setDataOnPatchHierarchy(exact_cc_p_idx, exact_cc_p_var, patch_hierarchy, 0.0, false);

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
            level->deallocatePatchData(f_cc_idx);
            level->deallocatePatchData(e_cc_idx);
            level->deallocatePatchData(draw_un_idx);
            level->deallocatePatchData(draw_fn_idx);
            level->deallocatePatchData(draw_en_idx);
            level->deallocatePatchData(draw_us_idx);
            level->deallocatePatchData(draw_fs_idx);
            level->deallocatePatchData(draw_es_idx);
            level->deallocatePatchData(exact_cc_un_idx);
            level->deallocatePatchData(exact_cc_us_idx);
            level->deallocatePatchData(exact_cc_p_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
