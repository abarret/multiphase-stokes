#include "multiphase/FullFACPreconditioner.h"
#include "multiphase/MultiphaseStaggeredStokesBoxRelaxationFACOperator.h"
#include "multiphase/MultiphaseStaggeredStokesOperator.h"
#include "multiphase/utility_functions.h"

#include <ibamr/StaggeredStokesSolverManager.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/app_namespaces.h>

#include "ibtk/CartCellDoubleQuadraticRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LinearOperator.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonSolverManager.h>
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

#include <chrono>

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

        // State variables: Pressure.
        Pointer<CellVariable<NDIM, double>> p_cc_var = new CellVariable<NDIM, double>("p");

        // variable coefficient: theta_n and thn^2 + ths^2
        Pointer<CellVariable<NDIM, double>> thn_cc_var = new CellVariable<NDIM, double>("thn");
        Pointer<SideVariable<NDIM, double>> thn_ths_sq_var = new SideVariable<NDIM, double>("thn_ths_sq_sc");

        // Results of laplace operator 
        Pointer<CellVariable<NDIM, double>> f_cc_var = new CellVariable<NDIM, double>("f");

        // Error terms.
        Pointer<CellVariable<NDIM, double>> e_cc_var = new CellVariable<NDIM, double>("e");

        // Register patch data indices...
        const int p_cc_idx = var_db->registerVariableAndContext(p_cc_var, ctx, IntVector<NDIM>(1));
        const int thn_cc_idx =
            var_db->registerVariableAndContext(thn_cc_var, ctx, IntVector<NDIM>(1)); // 1 layer of ghost cells
        const int thn_ths_sq_idx = var_db->registerVariableAndContext(thn_ths_sq_var, ctx, IntVector<NDIM>(1)); 
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity("Pressure", "SCALAR", p_cc_idx);
        visit_data_writer->registerPlotQuantity("Thn", "SCALAR", thn_cc_idx);
        visit_data_writer->registerPlotQuantity("RHS_P", "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity("error_p", "SCALAR", e_cc_idx);

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
            level->allocatePatchData(p_cc_idx, 0.0);
            level->allocatePatchData(thn_cc_idx, 0.0);
            level->allocatePatchData(thn_ths_sq_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("p", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(p_cc_var, p_cc_idx, h_cc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        // Setup pressures function.
        muParserCartGridFunction p_fcn("p", app_initializer->getComponentDatabase("p"), grid_geometry);

        // Set up Thn function
        muParserCartGridFunction thn_fcn("thn", app_initializer->getComponentDatabase("thn"), grid_geometry);

        // Setup exact solution functions
        muParserCartGridFunction f_p_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

        f_p_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);
        thn_fcn.setDataOnPatchHierarchy(thn_cc_idx, thn_cc_var, patch_hierarchy, 0.0);
        p_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);
        
        // This computes Thn^2 + Ths^2
        SAMRAI::hier::IntVector<NDIM> xp(1, 0), yp(0, 1);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {   
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx(); // dx[0] -> x, dx[1] -> y
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> thn_data = patch->getPatchData(thn_cc_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> thn_ths_sq_data = patch->getPatchData(thn_ths_sq_idx); 

                for(int axis = 0; axis < NDIM; axis++) 
                {
                    for (SAMRAI::pdat::SideIterator<NDIM> si(patch->getBox(), axis); si; si++) 
                    {
                        const SAMRAI::pdat::SideIndex<NDIM>& idx = si(); // if axis = 0, (i-1/2,j)

                        SAMRAI::pdat::CellIndex<NDIM> idx_c_low = idx.toCell(0);   // (i-1,j)
                        SAMRAI::pdat::CellIndex<NDIM> idx_c_up = idx.toCell(1);    // (i,j)
                        // thn at sides
                        double thn_lower = 0.5 * ((*thn_data)(idx_c_low) + (*thn_data)(idx_c_up)); // thn(i-1/2,j)
                        double thn_sq = thn_lower * thn_lower;
                        double ths_sq = convertToThs(thn_lower) * convertToThs(thn_lower);
                        (*thn_ths_sq_data)(idx) = thn_sq + ths_sq;
                    }
                }
            }
        }

        // Setup the Poisson solver.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCZero();
        poisson_spec.setDPatchDataId(thn_ths_sq_idx); 
        RobinBcCoefStrategy<NDIM>* bc_coef = nullptr;
        CCLaplaceOperator laplace_op("laplace_op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoef(bc_coef);
        laplace_op.initializeOperatorState(u_vec, f_vec);
        
        // Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
        // stokes_op->setPhysicalBoundaryHelper(bc_helper);
        // stokes_op->setThnIdx(thn_cc_idx);

        string solver_type = input_db->getString("solver_type");
        Pointer<Database> solver_db = input_db->getDatabase("solver_db");
        string precond_type = input_db->getString("precond_type");
        Pointer<Database> precond_db = input_db->getDatabase("precond_db");
        Pointer<PoissonSolver> poisson_solver = CCPoissonSolverManager::getManager()->allocateSolver(
            solver_type, "poisson_solver", solver_db, "", precond_type, "poisson_precond", precond_db, "");
        poisson_solver->setPoissonSpecifications(poisson_spec);
        poisson_solver->setPhysicalBcCoef(bc_coef);
        poisson_solver->initializeSolverState(u_vec, f_vec);

        // Solve -L*u = f.
        u_vec.setToScalar(0.0);
        poisson_solver->solveSystem(u_vec, f_vec);

        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(0, patch_hierarchy->getFinestLevelNumber());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        // This subtracts the projection of solution on the nullspace.
        // The nullspace consists of constants.
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
                patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        double integral = hier_cc_data_ops.integral(p_cc_idx, wgt_cc_idx);
        hier_cc_data_ops.addScalar(p_cc_idx, p_cc_idx, -1.0 * integral);

        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false),  // numerical
                       Pointer<SAMRAIVectorReal<NDIM, double>>(&e_vec, false)); // analytical
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";

        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        pout << "Error in p :\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(e_cc_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(e_cc_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(e_cc_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        std::ofstream out("output");
        out << "|e|_2  = " << e_vec.L2Norm() << "\n";
        out << "|e|_1  = " << e_vec.L1Norm() << "\n";
        out << "|e|_oo = " << e_vec.maxNorm() << "\n\n";

        out << "Error in p :\n"
            << "  L1-norm:  " << hier_cc_data_ops.L1Norm(e_cc_idx, wgt_cc_idx) << "\n"
            << "  L2-norm:  " << hier_cc_data_ops.L2Norm(e_cc_idx, wgt_cc_idx) << "\n"
            << "  max-norm: " << hier_cc_data_ops.maxNorm(e_cc_idx, wgt_cc_idx) << "\n"
            << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        // Compute discrete divergence.
        // Fill ghost cells for theta
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_fill_specs = { ITC(
            thn_cc_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR") };
        HierarchyGhostCellInterpolation ghost_fill;
        ghost_fill.initializeOperatorState(ghost_fill_specs, patch_hierarchy);
        ghost_fill.fillData(0.0);

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 1, 0.0);
        // Deallocate level data
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(p_cc_idx);
            level->deallocatePatchData(thn_cc_idx);
            level->deallocatePatchData(thn_ths_sq_idx);
            level->deallocatePatchData(f_cc_idx);
            level->deallocatePatchData(e_cc_idx);
        }

        // Write timer data
        TimerManager::getManager()->print(plog);
    } // cleanup dynamically allocated objects prior to shutdown
} // main
