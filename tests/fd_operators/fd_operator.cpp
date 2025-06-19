#include "multiphase/MultiphaseStaggeredStokesOperator.h"
#include "multiphase/fd_operators.h"

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

enum class TestType
{
    AccumulateMomentumCC,
    AccumulateMomentumVD,
    AccumulateMomentumCCNoPressure,
    AccumulateMomentumVDNoPressure,
    Coincompressibility,
    MultiphaseGradient,
    UNKNOWN
};
// Convert enum to string
inline std::string
to_string(TestType type)
{
    switch (type)
    {
    case TestType::AccumulateMomentumCC:
        return "AccumulateMomentumCC";
    case TestType::AccumulateMomentumVD:
        return "AccumulateMomentumVD";
    case TestType::AccumulateMomentumCCNoPressure:
        return "AccumulateMomentumCCNoPressure";
    case TestType::AccumulateMomentumVDNoPressure:
        return "AccumulateMomentumVDNoPressure";
    case TestType::Coincompressibility:
        return "Coincompressibility";
    case TestType::MultiphaseGradient:
        return "MultiphaseGradient";
    case TestType::UNKNOWN:
        return "UNKNOWN";
    default:
        return "UNKNOWN";
    }
}

// Convert string to enum
inline TestType
test_type_from_string(const std::string& str)
{
    if (str == "AccumulateMomentumCC")
        return TestType::AccumulateMomentumCC;
    else if (str == "AccumulateMomentumVD")
        return TestType::AccumulateMomentumVD;
    else if (str == "AccumulateMomentumCCNoPressure")
        return TestType::AccumulateMomentumCCNoPressure;
    else if (str == "AccumulateMomentumVDNoPressure")
        return TestType::AccumulateMomentumVDNoPressure;
    else if (str == "Coincompressibility")
        return TestType::Coincompressibility;
    else if (str == "MultiphaseGradient")
        return TestType::MultiphaseGradient;
    else if (str == "UNKNOWN")
        return TestType::UNKNOWN;
    else
        return TestType::UNKNOWN;
}

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

        // Create the boundary condition objects
        std::vector<RobinBcCoefStrategy<NDIM>*> un_bc_coefs(NDIM, nullptr), us_bc_coefs(NDIM, nullptr);
        std::vector<RobinBcCoefStrategy<NDIM>*> p_bc_coefs(NDIM, nullptr);
        RobinBcCoefStrategy<NDIM>* thn_bc_coef = nullptr;

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // State variables: Velocity and pressure.
        Pointer<SideVariable<NDIM, double>> un_var = new SideVariable<NDIM, double>("un");
        Pointer<SideVariable<NDIM, double>> us_var = new SideVariable<NDIM, double>("us");
        Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");

        // variable coefficient: theta_n
        Pointer<CellVariable<NDIM, double>> thn_var = new CellVariable<NDIM, double>("thn");

        // Drag coefficient
        Pointer<SideVariable<NDIM, double>> xi_var = new SideVariable<NDIM, double>("xi");

        // Results of operator "forces" and "divergence"
        Pointer<SideVariable<NDIM, double>> f_un_var = new SideVariable<NDIM, double>("f_un");
        Pointer<SideVariable<NDIM, double>> f_us_var = new SideVariable<NDIM, double>("f_us");
        Pointer<CellVariable<NDIM, double>> f_var = new CellVariable<NDIM, double>("f");

        // Error terms.
        Pointer<SideVariable<NDIM, double>> e_un_var = new SideVariable<NDIM, double>("e_un");
        Pointer<SideVariable<NDIM, double>> e_us_var = new SideVariable<NDIM, double>("e_us");
        Pointer<CellVariable<NDIM, double>> e_var = new CellVariable<NDIM, double>("e");

        // Create volume fraction manager
        auto thn_manager = std::make_unique<VolumeFractionDataManager>(
            "ThnManager", thn_var, ctx, thn_bc_coef, 1.0e-5 /*regularize_thn*/);

        // Register patch data indices...
        const int un_idx = var_db->registerVariableAndContext(un_var, ctx, IntVector<NDIM>(1));
        const int us_idx = var_db->registerVariableAndContext(us_var, ctx, IntVector<NDIM>(1));
        const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
        const int xi_idx = var_db->registerVariableAndContext(xi_var, ctx, IntVector<NDIM>(1));
        const int f_idx = var_db->registerVariableAndContext(f_var, ctx, IntVector<NDIM>(1));
        const int f_un_idx = var_db->registerVariableAndContext(f_un_var, ctx, IntVector<NDIM>(1));
        const int f_us_idx = var_db->registerVariableAndContext(f_us_var, ctx, IntVector<NDIM>(1));
        const int e_idx = var_db->registerVariableAndContext(e_var, ctx, 0);
        const int e_un_idx = var_db->registerVariableAndContext(e_un_var, ctx, 0);
        const int e_us_idx = var_db->registerVariableAndContext(e_us_var, ctx, 0);
        std::set<int> patch_idxs{ un_idx, us_idx, p_idx, xi_idx, f_idx, f_un_idx, f_us_idx, e_idx, e_un_idx, e_us_idx };

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
        allocate_patch_data(patch_idxs, patch_hierarchy, 0.0, 0, patch_hierarchy->getFinestLevelNumber());

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

        un_fcn.setDataOnPatchHierarchy(un_idx, un_var, patch_hierarchy, 0.0);
        us_fcn.setDataOnPatchHierarchy(us_idx, us_var, patch_hierarchy, 0.0);
        p_fcn.setDataOnPatchHierarchy(p_idx, p_var, patch_hierarchy, 0.0);
        thn_manager->updateVolumeFraction(thn_fcn, patch_hierarchy, 0.0, TimePoint::CURRENT_TIME);

        bool using_var_xi = input_db->getBool("USING_VAR_XI");

        // Setup the stokes operator
        MultiphaseParameters params;
        if (using_var_xi)
        {
            params.xi_idx = xi_idx;
            muParserCartGridFunction xi_fcn("xi", app_initializer->getComponentDatabase("xi_fcn"), grid_geometry);
            xi_fcn.setDataOnPatchHierarchy(xi_idx, xi_var, patch_hierarchy, 0.0);
        }
        else
        {
            params.xi = input_db->getDouble("XI");
            params.nu_n = params.nu_s = input_db->getDouble("NU");
        }
        params.eta_n = input_db->getDouble("ETAN");
        params.eta_s = input_db->getDouble("ETAS");
        params.lambda_n = params.eta_n;
        params.lambda_s = params.eta_s;

        // Make sure boundaries are filled in correctly
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comps;
        ghost_cell_comps.push_back(ITC(un_idx, "NONE", true, "NONE", "LINEAR", false, un_bc_coefs));
        ghost_cell_comps.push_back(ITC(us_idx, "NONE", true, "NONE", "LINEAR", false, us_bc_coefs));
        ghost_cell_comps.push_back(ITC(p_idx, "CONSERVATIVE_LINEAR_REFINE", true, "NONE", "LINEAR", false, p_bc_coefs));
        HierarchyGhostCellInterpolation hier_ghost_fill;
        hier_ghost_fill.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_fill.fillData(0.0);

        TestType test_type = test_type_from_string(input_db->getString("TEST_TYPE"));
        switch (test_type)
        {
        case TestType::MultiphaseGradient:
            multiphase_grad_on_hierarchy(
                *patch_hierarchy, f_un_idx, f_us_idx, thn_manager->getSideIndex(), p_idx, 1.0, false);
            break;
        case TestType::Coincompressibility:
            applyCoincompressibility(*patch_hierarchy, f_idx, un_idx, us_idx, thn_manager->getSideIndex(), 1.0);
            break;
        case TestType::AccumulateMomentumVDNoPressure:
            accumulateMomentumWithoutPressureVariableDrag(
                *patch_hierarchy, f_un_idx, f_us_idx, un_idx, us_idx, *thn_manager, params, 1.0, 1.0);
            break;
        case TestType::AccumulateMomentumVD:
            accumulateMomentumForcesVariableDrag(
                *patch_hierarchy, f_un_idx, f_us_idx, p_idx, un_idx, us_idx, *thn_manager, params, 1.0, 1.0, 1.0);
            break;
        case TestType::AccumulateMomentumCC:
            accumulateMomentumForcesConstantCoefficient(
                *patch_hierarchy, f_un_idx, f_us_idx, p_idx, un_idx, us_idx, *thn_manager, params, 1.0, 1.0, 1.0);
            break;
        case TestType::AccumulateMomentumCCNoPressure:
            accumulateMomentumWithoutPressureConstantCoefficient(
                *patch_hierarchy, f_un_idx, f_us_idx, un_idx, us_idx, *thn_manager, params, 1.0, 1.0);
            break;
        case TestType::UNKNOWN:
        default:
            TBOX_ERROR("Unknown test type\n");
        }

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(0, patch_hierarchy->getFinestLevelNumber());
        // just computes quadrature weights
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        switch (test_type)
        {
        case TestType::AccumulateMomentumVDNoPressure:
        case TestType::AccumulateMomentumCCNoPressure:
        case TestType::AccumulateMomentumCC:
        case TestType::AccumulateMomentumVD:
        case TestType::MultiphaseGradient:
        {
            f_un_fcn.setDataOnPatchHierarchy(e_un_idx, e_un_var, patch_hierarchy, 0.0);
            f_us_fcn.setDataOnPatchHierarchy(e_us_idx, e_us_var, patch_hierarchy, 0.0);
            hier_sc_data_ops.subtract(e_un_idx, e_un_idx, f_un_idx);
            hier_sc_data_ops.subtract(e_us_idx, e_us_idx, f_us_idx);
            pout << "Error in RHS_un :\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_un_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_un_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << hier_sc_data_ops.maxNorm(e_un_idx, wgt_sc_idx) << "\n";

            pout << "Error in RHS_us :\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(e_us_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << hier_sc_data_ops.L2Norm(e_us_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << hier_sc_data_ops.maxNorm(e_us_idx, wgt_sc_idx) << "\n";
            break;
        }
        case TestType::Coincompressibility:
        {
            f_p_fcn.setDataOnPatchHierarchy(e_idx, e_var, patch_hierarchy, 0.0);
            hier_cc_data_ops.subtract(e_idx, e_idx, f_idx);
            pout << "Error in RHS_p :\n"
                 << "  L1-norm:  " << hier_cc_data_ops.L1Norm(e_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm:  " << hier_cc_data_ops.L2Norm(e_idx, wgt_cc_idx) << "\n"
                 << "  max-norm: " << hier_cc_data_ops.maxNorm(e_idx, wgt_cc_idx) << "\n";
            break;
        }
        case TestType::UNKNOWN:
        default:
            TBOX_ERROR("Unknown test type\n");
        }

        // Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
        // for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        // {
        //     Pointer<Patch<NDIM>> patch = level->getPatch(p());
        //     Pointer<SideData<NDIM, double>> un_data = patch->getPatchData(un_idx);
        //     Pointer<SideData<NDIM, double>> Gun_data = patch->getPatchData(f_un_idx);
        //     const Box<NDIM>& patch_box = patch->getBox();
        //     const Box<NDIM>& ghost_box = un_data->getGhostBox();
        //     for (SideIterator<NDIM> si(ghost_box, 0); si; si++)
        //     {
        //         const SideIndex<NDIM>& idx = si();
        //         if (!(patch_box.contains(idx.toCell(1)) && patch_box.contains(idx.toCell(0))))
        //         {
        //             pout << "On index " << idx << " value is " << (*un_data)(idx) << "\n";
        //         }
        //         if ((!patch_box.contains(idx.toCell(0))) != !(patch_box.contains(idx.toCell(1))))
        //         {
        //             pout << "Gun value is " << (*Gun_data)(idx) << "\n";
        //         }
        //     }
        // }

        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        // Deallocate level data
        deallocate_patch_data(patch_idxs, patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        thn_manager->deallocateData();
    } // cleanup dynamically allocated objects prior to shutdown
} // main
