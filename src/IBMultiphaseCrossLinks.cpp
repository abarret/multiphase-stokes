/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBBeamForceSpec.h"
#include "ibamr/IBSpringForceFunctions.h"
#include "ibamr/IBSpringForceSpec.h"
#include "ibamr/IBTargetPointForceSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscmat.h"
#include "petscvec.h"
#include <petsclog.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Local includes
#include "IBMultiphaseCrossLinks.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMultiphaseCrossLinks::IBMultiphaseCrossLinks(LDataManager* other_data_manager, const double kappa)
    : d_other_data_manager(other_data_manager), d_kappa(kappa)
{
    // intentionally blank
    return;
} // IBMultiphaseCrossLinks

void
IBMultiphaseCrossLinks::initializeLevelData(const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                            const int level_number,
                                            const double init_data_time,
                                            const bool initial_time,
                                            LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    // We don't need ghost points or any initialization steps here.
    d_is_initialized.resize(std::max(level_number + 1, static_cast<int>(d_is_initialized.size())));
    d_is_initialized[level_number] = true;
    return;
} // initializeLevelData

void
IBMultiphaseCrossLinks::computeLagrangianForce(Pointer<LData> F_data,
                                               Pointer<LData> X_data,
                                               Pointer<LData> U_data,
                                               const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                               const int level_number,
                                               const double data_time,
                                               LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    computeCrossSpringForce(F_data, X_data, U_data, hierarchy, level_number, data_time, l_data_manager);
    return;
} // computeLagrangianForce

void
IBMultiphaseCrossLinks::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const Pointer<PatchHierarchy<NDIM>> /*hierarchy*/,
    const int level_number,
    LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    TBOX_ERROR("Not currently implemented\n");
    return;
} // computeLagrangianForceJacobianNonzeroStructure

void
IBMultiphaseCrossLinks::computeLagrangianForceJacobian(Mat& J_mat,
                                                       MatAssemblyType assembly_type,
                                                       const double X_coef,
                                                       Pointer<LData> X_data,
                                                       const double U_coef,
                                                       Pointer<LData> /*U_data*/,
                                                       const Pointer<PatchHierarchy<NDIM>> /*hierarchy*/,
                                                       const int level_number,
                                                       const double /*data_time*/,
                                                       LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    TBOX_ERROR("Not currently implemented\n");
    return;
} // computeLagrangianForceJacobian

double
IBMultiphaseCrossLinks::computeLagrangianEnergy(Pointer<LData> /*X_data*/,
                                                Pointer<LData> /*U_data*/,
                                                const Pointer<PatchHierarchy<NDIM>> /*hierarchy*/,
                                                const int level_number,
                                                const double /*data_time*/,
                                                LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return 0.0;

    // Compute the energy.
    TBOX_ERROR("not currently implemented\n");
    return std::numeric_limits<double>::quiet_NaN();
} // computeLagrangianEnergy

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBMultiphaseCrossLinks::computeCrossSpringForce(Pointer<LData> F_data,
                                                Pointer<LData> X_data,
                                                Pointer<LData> /*U_data*/,
                                                const Pointer<PatchHierarchy<NDIM>> /*hierarchy*/,
                                                const int level_number,
                                                const double /*data_time*/,
                                                LDataManager* const l_data_manager)
{
    double max_stretch = 0.0;
    // Pull out relevant data.
    double* const F_node = F_data->getLocalFormVecArray()->data();
    const double* const X_node = X_data->getLocalFormVecArray()->data();

    // We also need the other sites data
    Pointer<LData> X_other_data = d_other_data_manager->getLData(d_other_data_manager->POSN_DATA_NAME, level_number);
    double* const X_other_node = X_other_data->getLocalFormVecArray()->data();

    // Loop through local indices and compute force
    const std::vector<LNode*>& nodes = l_data_manager->getLMesh(level_number)->getLocalNodes();

    // Note that although the initial Lagrangian indexing is the same, slight differences in structures can mean that
    // the PETSc ordering can be different. Here we gather the Lagrangian indices on which we are calculating forces and
    // determine the PETSc ordering of the other structural data.
    // TODO: Remove the assumption that the initial Lagrangian indexing is the same. We need to store a map between
    // initial Lagrangian indices.
#ifndef NDEBUG
    TBOX_ASSERT(l_data_manager->getNumberOfNodes(level_number) >= nodes.size());
#endif
    std::vector<int> o_petsc_idxs(l_data_manager->getNumberOfNodes(level_number), -1);
    for (int i = 0; i < nodes.size(); ++i)
        o_petsc_idxs[nodes[i]->getLagrangianIndex()] = nodes[i]->getLagrangianIndex();
    d_other_data_manager->mapLagrangianToPETSc(o_petsc_idxs, level_number);

    for (const auto& node : nodes)
    {
        const int lag_idx = node->getLagrangianIndex();
        const int petsc_idx = node->getLocalPETScIndex();
        const int o_petsc_idx = o_petsc_idxs[lag_idx];

        Eigen::Map<const VectorNd> X(&X_node[petsc_idx * NDIM]);
        Eigen::Map<const VectorNd> X_o(&X_other_node[o_petsc_idx * NDIM]);
        Eigen::Map<VectorNd> F(&F_node[petsc_idx * NDIM]);
        F += -d_kappa * (X - X_o);
        if (d_log_point_distances) max_stretch = std::max(max_stretch, (X - X_o).norm());
    }

    if (d_log_point_distances) plog << "max_stretch: " << max_stretch << "\n";

    F_data->restoreArrays();
    X_data->restoreArrays();
    X_other_data->restoreArrays();
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
