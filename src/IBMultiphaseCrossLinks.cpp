/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBBeamForceSpec.h"
#include "ibamr/IBSpringForceFunctions.h"
#include "ibamr/IBSpringForceSpec.h"
#include "ibamr/IBTargetPointForceSpec.h"
#include "ibamr/app_namespaces.h" // IWYU pragma: keep

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

#include "multiphase/IBMultiphaseCrossLinks.h"

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace multiphase
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMultiphaseCrossLinks::IBMultiphaseCrossLinks(Pointer<IBMethod> ibn_ops,
                                               Pointer<IBMethod> ibs_ops,
                                               Pointer<PatchHierarchy<NDIM>> hierarchy,
                                               const double kappa,
                                               const double eta,
                                               std::function<VectorNd(double)> vel_fcn)
    : MultiphaseCrossLinksStrategy(),
      d_ibn_ops(ibn_ops),
      d_ibs_ops(ibs_ops),
      d_kappa(kappa),
      d_eta(eta),
      d_vel_fcn(vel_fcn),
      d_hierarchy(hierarchy)
{
    // intentionally blank
    return;
} // IBMultiphaseCrossLinks

void
IBMultiphaseCrossLinks::doComputeLagrangianForce(const double time, IBTK::TimePoint time_pt)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    LDataManager* ibn_manager = d_ibn_ops->getLDataManager();
    LDataManager* ibs_manager = d_ibs_ops->getLDataManager();
    d_Fn_data.resize(finest_ln - coarsest_ln + 1);
    d_Fs_data.resize(finest_ln - coarsest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        static const bool maintain_data = false;
        d_Fn_data[ln] = ibn_manager->createLData("Fn_data", ln, NDIM, maintain_data);
        d_Fs_data[ln] = ibs_manager->createLData("Fs_data", ln, NDIM, maintain_data);
    }

    std::vector<Pointer<LData>>* Xn_data;
    std::vector<Pointer<LData>>* Un_data;
    bool* Xn_needs_ghost_fill;
    d_ibn_ops->getPositionData(&Xn_data, &Xn_needs_ghost_fill, time_pt);
    d_ibn_ops->getVelocityData(&Un_data, time_pt);

    std::vector<Pointer<LData>>* Xs_data;
    std::vector<Pointer<LData>>* Us_data;
    bool* Xs_needs_ghost_fill;

    d_ibs_ops->getPositionData(&Xs_data, &Xs_needs_ghost_fill, time_pt);
    d_ibs_ops->getVelocityData(&Us_data, time_pt);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        doComputeLagrangianForce(
            d_Fn_data[ln], d_Fs_data[ln], (*Xn_data)[ln], (*Un_data)[ln], (*Xs_data)[ln], (*Us_data)[ln], time, ln);
    }

    return;
} // computeLagrangianForce

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBMultiphaseCrossLinks::doSpreadForce(const int f_data_idx,
                                      RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                                      const std::vector<Pointer<RefineSchedule<NDIM>>>& f_prolongation_scheds,
                                      const double data_time,
                                      const bool spread_network,
                                      IBTK::TimePoint time_pt)
{
    std::vector<Pointer<LData>>* X_data;
    bool* X_needs_ghost_fill;
    static bool F_needs_ghost_fill = true;
    if (spread_network)
    {
        d_ibn_ops->getPositionData(&X_data, &X_needs_ghost_fill, time_pt);
        d_ibn_ops->getLDataManager()->spread(f_data_idx,
                                             d_Fn_data,
                                             *X_data,
                                             f_phys_bdry_op,
                                             f_prolongation_scheds,
                                             data_time,
                                             F_needs_ghost_fill,
                                             *X_needs_ghost_fill);
    }
    else
    {
        d_ibs_ops->getPositionData(&X_data, &X_needs_ghost_fill, time_pt);
        d_ibs_ops->getLDataManager()->spread(f_data_idx,
                                             d_Fs_data,
                                             *X_data,
                                             f_phys_bdry_op,
                                             f_prolongation_scheds,
                                             data_time,
                                             F_needs_ghost_fill,
                                             *X_needs_ghost_fill);
    }
}

void
IBMultiphaseCrossLinks::doComputeLagrangianForce(Pointer<LData>& Fn_data,
                                                 Pointer<LData>& Fs_data,
                                                 Pointer<LData>& Xn_data,
                                                 Pointer<LData>& Un_data,
                                                 Pointer<LData>& Xs_data,
                                                 Pointer<LData>& Us_data,
                                                 const double time,
                                                 const int ln)
{
    LDataManager* ibn_manager = d_ibn_ops->getLDataManager();
    LDataManager* ibs_manager = d_ibs_ops->getLDataManager();
    if (ibn_manager->getNumberOfLocalNodes(ln) == 0) return;
    double max_stretch = 0.0;
    // Pull out relevant data.
    double* const Fn_node = Fn_data->getLocalFormVecArray()->data();
    const double* const Xn_node = Xn_data->getLocalFormVecArray()->data();
    const double* const Un_node = Un_data->getLocalFormVecArray()->data();

    // We also need the other sites data
    double* const Fs_node = Fs_data->getLocalFormVecArray()->data();
    const double* const Xs_node = Xs_data->getLocalFormVecArray()->data();
    const double* const Us_node = Us_data->getLocalFormVecArray()->data();

    // Loop through local indices and compute forces on the network's nodes
    const std::vector<LNode*>& nodes = ibn_manager->getLMesh(ln)->getLocalNodes();

    // Note that although the initial Lagrangian indexing is the same, slight differences in structures can mean that
    // the PETSc ordering can be different. Here we gather the Lagrangian indices on which we are calculating forces and
    // determine the PETSc ordering of the other structural data.
    // TODO: Remove the assumption that the initial Lagrangian indexing is the same. We need to store a map between
    // initial Lagrangian indices.
    // TODO: We also assume that network and solvent pairs of indices are available on the same process.
#ifndef NDEBUG
    TBOX_ASSERT(ibs_manager->getNumberOfNodes(ln) >= nodes.size());
#endif
    std::vector<int> solvent_petsc_idxs(ibs_manager->getNumberOfNodes(ln), -1);
    for (size_t i = 0; i < nodes.size(); ++i)
        solvent_petsc_idxs[nodes[i]->getLagrangianIndex()] = nodes[i]->getLagrangianIndex();
    ibs_manager->mapLagrangianToPETSc(solvent_petsc_idxs, ln);

    for (const auto& node : nodes)
    {
        const int lag_idx = node->getLagrangianIndex();
        const int network_petsc_idx = node->getLocalPETScIndex();
        const int solvent_petsc_idx = solvent_petsc_idxs[lag_idx];

        Eigen::Map<const VectorNd> Xn(&Xn_node[network_petsc_idx * NDIM]);
        Eigen::Map<const VectorNd> Xs(&Xs_node[solvent_petsc_idx * NDIM]);
        Eigen::Map<const VectorNd> Un(&Un_node[network_petsc_idx * NDIM]);
        Eigen::Map<const VectorNd> Us(&Us_node[solvent_petsc_idx * NDIM]);
        Eigen::Map<VectorNd> Fn(&Fn_node[network_petsc_idx * NDIM]);
        Eigen::Map<VectorNd> Fs(&Fs_node[solvent_petsc_idx * NDIM]);
        VectorNd vel = d_vel_fcn(time);
        Fn = d_kappa * (Xs - Xn) + d_eta * (vel - Un);
        Fs = d_kappa * (Xn - Xs) + d_eta * (vel - Us);
        max_stretch = std::max(max_stretch, (Xn - Xs).norm());
    }

    pout << "max_stretch: " << max_stretch << "\n";

    Fn_data->restoreArrays();
    Xn_data->restoreArrays();
    Fs_data->restoreArrays();
    Xs_data->restoreArrays();
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace multiphase

//////////////////////////////////////////////////////////////////////////////
