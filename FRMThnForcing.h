// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_FRMThnForcing
#define included_FRMThnForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class FRMThnForcing provides forcing for the momentum equations
 * based on the Boussinesq approximation to the variable-density incompressible
 * Navier-Stokes equations.
 */
class FRMThnForcing : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    FRMThnForcing(Pointer<Variable<NDIM>> thn_var,
                  Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
                  bool using_thn);

    /*!
     * \brief Empty destructor.
     */
    ~FRMThnForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete FRMThnForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 Pointer<Variable<NDIM>> var,
                                 Pointer<PatchHierarchy<NDIM>> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<Variable<NDIM>> var,
                        Pointer<Patch<NDIM>> patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM>> patch_level = Pointer<PatchLevel<NDIM>>(NULL));

    //\}

private:
    FRMThnForcing();

    FRMThnForcing(const FRMThnForcing& from);

    FRMThnForcing& operator=(const FRMThnForcing& that);

    Pointer<Variable<NDIM>> d_thn_var;
    Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;
    bool d_using_thn = true;
};

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_FRMThnForcing
