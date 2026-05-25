#ifndef included_multiphase_MultiphaseStaggeredStokesBlockFACOperator
#define included_multiphase_MultiphaseStaggeredStokesBlockFACOperator

#include <ibtk/config.h>

#include "multiphase/MultiphaseParameters.h"
#include "multiphase/VolumeFractionDataManager.h"

#include "ibtk/FACPreconditionerStrategy.h"

#include "PoissonSpecifications.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
template <int DIM>
class PatchHierarchy;
}
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
template <int DIM, class TYPE>
class SideVariable;
}
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
template <int DIM, class TYPE>
class SAMRAIVectorReal;
}
namespace xfer
{
template <int DIM>
class RefinePatchStrategy;
}
} // namespace SAMRAI

namespace IBTK
{
class PoissonFACPreconditionerStrategy;
}

namespace multiphase
{
class MultiphaseStaggeredStokesBoxRelaxationFACOperator;
class MultiphaseStaggeredVelocityAsymmetricFACOperator;
class MultiphaseStaggeredVelocityBlockFACOperator;

/*!
 * \brief FAC strategy for the full staggered Stokes system that uses the
 * approximate block factorization from
 * MultiphaseStaggeredStokesBlockPreconditioner as the level smoother.
 *
 * The outer FAC operations (restriction, prolongation, and residual
 * computation) are reused from the existing box-relaxation operator. The
 * smoother itself is replaced by a level-local block application:
 *
 * 1. solve an approximate pressure Poisson problem,
 * 2. apply the commutator approximation to the Schur complement,
 * 3. smooth the velocity block with the existing velocity FAC strategy.
 */
class MultiphaseStaggeredStokesBlockFACOperator : public IBTK::FACPreconditionerStrategy
{
public:
    /*!
     * \brief Construct a block-smoother FAC strategy.
     *
     * \param object_name SAMRAI object name.
     * \param default_options_prefix PETSc options prefix base.
     * \param params Discretization and material parameters.
     * \param thn_manager Accessor for hierarchy volume-fraction data.
     */
    MultiphaseStaggeredStokesBlockFACOperator(const std::string& object_name,
                                              const std::string& default_options_prefix,
                                              const MultiphaseParameters& params,
                                              const std::unique_ptr<VolumeFractionDataManager>& thn_manager);

    /*!
     * \brief Destructor.
     */
    ~MultiphaseStaggeredStokesBlockFACOperator();

    const std::unique_ptr<VolumeFractionDataManager>& getVolumeFractionManager() const
    {
        return d_thn_manager;
    }

    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& un_bc_coefs,
                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& us_bc_coefs,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    void setToZero(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error, int level_num) override;
    void restrictResidual(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                          int dest_level_num) override;
    void prolongError(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                      int dest_level_num) override;
    void prolongErrorAndCorrect(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& source,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dest,
                                int dest_level_num) override;

    void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps) override;

    bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                            int coarsest_ln) override;

    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num) override;

    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs) override;

    void deallocateOperatorState() override;

    /*!
     * \brief Set the diagonal stabilization coefficients used in the Stokes
     * block.
     */
    void setCandDCoefficients(double C, double D);

    /*!
     * \brief Select the symmetric or asymmetric velocity block smoother.
     */
    void setUsingSymmetric(bool using_symmetric);

private:
    /*!
     * \brief Build the scalar pressure smoother used by the Schur
     * approximation.
     *
     * The pressure stage must understand coarse-fine interfaces, so this uses
     * a cell-centered Poisson FAC smoother rather than a bare patch-level
     * solver.
     */
    void initializePressureSmoother();

    /*!
     * \brief Initialize the velocity-block smoother used in the final block
     * solve stage.
     */
    void initializeVelocitySmoother(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                    const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs);

    /*!
     * \brief Compute the `thn^2 + ths^2` side coefficient used by the pressure
     * Schur approximation on the active hierarchy range.
     */
    void updateThnThsSqData();

    /*!
     * \brief Fill patch-boundary and coarse-fine ghost values for the
     * intermediate pressure-gradient data on one FAC level.
     *
     * These side-centered scratch fields are consumed by the momentum block,
     * so their coarse-fine ghost cells must be extended the same way as the
     * existing velocity FAC smoothers.
     */
    void fillVelocityScratchData(int level_num);

    SAMRAI::tbox::Pointer<IBTK::PoissonFACPreconditionerStrategy> d_pressure_smoother;
    SAMRAI::solv::PoissonSpecifications d_pressure_coefs;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_pressure_solver_db;

    const MultiphaseParameters& d_params;
    const std::unique_ptr<VolumeFractionDataManager>& d_thn_manager;
    SAMRAI::tbox::Pointer<MultiphaseStaggeredStokesBoxRelaxationFACOperator> d_delegate;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    SAMRAI::tbox::Pointer<MultiphaseStaggeredVelocityBlockFACOperator> d_velocity_block_smoother;
    SAMRAI::tbox::Pointer<MultiphaseStaggeredVelocityAsymmetricFACOperator> d_velocity_asymmetric_smoother;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_un_corr_var, d_us_corr_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_p_corr_var, d_p_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_un_rhs_var, d_us_rhs_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_un_scr_var, d_us_scr_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_fn_scr_var, d_fs_scr_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_thn_ths_sq_var;

    int d_un_corr_idx = IBTK::invalid_index, d_us_corr_idx = IBTK::invalid_index;
    int d_p_corr_idx = IBTK::invalid_index, d_p_rhs_idx = IBTK::invalid_index;
    int d_un_rhs_idx = IBTK::invalid_index, d_us_rhs_idx = IBTK::invalid_index;
    int d_un_scr_idx = IBTK::invalid_index, d_us_scr_idx = IBTK::invalid_index;
    int d_fn_scr_idx = IBTK::invalid_index, d_fs_scr_idx = IBTK::invalid_index;
    int d_thn_ths_sq_idx = IBTK::invalid_index;

    double d_C = std::numeric_limits<double>::quiet_NaN();
    double d_D = std::numeric_limits<double>::quiet_NaN();
    bool d_using_symmetric = true;

    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_un_bc_coefs, d_us_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef = nullptr;
};
} // namespace multiphase

#endif
