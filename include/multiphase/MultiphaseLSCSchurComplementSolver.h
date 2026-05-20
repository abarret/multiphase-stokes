#ifndef included_MultiphaseLSCSchurComplementSolver
#define included_MultiphaseLSCSchurComplementSolver

#include <ibtk/GeneralSolver.h>

namespace multiphase
{
/*!
 * \brief Class MultiphaseLSCSchurComplementSolver defines the GeneralSolver interface for solving the LSC Schur
 * complement system associated with the multiphase Stokes block formulation.
 *
 * The linear system is
 * S*x = b, with S = (G^T*G)*(G^T*F*G)^-1*(G^T*G),
 * and concrete implementations solve it approximately as the pressure stage of the block preconditioner.
 */
class MultiphaseLSCSchurComplementSolver : IBTK::GeneralSolver
{
public:
    MultiphaseLSCSchurComplementSolver(std::string object_name);

    virtual ~MultiphaseLSCSchurComplementSolver();

    /*!
     * Solve G^T*G(G^T*F*G)^-1(G^T*G)x = b
     */
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * Initialize the solver state
     */
    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * Deallocate all temporary data used to solve the linear system.
     */
    void deallocateSolverState() override;
};
} // namespace multiphase
#endif
