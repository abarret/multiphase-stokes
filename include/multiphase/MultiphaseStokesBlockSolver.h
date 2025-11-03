#ifndef included_MultiphaseStokesBlockSolver
#define included_MultiphaseStokesBlockSolver

#include <ibtk/GeneralSolver.h>

namespace multiphase
{
class MultiphaseStokesBlockSolver : IBTK::GeneralSolver
{
public:
    MultiphaseStokesBlockSolver(std::string object_name);

    virtual ~MultiphaseStokesBlockSolver();

    /*!
     * Solve Ax = b with A the multiphase stokes block.
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
