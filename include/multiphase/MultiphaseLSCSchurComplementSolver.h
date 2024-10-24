#ifndef included_MultiphaseLSCSchurComplementSolver
#define included_MultiphaseLSCSchurComplementSolver

#include <ibtk/GeneralSolver.h>

namespace multiphase
{
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
