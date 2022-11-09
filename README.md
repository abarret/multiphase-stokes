# multiphase-stokes
Contains an operator class for variable-coefficient two-phase staggered stokes. See *VCTwoFluidStaggeredStokesOperator.cpp*.

There are two driver programs:
+ **apply.cpp**: Used to test that the operator is correct. This code applies the operator to user-defined network velocity, solvent velocity and pressure. The result is checked against the analytical solution to ensure the operator is 2nd order accurate.
+ **solver.cpp**: This solves the coupled two-phase system using a krylov method (GMRES). The krylov method requires the operator to be passed to it. 
+ **smooth.cpp**: This uses the stand-alone box relaxation scheme with gauss-seidel type iterations to solve the fluid equations

## FAC Precondioner Stragegy: VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp

Contains routines for smoothening the error using box relaxation, computing residual, prolongation, restriction and solving on coarsest level.
