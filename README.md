# Staggered Stokes Operator Test

The purpose of this repository is to test the application of the Staggered Stokes Operator when physical boundary conditions are employed. The physical boundary conditions employed here correspond to that of 2-D plane Poiseulle flow as can be seen in the input file: input2D. The MAC scheme is used to discretize the velocities and pressure. 

In the main.cpp file, the Staggered Stokes Operator is called and applied to the exact solution vector of 2-D plane Poiseulle flow (u_vec in the main.cpp file). The resulting vector (f_vec in the main.cpp file) is compared against what we'd expect to get after applying the Staggered Stokes Operator (e_vec in the main.cpp file) and the L1, L2, and L_infty errors are calculated. The resulting error is then interpolated onto the cell centers so that one may visualize the results in Visit. 

An acccompanying MATLAB script is provided which computes the Staggered Stokes Operator and interpolates the results onto the cell centers so that the user compare the results in Visit directly with the results obtained using MATLAB. 

In order to generate the makefile and executable (main2d) associated with this test, we use Cmake, and so a CMakeLists.txt file is provided. To generate the makefile simply run: cmake -DIBAMR_ROOT=<path to IBAMR install> . Depending on the users setup, one may need to include their paths to their compilers and their compiler flags. 

Once the makefile is generated, one can then run make to generate the executable main2d. 

