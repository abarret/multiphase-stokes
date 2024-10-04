# Adaptive Mesh Refinement for Multiphase Fluids

This repository contains the source code to simulate multiphase flows on adaptively refined grids. The solver implemented here is second order accuracte in $L^1$, $L^2$, $L^\infty$ norms for all solution variables of the Newtonian and non-Newtonian models. Our experiments demonstrate the solver is stable provided $\Delta t$ satisfies the imposed CFL condition. The solver can accurately resolve sharp gradients in the solution and, with the multigrid preconditioner, its behavior is independent of grid spacing. This AMR solver offers a major cost savings benefits by using refined grids in areas of interest and using coarser grids elsewhere.

# Dependencies

You will need to [build and install IBAMR](https://ibamr.github.io/linux) and all its required dependencies. We recommend compiling the optimized builds of the dependencies if you simply want to try out some simulations. Additionally, you should compile IBAMR with CMake since we provide a CMakeLists.txt to compile the source code in this repository. 

The simulations generated can be viewed by installing [VISit](https://visit-dav.github.io/visit-website/) on your machine.

# Compiling the code

1. In the project’s root directory, create a build directory so that all object files and executables are stored separately from the source code.

```
mkdir build
```

2. Enter the build directory and run cmake to configure and generate build files for our `multiphaseLib` library, and link it to IBAMR. The command below uses the `O3` compiler optimization flag to enable level 3 optimizations. 


```
cd build
cmake /path/to/CMakeLists.txt \
    -DIBAMR_ROOT=/path/to/ibamr/install/dir \
    -DCMAKE_CXX_FLAGS="-O3"
```
Note: If cmake fails, then make sure you're using the same compilers that were used for compiling IBAMR. You may provide a path to your compilers by passing flags for `-DCMAKE_CXX_COMPILER`, `-DCMAKE_C_COMPILER` and `-DCMAKE_Fortran_COMPILER`.

3. Now you can build and compile the source code using `make`. Targets can be found by viewing the CMakeLists.txt file in a particular subdirectory, e.g., in `/frm` there are two targets, `frm` and `frm_convergence`. You can compile the target `frm` as follows

```
make frm
```

4. Finally, you can run the executable produced by make. For example, you can execute the `frm` executable by giving the input filename and other command-line arguments as follows

```
cd frm
./frm input2d \
    -solver_ksp_initial_guess_nonzero FALSE \
    -solver_ksp_monitor_true_residual \
    -solver_ksp_type fgmres \
    -solver_ksp_max_it 20 \
    -solver_ksp_rtol 1e-14 \
```

5. While the code is being executed, you may view the log files to monitor progress. After the run is complete, you may use VISit to plot the simulations results.

# Additional Support

If you run into issues compiling or running this code, you may contact nagdabindi@gmail.com or barrett@math.utah.edu for additional support. 

# Acknowledgements

We are grateful to receive support for this work from NIH grants 1R01GM131408, U01HL143336, and R01HL157631 as well as NSF OAC 1931516.
