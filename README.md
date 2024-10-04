# Adaptive Mesh Refinement for Multiphase Fluids

This repository contains the source code to simulate multiphase flows on adaptively refined grids. The solver implemented here is second order accuracte in $L^1$, $L^2$, $L^\infty$ norms for all solution variables of the Newtonian and non-Newtonian models. Our experiments demonstrate the solver is stable provided $\Delta t$ satisfies the imposed CFL condition. The solver can accurately resolve sharp gradients in the solution and uses an efficient multigrid preconditioner to solve the resulting saddle point system. This AMR solver offers substantial cost savings by using refined grids in areas of interest and using coarser grids elsewhere. 

For additional details on the numerical solver and simulation results, please refer to our paper: [Adaptive Mesh Refinement for Two-Phase Viscoelastic Fluid Mixture Models](https://arxiv.org/abs/2409.19974).

# Dependencies

You will need to [build and install IBAMR](https://ibamr.github.io/building) and all its required dependencies. This library uses CMake, and IBAMR must also be compiled with CMake. Unless you intend to actively develop the code, we suggest compiling optimized versions of the dependencies.

Note that IBAMR and its dependencies can be built automatically with `autoibamr`, and we are compatible with this method of building: https://github.com/ibamr/autoibamr.

The simulations generated can be viewed by installing [VISit](https://visit-dav.github.io/visit-website/) on your machine.

# Compiling the code

1. In the project’s root directory, create a build directory so that all object files and executables are stored separately from the source code.

```
mkdir build
```

2. In the build directory, run `cmake` to configure and generate makefiles for the library and examples. For optimized builds, we suggest using the compilation flags `-O3 -march=native`. Note that the flag `-march=native` enables hardware-specific flags that are unsuitable for cross-platform builds.


```
cd build
cmake .. \
    -DIBAMR_ROOT=/path/to/ibamr/install/dir \
    -DCMAKE_CXX_FLAGS="-O3 -march=native"
```

3. Now you can build and compile the source code using `make`. Targets can be found by viewing the CMakeLists.txt file in a particular subdirectory, e.g., in `/frm` there are two targets, `frm` and `frm_convergence`. You can compile the target `frm` as follows

```
make frm
```

4. Finally, you can run the executable produced by `make`. For example, you can execute the `frm` executable by giving the input filename and other command-line arguments as follows

```
cd frm
./frm input2d \
    -solver_ksp_initial_guess_nonzero FALSE \
    -solver_ksp_monitor_true_residual \
    -solver_ksp_type fgmres \
    -solver_ksp_max_it 20 \
    -solver_ksp_rtol 1e-14 \
```

5. While the code is being executed, you may view the log files to monitor progress. After the run is complete, you may use VisIt to plot the simulation results.

# Additional Support

If you run into issues compiling or running this code, you may open a disussion on GitHub for additional support.

# Acknowledgements

We are grateful to receive support for this work from NIH grants 1R01GM131408, U01HL143336, and R01HL157631 as well as NSF OAC 1931516.
