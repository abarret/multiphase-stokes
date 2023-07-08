CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

ENABLE_LANGUAGE(Fortran)
PROJECT(MultiphaseStokes LANGUAGES C CXX Fortran)
FIND_PACKAGE(IBAMR REQUIRED)

SET(MULTIPHASE_SRC
    src/FullFACPreconditioner.cpp
    src/INSVCTwoFluidStaggeredHierarchyIntegrator.cpp
    src/VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp
    src/VCTwoFluidStaggeredStokesOperator.cpp
    src/RBGS.f
    )

ADD_LIBRARY(MultiphaseLib)
TARGET_SOURCES(MultiphaseLib PRIVATE ${MULTIPHASE_SRC})
TARGET_INCLUDE_DIRECTORIES(MultiphaseLib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/include/>)
TARGET_LINK_LIBRARIES(MultiphaseLib PUBLIC IBAMR::IBAMR2d)

ADD_EXECUTABLE(apply2d apply.cpp)
TARGET_LINK_LIBRARIES(apply2d MultiphaseLib)

ADD_EXECUTABLE(solver2d solver.cpp)
TARGET_LINK_LIBRARIES(solver2d MultiphaseLib)

ADD_EXECUTABLE(smooth2d smooth.cpp)
TARGET_LINK_LIBRARIES(smooth2d MultiphaseLib)

ADD_EXECUTABLE(fluid_multigrid multigrid.cpp)
TARGET_LINK_LIBRARIES(fluid_multigrid MultiphaseLib)

ADD_EXECUTABLE(time_stepping time_stepping.cpp)
TARGET_LINK_LIBRARIES(time_stepping MultiphaseLib)

ADD_EXECUTABLE(time_stepping_thn time_stepping_thn.cpp)
TARGET_LINK_LIBRARIES(time_stepping_thn MultiphaseLib)

ADD_EXECUTABLE(four_roll_mill four_roll_mill.cpp)
TARGET_LINK_LIBRARIES(four_roll_mill MultiphaseLib)
