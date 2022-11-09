CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

ENABLE_LANGUAGE(Fortran)
PROJECT(MultiphaseStokes LANGUAGES C CXX Fortran)
FIND_PACKAGE(IBAMR REQUIRED)

SET(CXX_SRC
    VCTwoFluidStaggeredStokesOperator.cpp
    VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp
    )

SET(F_SRC
    RBGS.f
    )

ADD_EXECUTABLE(apply2d apply.cpp)
TARGET_SOURCES(apply2d PRIVATE VCTwoFluidStaggeredStokesOperator.cpp)
TARGET_LINK_LIBRARIES(apply2d IBAMR::IBAMR2d)

ADD_EXECUTABLE(solver2d solver.cpp)
TARGET_SOURCES(solver2d PRIVATE VCTwoFluidStaggeredStokesOperator.cpp)
TARGET_LINK_LIBRARIES(solver2d IBAMR::IBAMR2d)

ADD_EXECUTABLE(smooth2d smooth.cpp)
TARGET_SOURCES(smooth2d PRIVATE VCTwoFluidStaggeredStokesOperator.cpp VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp RBGS.f)
TARGET_LINK_LIBRARIES(smooth2d IBAMR::IBAMR2d)

ADD_EXECUTABLE(fluid_multigrid multigrid.cpp)
TARGET_SOURCES(fluid_multigrid PRIVATE VCTwoFluidStaggeredStokesOperator.cpp VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp RBGS.f)
TARGET_LINK_LIBRARIES(fluid_multigrid IBAMR::IBAMR2d)

ADD_SUBDIRECTORY(tests)