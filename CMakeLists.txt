CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

ENABLE_LANGUAGE(Fortran)
PROJECT(MultiphaseStokes LANGUAGES C CXX Fortran)
FIND_PACKAGE(IBAMR REQUIRED)

SET(MULTIPHASE_SRC
    src/CFMultiphaseOldroydB.cpp
    src/FullFACPreconditioner.cpp
    src/MultiphaseConvectiveManager.cpp
    src/MultiphaseStaggeredHierarchyIntegrator.cpp
    src/MultiphaseStaggeredStokesBoxRelaxationFACOperator.cpp
    src/MultiphaseStaggeredStokesOperator.cpp
    src/MultiphaseStaggeredStokesBlockFACOperator.cpp
    src/MultiphaseStaggeredStokesBlockOperator.cpp
    src/MultiphaseStaggeredStokesBlockPreconditioner.cpp
    src/RBGS.f
    src/VRBGS.f
    include/multiphase/private/FDOPS.f
    )

ADD_LIBRARY(MultiphaseLib)
TARGET_SOURCES(MultiphaseLib PRIVATE ${MULTIPHASE_SRC})
TARGET_INCLUDE_DIRECTORIES(MultiphaseLib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/include/>)
TARGET_LINK_LIBRARIES(MultiphaseLib PUBLIC IBAMR::IBAMR2d)

SET(DIRS smooth frm physical_bc convective_op)
FOREACH(_dir ${DIRS})
   ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()

ADD_SUBDIRECTORY(tests)
