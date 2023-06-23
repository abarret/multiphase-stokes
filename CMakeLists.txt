CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

ENABLE_LANGUAGE(Fortran)
PROJECT(MultiphaseStokes LANGUAGES C CXX Fortran)
FIND_PACKAGE(IBAMR REQUIRED)

SET(MULTIPHASE_SRC
    src/CFMultiphaseOldroydB.cpp
    src/FullFACPreconditioner.cpp
    src/IBMultiphaseEulerianForceFunction.cpp
    src/IBMultiphaseHierarchyIntegrator.cpp
    src/INSVCTwoFluidConvectiveManager.cpp
    src/INSVCTwoFluidStaggeredHierarchyIntegrator.cpp
    src/VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp
    src/VCTwoFluidStaggeredStokesOperator.cpp
    src/RBGS.f
    )

ADD_LIBRARY(MultiphaseLib)
TARGET_SOURCES(MultiphaseLib PRIVATE ${MULTIPHASE_SRC})
TARGET_INCLUDE_DIRECTORIES(MultiphaseLib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/include/>)
TARGET_LINK_LIBRARIES(MultiphaseLib PUBLIC IBAMR::IBAMR2d)

SET(DIRS apply smooth multigrid time_stepping time_stepping_thn frm convective_op ib_ex)
FOREACH(_dir ${DIRS})
   ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
