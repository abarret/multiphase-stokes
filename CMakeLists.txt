CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

PROJECT(MultiphaseStokes LANGUAGES C CXX Fortran)

INCLUDE(GNUInstallDirs)
INCLUDE(CMakePackageConfigHelpers)

FIND_PACKAGE(IBAMR REQUIRED)

SET(MULTIPHASE_SRC
    src/FullFACPreconditioner.cpp
    src/IBMultiphaseCrossLinks.cpp
    src/IBMultiphaseEulerianForceFunction.cpp
    src/IBMultiphaseHierarchyIntegrator.cpp
    src/INSVCTwoFluidStaggeredHierarchyIntegrator.cpp
    src/ScaleStress.cpp
    src/StressRelaxation.cpp
    src/VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp
    src/VCTwoFluidStaggeredStokesOperator.cpp
    src/RBGS.f
    )

ADD_LIBRARY(MultiphaseLib)
TARGET_SOURCES(MultiphaseLib PRIVATE ${MULTIPHASE_SRC})
TARGET_INCLUDE_DIRECTORIES(MultiphaseLib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/include/>)
TARGET_LINK_LIBRARIES(MultiphaseLib PUBLIC IBAMR::IBAMR2d)

INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/multiphase DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
INSTALL(TARGETS MultiphaseLib EXPORT MultiphaseTargets)
INSTALL(EXPORT "MultiphaseTargets" FILE "MultiphaseTargets.cmake" NAMESPACE Multiphase:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/multiphase)

WRITE_BASIC_PACKAGE_VERSION_FILE(
	"${PROJECT_BINARY_DIR}/cmake/MultiphaseConfigVersion.cmake"
	VERSION 0.0.0
	COMPATIBILITY AnyNewerVersion)

SET(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/multiphase)
EXPORT(EXPORT "MultiphaseTargets" FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/MultiphaseTargets.cmake NAMESPACE Multiphase::)
CONFIGURE_PACKAGE_CONFIG_FILE(
	${CMAKE_SOURCE_DIR}/cmake/MultiphaseConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/cmake/MultiphaseConfig.cmake
	INSTALL_DESTINATION ${INSTALL_CONFIGDIR})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/MultiphaseConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/MultiphaseConfigVersion.cmake
	DESTINATION ${INSTALL_CONFIGDIR})

SET(DIRS apply smooth multigrid time_stepping time_stepping_thn frm ib_ex ib_peristalsis)
FOREACH(_dir ${DIRS})
   ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
