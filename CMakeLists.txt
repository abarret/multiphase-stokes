CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

PROJECT(MultiphaseStokes LANGUAGES C CXX Fortran)
FIND_PACKAGE(IBAMR REQUIRED)

SET(CXX_SRC
    VCTwoFluidStaggeredStokesOperator.cpp
    VCTwoFluidStaggeredStokesBoxRelaxationFACOperator.cpp
    )

SET(CMAKE_CXX_COMPILER
   /export/spack/linux-centos7-x86_64/gcc-7.5.0/mpich-3.4.2-wj5qtlnmxknr7oycz5wblj3iwsbpgc3m/bin/mpicxx
   )

ADD_EXECUTABLE(main2d main.cpp)
ADD_EXECUTABLE(solver2d solver.cpp)
TARGET_SOURCES(main2d PRIVATE ${CXX_SRC})
TARGET_LINK_LIBRARIES(main2d IBAMR::IBAMR2d)
TARGET_SOURCES(solver2d PRIVATE ${CXX_SRC})
TARGET_LINK_LIBRARIES(solver2d IBAMR::IBAMR2d)
CONFIGURE_FILE(input2d_four ${CMAKE_CURRENT_BINARY_DIR}/input2d_four COPYONLY)
