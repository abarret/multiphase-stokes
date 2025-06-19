## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

## Note this file is mostly copied from IBAMR's repository.

ADD_CUSTOM_TARGET(tests)

# We use the list of test directories in two ways:
# 1. each test is added to a target tests-dir so that, e.g., 'make tests-IBFE'
#    only compiles the IBFE tests.
# 2. all input and output files in these source directories are symlinked into
#    their corresponding build directories.
SET(TEST_DIRECTORIES apply multigrid time_stepping time_stepping_thn frm convective_op physical_bc fd_operators)

FOREACH(_dir ${TEST_DIRECTORIES})
  ADD_CUSTOM_TARGET("tests-${_dir}")
  ADD_DEPENDENCIES(tests "tests-${_dir}")
ENDFOREACH()

# Convenience macro that sets up an executable target which links against
# _target_link. For example, if the inputs are Foo, bar.cpp, and quux then we
# create a target Foo_bar in directory Foo which depends on quux.
MACRO(SETUP _dir _src _target_link)
  GET_FILENAME_COMPONENT(_out_name "${_src}" NAME_WE)
  SET(_target "tests-${_dir}_${_out_name}")
  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL "${_dir}/${_src}")
  SET_TARGET_PROPERTIES(${_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/${_dir}"
    OUTPUT_NAME
    ${_out_name}
    )
  TARGET_COMPILE_OPTIONS(${_target} PUBLIC -DSOURCE_DIR="${CMAKE_SOURCE_DIR}/tests/${_dir}/")
  TARGET_LINK_LIBRARIES(${_target} PRIVATE "${_target_link}")
  ADD_DEPENDENCIES("tests-${_dir}" ${_target})
ENDMACRO()

# Like the last one, but for 2d targets. For example, if the inputs are Foo and
# bar.cpp then we create a target Foo_bar_2d in directory Foo which depends on
# IBAMR2d.
MACRO(SETUP_2D _dir _src)
  GET_FILENAME_COMPONENT(_dest "${_src}" NAME_WE)
  SET(_out_name "${_dest}_2d")
  SET(_target "tests-${_dir}_${_out_name}")
  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL "${_dir}/${_src}")
  SET_TARGET_PROPERTIES(${_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/${_dir}"
    OUTPUT_NAME
    ${_out_name}
    )
  TARGET_COMPILE_OPTIONS(${_target} PUBLIC -DSOURCE_DIR="${CMAKE_SOURCE_DIR}/tests/${_dir}/")
  TARGET_LINK_LIBRARIES(${_target} PRIVATE MultiphaseLib)
  ADD_DEPENDENCIES("tests-${_dir}" ${_target})
ENDMACRO()

# Like the last one, but for 3d targets. For example, if the inputs are Foo and
# bar.cpp then we create a target Foo_bar_3d in directory Foo which depends on
# IBAMR3d.
MACRO(SETUP_3D _dir _src)
  GET_FILENAME_COMPONENT(_dest "${_src}" NAME_WE)
  SET(_out_name "${_dest}_3d")
  SET(_target "tests-${_dir}_${_out_name}")
  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL "${_dir}/${_src}")
  SET_TARGET_PROPERTIES(${_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/${_dir}"
    OUTPUT_NAME
    ${_out_name}
    )
  TARGET_COMPILE_OPTIONS(${_target} PUBLIC -DSOURCE_DIR="${CMAKE_SOURCE_DIR}/tests/${_dir}/")
  TARGET_LINK_LIBRARIES(${_target} PRIVATE MultiphaseLib)
  ADD_DEPENDENCIES(tests-${_dir} ${_target})
ENDMACRO()

# Inform CMake of tests. This should have code such as SETUP_2D(<folder> source_code.cpp)

ADD_CUSTOM_COMMAND(TARGET tests
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/attest ${CMAKE_BINARY_DIR}/attest)

# Set up the input and output files. Since the input and output files aren't
# really used by the build system we use a shell script to find them every time
# 'make test' is run rather than evaluating the glob when cmake generates the
# build system.
FOREACH(_dir ${TEST_DIRECTORIES})
  ADD_CUSTOM_COMMAND(TARGET "tests-${_dir}"
    POST_BUILD
    COMMAND bash ${CMAKE_SOURCE_DIR}/tests/link-test-files.sh
    ${CMAKE_SOURCE_DIR}/tests/${_dir} ${CMAKE_BINARY_DIR}/tests/${_dir}
    VERBATIM)
ENDFOREACH()

SETUP_2D(apply apply.cpp)
SETUP_2D(apply applyStokesBlock.cpp)
SETUP_2D(multigrid multigrid.cpp)
SETUP_2D(multigrid multigridStokesBlock.cpp)
SETUP_2D(multigrid pressureSolve.cpp)
SETUP_2D(time_stepping time_stepping.cpp)
SETUP_2D(time_stepping_thn time_stepping_thn.cpp)
SETUP_2D(frm four_roll_mill.cpp)
SETUP_2D(convective_op convective_op.cpp)
SETUP_2D(physical_bc apply.cpp)
SETUP_2D(physical_bc channel.cpp)
SETUP_2D(physical_bc solve.cpp)
SETUP_2D(fd_operators fd_operator.cpp)

# Find numdiff, if possible (we only need it for tests so its not essential that
# we find it now)
FIND_PROGRAM(NUMDIFF_EXECUTABLE NAMES numdiff HINTS ${NUMDIFF_ROOT} PATH_SUFFIXES bin)

IF ("${NUMDIFF_EXECUTABLE}" STREQUAL "NUMDIFF_EXECUTABLE-NOTFOUND")
  MESSAGE(WARNING "\
The configuration script was not able to locate numdiff. If you want to run \
the test suite you will need to either edit attest.conf, specify the path to \
numdiff to attest, or rerun CMake with the argument NUMDIFF_ROOT specifying \
numdiff's root installation directory.")
  # clear the value so that attest.conf doesn't contain an invalid path
  SET(NUMDIFF_EXECUTABLE "")
ENDIF()

# Set up the default attest configuration file:
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/tests/attest.conf.in
  ${CMAKE_BINARY_DIR}/attest.conf)
