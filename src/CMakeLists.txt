# set version number
set (linearCSD_VERSION_MAJOR 1)
set (linearCSD_VERSION_MINOR 0)
set (linearCSD_VERSION_PATCH 0)

message (STATUS "Cmake version ${CMAKE_VERSION}")
message (STATUS "Configuring linearCSD ${linearCSD_VERSION_MAJOR}.${linearCSD_VERSION_MINOR}.${linearCSD_VERSION_PATCH}")
message (STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}")
message (STATUS "Include path: ${CMAKE_INCLUDE_PATH}")
message (STATUS "Library path: ${CMAKE_LIBRARY_PATH}")
message (STATUS "Prefix path: ${CMAKE_PREFIX_PATH}")
message (STATUS "Module path: ${CMAKE_MODULE_PATH}")
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message (STATUS "linearFSI library directory: ${LINEARFSI_LIB_DIR}")

if (LINEARFSI_LIB_DIR)
  set (CMAKE_PREFIX_PATH ${LINEARFSI_LIB_DIR})
else ()
  message(WARNING "LINEARFSI_LIB_DIR is not set, this is probably a mistake")
endif()

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# add include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${LINEARFSI_LIB_DIR}/linearFSI)
include_directories(${LINEARFSI_LIB_DIR}/libUG_lite)
include_directories(${LINEARFSI_LIB_DIR}/libFSI)
include_directories(${PROJECT_BINARY_DIR})

# configure header file to pass cmake variables
configure_file (
  "${CMAKE_SOURCE_DIR}/include/macros.hpp.in"
  "${PROJECT_BINARY_DIR}/macros.hpp"
  )

# set all source files
set(sources
  inputs.cpp
  main.cpp
  utility.cpp
  vec3.cpp
  )

# create targets for executable
add_executable (linearCSD ${sources})

# set c++ standard for all targets
set_property (TARGET linearCSD PROPERTY CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
message (STATUS "Using C++17")

# add compiler specific flags
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(STATUS "Changing stdlib for Clang")
  set(COMPILER_SPEC_FLAGS "-stdlib=libc++")
endif()

# set additional c++ flags for all build types
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_SPEC_FLAGS} -Wall -march=native -O3 -DNDEBUG -DMPI_NO_CPPBIND -DOMPI_SKIP_MPICXX")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${COMPILER_SPEC_FLAGS} -Wall -O0 -ggdb -pg -coverage")
set (CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -pg")

# set compiler flags based on build type
if (CMAKE_BUILD_TYPE MATCHES "release")
  message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_RELEASE})
elseif (CMAKE_BUILD_TYPE MATCHES "debug")
  message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_DEBUG})
elseif (CMAKE_BUILD_TYPE MATCHES "profile")
  message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_PROFILE})
else ()
  message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS})
endif ()

# add mapper libraries
target_link_libraries (linearCSD -L${LINEARFSI_LIB_DIR}/linearFSI)
target_link_libraries (linearCSD -lCSE)
target_link_libraries (linearCSD -lgfortran)
# needed for gcc < 9
target_link_libraries (linearCSD stdc++fs)

# determine if MPI_DIR has been set
if (MPI_DIR)
  message (STATUS "Looking for MPI here: ${MPI_DIR}")
  set (CMAKE_PREFIX_PATH ${MPI_DIR})
endif()

# get MPI and link to all targets
if (CRAY)
  message (STATUS "Compiling on Cray system, not using find_package for MPI")
  set (MPI_INCLUDE_PATH ${MPI_DIR}/include)
  message (STATUS "Using MPI include path: ${MPI_INCLUDE_PATH}")
  include_directories(SYSTEM "${MPI_INCLUDE_PATH}")
else ()
  find_package(MPI REQUIRED)
  include_directories(SYSTEM "${MPI_INCLUDE_PATH}")
  target_link_libraries (linearCSD ${MPI_C_LIBRARIES})
endif ()

# determine if PETSC_DIR has been set
if (PETSc_DIR)
  message (STATUS "Looking for PETSc here: ${PETSc_DIR}")
  set (CMAKE_PREFIX_PATH ${PETSC_DIR})
endif()

# get PETSc and link to all targets
find_package(PETSc REQUIRED)
include_directories(SYSTEM "${PETSc_INCLUDE_DIRS}")
target_link_libraries (linearCSD ${PETSc_LIBRARIES})

# install executable
install (TARGETS linearCSD RUNTIME DESTINATION bin)
