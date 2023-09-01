# Set minimum CMake version
cmake_minimum_required(VERSION 3.15)

# Set project versions: major, minor and patch
set(TEST_VERSION_MAJOR "0")
set(TEST_VERSION_MINOR "1")
set(TEST_VERSION_PATCH "0")
set(TEST_VERSION "${TEST_VERSION_MAJOR}.${TEST_VERSION_MINOR}.${TEST_VERSION_PATCH}")

# Options
option(UNIT_TESTING "Enable unit tests" OFF)
option(USE_MPI "Compile using MPI" ON)

option(USE_PCPOWER "Compiling in PC-POWER" OFF)
option(USE_MN "Compiling in MareNostrum" OFF)
option(USE_CRAY "Compiling with CrayMpi" OFF)
option(SHARED_LIBS "Compile using shared libraries" OFF)

option(USE_GPU "Compile using GPU" OFF)
option(USE_MEM_MANAGED "Compile using Mem managed" OFF)

# Folder with files configuring extra CMake options
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set MPI Compilers if required
if(USE_CRAY)
    message("Compiling with CrayMpi compilers")
    set(ENV{CC} "cc")
    set(ENV{CXX} "CC")
    set(ENV{FC} "ftn")
else()    
    set(ENV{CC} "mpicc")
    set(ENV{CXX} "mpicxx")
    set(ENV{FC} "mpif90")
endif()

# Set project name, version and languages
project(App_TMU LANGUAGES C CXX Fortran)

# Configure MPI module
include(mpi)

# Configure compiler options
include(compilerOps)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

add_subdirectory(src)

# Enable testing
if(UNIT_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
