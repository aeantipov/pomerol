cmake_minimum_required(VERSION 3.1.0)

project(pomerol CXX)
set(POMEROL_VERSION 2.0)
message(STATUS "Configuring ${PROJECT_NAME} ${POMEROL_VERSION}")

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Disable in-source builds
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# RPATH fix
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Library linking type
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Compiler quirks
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR
   "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
        message(STATUS "Adding -stdlib=libc++ flag")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()

# Enable/disable and find OpenMP
option(POMEROL_USE_OPENMP "Use OpenMP" TRUE)
if(POMEROL_USE_OPENMP)
    find_package(OpenMP)
else()
    message(STATUS "OpenMP disabled")
endif(POMEROL_USE_OPENMP)

#
# Dependencies
#

# Eigen 3
find_package (Eigen3 3.1 REQUIRED)
message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR})

# libcommute
find_package(libcommute 0.7.0 CONFIG REQUIRED)
message(STATUS "Found libcommute version ${libcommute_VERSION}")

# MPI
find_package (MPI 3.0 REQUIRED)
message(STATUS "MPI includes: ${MPI_CXX_INCLUDE_PATH}")
message(STATUS "MPI CXX libs: ${MPI_CXX_LIBRARIES}")
message(STATUS "MPI flags: ${MPI_CXX_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS}")

# Boost
find_package (Boost 1.54.0 REQUIRED)
message(STATUS "Boost includes: ${Boost_INCLUDE_DIRS}" )

# Main library
add_subdirectory(include)
add_subdirectory(src)

# Build executables
option(Progs "Build executables" OFF)
if(Progs)
    add_subdirectory(prog)
endif(Progs)

# Enable unit tests
option(Testing "Enable testing" ON)
if (Testing)
    enable_testing()
    add_subdirectory(test)
endif(Testing)

# Build documentation
option(Documentation "Build documentation" ON)
if(Documentation)
    add_subdirectory(doc)
endif()

# Install pkg-config file
configure_file(pomerol.pc.in pomerol.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerol.pc
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

# Configure cmake config
configure_file(pomerolConfig.cmake.in pomerolConfig.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerolConfig.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pomerol)

# Configure lmod file
configure_file(pomerol.lmod.in pomerol.lmod @ONLY)
