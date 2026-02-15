cmake_minimum_required(VERSION 3.18)
project(parabem LANGUAGES C CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_WITH_MPI "Enable OpenMP support" ON)

if(WIN32)
  add_compile_definitions(_DISABLE_EXTENDED_ALIGNED_STORAGE)
endif()

if(POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()
if(POLICY CMP0148)
  cmake_policy(SET CMP0148 NEW)
endif()

## python
if (BUILD_PYTHON_BINDINGS)
  set(PYBIND11_FINDPYTHON ON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 REQUIRED)
endif()

## eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

## mpi
if(BUILD_WITH_MPI)
  find_package(OpenMP REQUIRED)
endif(BUILD_WITH_MPI)

add_subdirectory(parabem)
add_subdirectory(src)
