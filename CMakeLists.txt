cmake_minimum_required(VERSION 3.0)
project(parabem)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # needed?

option(BUILD_PYTHON_BINDINGS "build the python bindings" ON)
option(BUILD_WITH_MPI "switch multithreading support" ON)
option(BUILD_PYTHON_BINDINGS "create python bindings" ON)

if (WIN32)
  add_definitions(-D_DISABLE_EXTENDED_ALIGNED_STORAGE)
endif()


## python
if (BUILD_PYTHON_BINDINGS)
  ## PYTHON
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  find_package(pybind11 REQUIRED)
endif()

## eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


## mpi
if(BUILD_WITH_MPI)
  find_package(OpenMP REQUIRED)
  ## OPENMP
  if(OPENMP_FOUND)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
          link_libraries(${OpenMP_libomp_LIBRARY})
          link_libraries(${OpenMP_libomp_LIBRARY})
      endif()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
  endif(OPENMP_FOUND)
endif(BUILD_WITH_MPI)

add_subdirectory(parabem)
add_subdirectory(src)
