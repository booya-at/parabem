cmake_minimum_required(VERSION 3.0)
project(parabem)
MESSAGE(STATUS "python version set: ${py}")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

OPTION(BUILD_PYTHON_BINDINGS "build the python bindings" ON)
OPTION(BUILD_WITH_MPI "switch multithreading support" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

if (WIN32)
  add_definitions(-D_DISABLE_EXTENDED_ALIGNED_STORAGE)
endif()

# OPTIONS
set (CMAKE_CXX_STANDARD 11)
OPTION(BUILD_PYTHON_BINDINGS "create python bindings" ON)


if (BUILD_PYTHON_BINDINGS)
  ## PYTHON
  
  set(Python_ADDITIONAL_VERSIONS ${py})
  find_package(PythonInterp ${py} REQUIRED)
  find_package(PythonLibs ${py} REQUIRED)

  if(EXISTS ${PROJECT_SOURCE_DIR}/lib/pybind11)
    add_subdirectory(${PROJECT_SOURCE_DIR}/lib/pybind11)
  else()
    find_package(pybind11 REQUIRED)
  endif()
  
  #set(PYBIND11_CPP_STANDARD -std=c++14)
endif()

# FIND DEPENDENCIES

## EIGEN
if(EXISTS ${PROJECT_SOURCE_DIR}/lib/eigen)
  add_subdirectory(${PROJECT_SOURCE_DIR}/lib/eigen)
else()
  find_package (Eigen3 3.3 REQUIRED NO_MODULE)
endif()

#INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR}/eigen3)     # eigen3

if (BUILD_WITH_MPI)
  find_package(OpenMP REQUIRED)
  ## OPENMP
  if(OPENMP_FOUND)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
          LINK_LIBRARIES(${OpenMP_libomp_LIBRARY})
          LINK_LIBRARIES(${OpenMP_libomp_LIBRARY})
      endif()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
  endif()
endif()

add_subdirectory(parabem)
add_subdirectory(src)

MESSAGE(STATUS "installing to ${python_libs}")
