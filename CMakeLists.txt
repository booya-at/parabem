cmake_minimum_required(VERSION 3.0)
project(paraBEM)

MESSAGE(STATUS "python version set: ${py}")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# OPTIONS
set (CMAKE_CXX_STANDARD 11)
OPTION(BUILD_PYTHON_BINDINGS "create python bindings" ON)

if (BUILD_PYTHON_BINDINGS)
  ## PYTHON
  set(Python_ADDITIONAL_VERSIONS ${py})
  find_package(PythonInterp ${py} REQUIRED)
  find_package(PythonLibs ${py} REQUIRED)
  find_package(pybind11 2.0 REQUIRED)
endif()

# FIND DEPENDENCIES

## EIGEN
find_package(Eigen3 REQUIRED)
MESSAGE(${EIGEN3_VERSION})
INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR})     # eigen3

find_package(OpenMP)


## OPENMP
if(OPENMP_FOUND)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        LINK_LIBRARIES(${OpenMP_libomp_LIBRARY})
        LINK_LIBRARIES(${OpenMP_libomp_LIBRARY})
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
endif()

add_subdirectory(paraBEM)
add_subdirectory(src)

MESSAGE(STATUS "installing to ${python_libs}")
