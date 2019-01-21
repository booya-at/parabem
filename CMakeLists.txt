cmake_minimum_required(VERSION 3.0)
project(paraBEM)
include_directories("/Users/simon/.nix-profile/include")
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

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(OpenMP_C "${CMAKE_C_COMPILER}")
    set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
    set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
  set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
endif()
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

add_subdirectory(paraBEM)
add_subdirectory(src)

MESSAGE(STATUS "installing to ${python_libs}")
