cmake_minimum_required(VERSION 2.8.4)
project(paraBEM)

MESSAGE(STATUS "python version set: ${py}")

SET(CMAKE_BUILD_TYPE REALEASE)
SET(CMAKE_CXX_FLAGS_REALEASE "-O3")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
ADD_DEFINITIONS(-std=c++11)

# OPTIONS
OPTION(BUILD_PYTHON_BINDINGS "create python bindings" ON)


# FIND DEPENDENCIES

## EIGEN
find_package(Eigen3)
MESSAGE(${EIGEN3_VERSION})
find_package(OpenMP)
INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR})     # eigen3

## OPENMP
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
endif()


add_subdirectory(paraBEM)
add_subdirectory(src)

MESSAGE(STATUS "installing to ${python_libs}")
