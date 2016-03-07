cmake_minimum_required(VERSION 2.8.4)
project(pm)

MESSAGE(STATUS "python version set: ${py}")

SET(CMAKE_BUILD_TYPE REALEASE)
SET(CMAKE_CXX_FLAGS_REALEASE "-O3")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
ADD_DEFINITIONS(-std=c++11)

#find necessary packages

set(Python_ADDITIONAL_VERSIONS ${py})
find_package(PythonInterp ${py})
find_package(PythonLibs ${py})
find_package(PyEigen)

IF (NOT boost)
  find_package( Boost COMPONENTS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
ELSE()
  find_package( Boost COMPONENTS ${boost})
ENDIF()

message($PYTHON_EXECUTABLE)
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c 
  "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )

MESSAGE(STATUS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")

MESSAGE(STATUS "installing to ${python_libs}")

find_package(Eigen3)
MESSAGE(${EIGEN3_VERSION})
find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
endif()

SET(python_install_dir ${python_libs}/ppm)

#include drectories
include_directories( ${Boost_INCLUDE_DIR})
include_directories( ${PYTHON_INCLUDE_DIRS})
include_directories( ${EIGEN3_INCLUDE_DIR})


add_subdirectory(ppm)
add_subdirectory(src)

install(DIRECTORY ppm DESTINATION ${python_libs} FILES_MATCHING
    PATTERN "*.py")

MESSAGE(STATUS "installing to ${python_libs}")
