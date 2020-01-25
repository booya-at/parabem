EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c
"from distutils.sysconfig import get_python_lib; import os; print(os.path.join(get_python_lib(), 'parabem'))"
OUTPUT_VARIABLE _PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )

string(REPLACE "\\" "/" PYTHON_INSTALL_DIR ${_PYTHON_INSTALL_DIR})


MESSAGE(STATUS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
MESSAGE(STATUS "installing to ${python_libs}")

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS})    # python.h

# Build the python shared library
# PYTHON_ADD_MODULE(_parabem SHARED parabem_ext.cpp)
pybind11_add_module(_parabem SHARED parabem_ext.cpp)
set_target_properties(_parabem PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/parabem)

set_target_properties(_parabem PROPERTIES PREFIX "")
# set_target_properties(_parabem PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
if( WIN32 )
    set_target_properties(_parabem PROPERTIES SUFFIX ".pyd")
endif()

target_link_libraries(_parabem PUBLIC pm)

install(DIRECTORY "../../parabem/" DESTINATION ${PYTHON_INSTALL_DIR} FILES_MATCHING PATTERN "*.py")
install(TARGETS _parabem DESTINATION ${PYTHON_INSTALL_DIR})
