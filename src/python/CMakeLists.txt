EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c
"from distutils.sysconfig import get_python_lib; import os; print(os.path.join(get_python_lib(), 'paraBEM'))"
OUTPUT_VARIABLE _PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )

string(REPLACE "\\" "/" PYTHON_INSTALL_DIR ${_PYTHON_INSTALL_DIR})


MESSAGE(STATUS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
MESSAGE(STATUS "installing to ${python_libs}")

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS})    # python.h

# Build the python shared library
# PYTHON_ADD_MODULE(_paraBEM SHARED paraBEM_ext.cpp)
pybind11_add_module(_paraBEM SHARED paraBEM_ext.cpp)
set_target_properties(_paraBEM PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/paraBEM)

set_target_properties(_paraBEM PROPERTIES PREFIX "")
# set_target_properties(_paraBEM PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
if( WIN32 )
    set_target_properties(_paraBEM PROPERTIES SUFFIX ".pyd")
endif()

target_link_libraries(_paraBEM PUBLIC pm)

install(DIRECTORY "../../paraBEM/" DESTINATION ${PYTHON_INSTALL_DIR} FILES_MATCHING PATTERN "*.py")
install(TARGETS _paraBEM DESTINATION ${PYTHON_INSTALL_DIR})
