## PYTHON
set(Python_ADDITIONAL_VERSIONS ${py})
find_package(PythonInterp ${py})
find_package(PythonLibs ${py})

EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c 
"from distutils.sysconfig import get_python_lib; print(get_python_lib())"
OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )
SET(PYTHON_INSTALL_DIR ${python_libs}/paraBEM)

MESSAGE(STATUS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
MESSAGE(STATUS "installing to ${python_libs}")

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS})    # python.h

# Build the python shared library
add_library(_paraBEM SHARED paraBEM_ext.cpp)
target_link_libraries(_paraBEM pm)
set_target_properties(_paraBEM PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/paraBEM)
set_target_properties(_paraBEM PROPERTIES PREFIX "")
set_target_properties(_paraBEM PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(_paraBEM ${PYTHON_LIBRARIES})

install(DIRECTORY "../../paraBEM/" DESTINATION ${PYTHON_INSTALL_DIR} FILES_MATCHING PATTERN "*.py")
install(TARGETS _paraBEM DESTINATION ${PYTHON_INSTALL_DIR})