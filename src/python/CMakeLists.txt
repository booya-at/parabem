execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
    "import sysconfig; print(sysconfig.get_paths()['purelib'])"
    OUTPUT_VARIABLE _PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )

string(REPLACE "\\" "/" PYTHON_INSTALL_DIR ${_PYTHON_INSTALL_DIR})


message(STATUS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
message(STATUS "installing to ${python_libs}")

include_directories( ${PYTHON_INCLUDE_DIRS})    # python.h

pybind11_add_module(parabem_cpp SHARED parabem_ext.cpp)

if( WIN32 )
    set_target_properties(parabem_cpp PROPERTIES SUFFIX ".pyd")
endif()

target_link_libraries(parabem_cpp PUBLIC pm)

install(DIRECTORY "../../parabem/" DESTINATION ${PYTHON_INSTALL_DIR}/parabem FILES_MATCHING PATTERN "*.py")
install(TARGETS parabem_cpp DESTINATION ${PYTHON_INSTALL_DIR})
