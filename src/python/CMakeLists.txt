if(DEFINED Python_SITEARCH AND NOT Python_SITEARCH STREQUAL "")
  set(PYTHON_INSTALL_DIR "${Python_SITEARCH}")
else()
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_paths()['purelib'])"
    OUTPUT_VARIABLE _PYTHON_INSTALL_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(PYTHON_INSTALL_DIR "${_PYTHON_INSTALL_DIR}")
endif()
string(REPLACE "\\" "/" PYTHON_INSTALL_DIR "${PYTHON_INSTALL_DIR}")

message(STATUS "python-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
message(STATUS "installing to ${PYTHON_INSTALL_DIR}")

# Avoid ambiguous build-tree rpath auto-generation with mixed system/env runtime libs.
set(CMAKE_SKIP_BUILD_RPATH TRUE)
pybind11_add_module(parabem_cpp MODULE parabem_ext.cpp)

if(WIN32)
  set_target_properties(parabem_cpp PROPERTIES SUFFIX ".pyd")
endif()

target_link_libraries(parabem_cpp PRIVATE pm)

install(DIRECTORY "../../parabem/" DESTINATION "${PYTHON_INSTALL_DIR}/parabem" FILES_MATCHING PATTERN "*.py")
install(TARGETS parabem_cpp DESTINATION "${PYTHON_INSTALL_DIR}")
