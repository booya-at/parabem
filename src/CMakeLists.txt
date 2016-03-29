# Build our library
add_library( pm SHARED element_influence.cpp panel3.cpp panel2.cpp case3.cpp case2.cpp  lifting_line.cpp)
set_target_properties(pm PROPERTIES PREFIX "")

PYTHON_ADD_MODULE(_ppm SHARED ppm_ext.cpp)

# target_link_libraries(_ppm ${PY_EIGEN_FILE_PATH} pm)
set_target_properties(_ppm PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
set_target_properties(_ppm PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/ppm)
set_target_properties(_ppm PROPERTIES PREFIX "")

add_executable(main main.cpp)
target_link_libraries(main pm)


install(TARGETS _ppm DESTINATION ${python_install_dir})
install(TARGETS pm DESTINATION /usr/lib/)
