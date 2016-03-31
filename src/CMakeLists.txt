# Build the c++ shared library
add_library(pm SHARED element_influence.cpp panel3.cpp panel2.cpp case2.cpp case3.cpp  lifting_line.cpp)

# Build the python shared library
add_library(_ppm SHARED ppm_ext.cpp)
target_link_libraries(_ppm pm)
set_target_properties(_ppm PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/ppm)
set_target_properties(_ppm PROPERTIES PREFIX "")

# add main
add_executable(main main.cpp)
target_link_libraries(main pm)

install(TARGETS _ppm DESTINATION ${python_install_dir})
install(TARGETS pm DESTINATION /usr/lib/)
