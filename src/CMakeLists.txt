# Build the c++ shared library
add_library(pm SHARED element_influence.cpp panel3.cpp panel2.cpp case2.cpp case3.cpp  lifting_line.cpp)

# Build the python shared library
add_library(_paraBEM SHARED paraBEM_ext.cpp)
target_link_libraries(_paraBEM pm)
set_target_properties(_paraBEM PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/paraBEM)
set_target_properties(_paraBEM PROPERTIES PREFIX "")

# add main
add_executable(main main.cpp)
target_link_libraries(main pm)

install(TARGETS _paraBEM DESTINATION ${python_install_dir})
install(TARGETS pm DESTINATION lib)
