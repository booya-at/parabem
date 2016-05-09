# Build the c++ shared library
include_directories( ${CMAKE_SOURCE_DIR}/src/headers)
set(HEADER_FILES ${CMAKE_SOURCE_DIR}/src/headers/case2.h
				 ${CMAKE_SOURCE_DIR}/src/headers/case3.h
				 ${CMAKE_SOURCE_DIR}/src/headers/element_influence.h
				 ${CMAKE_SOURCE_DIR}/src/headers/lifting_line.h
				 ${CMAKE_SOURCE_DIR}/src/headers/panel2.h
				 ${CMAKE_SOURCE_DIR}/src/headers/panel3.h
				 ${CMAKE_SOURCE_DIR}/src/headers/vector.h)

set(CPP_FILES element_influence.cpp
			  ${CMAKE_SOURCE_DIR}/src/panel3.cpp
			  ${CMAKE_SOURCE_DIR}/src/panel2.cpp
			  ${CMAKE_SOURCE_DIR}/src/case2.cpp
			  ${CMAKE_SOURCE_DIR}/src/case3.cpp 
			  ${CMAKE_SOURCE_DIR}/src/lifting_line.cpp)

add_library(pm SHARED ${CPP_FILES} ${HEADER_FILES})

# Build the python shared library
add_library(_paraBEM SHARED paraBEM_ext.cpp)
target_link_libraries(_paraBEM pm)
set_target_properties(_paraBEM PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/paraBEM)
set_target_properties(_paraBEM PROPERTIES PREFIX "")

install(TARGETS _paraBEM DESTINATION ${python_install_dir})
install(TARGETS pm DESTINATION lib)
