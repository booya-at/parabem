# Build the c++ shared library
include_directories( ${CMAKE_SOURCE_DIR}/src/headers)
set(HEADER_FILES ${CMAKE_SOURCE_DIR}/src/headers/case2.h
				 ${CMAKE_SOURCE_DIR}/src/headers/case3.h
				 ${CMAKE_SOURCE_DIR}/src/headers/element_influence.h
				 ${CMAKE_SOURCE_DIR}/src/headers/lifting_line.h
				 ${CMAKE_SOURCE_DIR}/src/headers/panel2.h
				 ${CMAKE_SOURCE_DIR}/src/headers/panel3.h
				 ${CMAKE_SOURCE_DIR}/src/headers/vector.h)

set(CPP_FILES ${CMAKE_SOURCE_DIR}/src/element_influence.cpp
			  ${CMAKE_SOURCE_DIR}/src/panel3.cpp
			  ${CMAKE_SOURCE_DIR}/src/panel2.cpp
			  ${CMAKE_SOURCE_DIR}/src/case2.cpp
			  ${CMAKE_SOURCE_DIR}/src/case3.cpp
			  ${CMAKE_SOURCE_DIR}/src/lifting_line.cpp)



IF(BUILD_PYTHON_BINDINGS)
	pybind11_add_module(parabem_cpp python/parabem_ext.cpp ${CPP_FILES} ${HEADER_FILES})
	target_link_libraries(parabem_cpp PUBLIC Eigen3::Eigen)
	install(DIRECTORY "../parabem/" DESTINATION ${PYTHON_INSTALL_DIR}/parabem FILES_MATCHING PATTERN "*.py")
else()
	add_library(pm STATIC ${CPP_FILES} ${HEADER_FILES})
	target_link_libraries (pm Eigen3::Eigen)
	install(TARGETS pm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
ENDIF(BUILD_PYTHON_BINDINGS)
