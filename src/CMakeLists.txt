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

add_library(pm STATIC ${CPP_FILES} ${HEADER_FILES})

install(TARGETS pm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)


IF(BUILD_PYTHON_BINDINGS)
	add_subdirectory(python)
ENDIF(BUILD_PYTHON_BINDINGS)
