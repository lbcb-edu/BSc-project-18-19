cmake_minimum_required (VERSION 2.6)
project (white_mapper)

set (white_mapper_VERSION_MAJOR 1)
set (white_mapper_VERSION_MINOR 0)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

option (white_mapper_build_tests "Build white_mapper unit tests" ON)

configure_file (
	"${PROJECT_SOURCE_DIR}/src/white_mapper.h.in"
	"${PROJECT_BINARY_DIR}/bin/white_mapper.h"
)

add_library(white_alignment STATIC src/white_alignment.cxx)
add_library(white_minimizers STATIC src/white_minimizers.cxx)

include_directories("${PROJECT_BINARY_DIR}/bin")

add_subdirectory (vendor/bioparser EXCLUDE_FROM_ALL)

add_executable (white_mapper src/white_mapper.cxx)
target_link_libraries (white_mapper bioparser white_alignment white_minimizers)

if (white_mapper_build_tests)
	include_directories(${PROJECT_SOURCE_DIR}/src)

	add_executable (white_mapper_test test/white_mapper_test.cxx)

	add_subdirectory (vendor/googletest/googletest EXCLUDE_FROM_ALL)

	target_link_libraries (white_mapper_test white_alignment white_minimizers gtest_main)
endif (white_mapper_build_tests)
