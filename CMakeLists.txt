cmake_minimum_required(VERSION 2.8)
project(blue_mapper)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(blue_mapper_build_tests "Build unit tests" OFF)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_library(blue_alignment STATIC src/blue_alignment.cpp)
add_library(blue_minimizers STATIC src/blue_minimizers.cpp)

add_executable(blue_mapper src/blue_mapper.cpp)

target_link_libraries(blue_mapper bioparser blue_alignment blue_minimizers)

if(blue_mapper_build_tests)
	include_directories(${PROJECT_SOURCE_DIR}/src)
	
	add_executable(blue_mapper_tests test/blue_mapper_tests.cpp)
	
	add_subdirectory(test/googletests EXCLUDE_FROM_ALL)
	
	target_link_libraries(blue_mapper_tests blue_alignment blue_minimizers gtest_main)
endif(blue_mapper_build_tests)
