cmake_minimum_required (VERSION 3.2)
project(orange_mapper)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(mapper_build_tests "Unit tests for mapper" OFF)

#The version number.
set (orange_mapper_VERSION_MAJOR 1)
set (orange_mapper_VERSION_MINOR 0)
set (orange_mapper_VERSION_PATCH 1)

set(CMAKE_CXX_STANDARD 11)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/OrangeConfig.h.in"
    "${PROJECT_BINARY_DIR}/OrangeConfig.h"
    )

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

add_library(orange_alignment ${PROJECT_SOURCE_DIR}/src/orange_alignment.cpp)
add_library(orange_minimizers ${PROJECT_SOURCE_DIR}/src/orange_minimizers.cpp)

if (mapper_build_tests)
    include_directories(${PROJECT_SOURCE_DIR}/src)
    add_executable(orange_mapper_unittest test/orange_mapper_tests.cpp)
    add_subdirectory(vendor/googletest EXCLUDE_FROM_ALL)
    target_link_libraries(orange_mapper_unittest orange_alignment orange_minimizers gtest_main)
endif (mapper_build_tests)

add_executable(orange_mapper src/orange_mapper.cpp)
target_link_libraries(orange_mapper bioparser orange_alignment orange_minimizers thread_pool pthread)
