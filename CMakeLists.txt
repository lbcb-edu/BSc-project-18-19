cmake_minimum_required(VERSION 3.2)
project(brown_mapper)

set(brown_mapper_VERSION_MAJOR 2)
set(brown_mapper_VERSION_MINOR 0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(mapper_build_tests "Build mapper unit tests" OFF)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/brown_mapper.hpp.in"
  "${PROJECT_BINARY_DIR}/brown_mapper.hpp"
)

add_library(brown_alignment STATIC src/brown_alignment.cc)
add_library(brown_minimizers STATIC src/brown_minimizers.cc)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

if (mapper_build_tests)
  include_directories(${PROJECT_SOURCE_DIR}/src)

  add_executable(brown_mapper_unittest test/brown_mapper_unittest.cc)

  add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)

  target_link_libraries(brown_mapper_unittest brown_alignment gtest_main brown_minimizers)
endif (mapper_build_tests)
add_executable(brown_mapper src/brown_mapper.cc)
target_link_libraries(brown_mapper bioparser thread_pool pthread brown_alignment brown_minimizers)
