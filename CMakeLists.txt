cmake_minimum_required (VERSION 3.5)

project (lbcb-mapper)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++1y")
set (source_dir "${PROJECT_SOURCE_DIR}/src/")

option(pink_mapper_build_tests "Build pink mapper unit tests" OFF)

file (GLOB source_files "${source_dir}/*.cpp")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_library(pink_alignment src/pink_alignment.cpp)

include_directories("${PROJECT_BINARY_DIR}/lib")

if (pink_mapper_build_tests)
   include_directories(${PROJECT_SOURCE_DIR}/src)

   add_executable(pink_mapper_test test/pink_mapper_test.cpp)

   add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)

   target_link_libraries(pink_mapper_test pink_alignment gtest_main)
else ()
   add_executable (pink_mapper ${source_files})

   target_link_libraries(pink_mapper bioparser pink_alignment)
endif (pink_mapper_build_tests)
   

