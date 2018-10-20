cmake_minimum_required (VERSION 3.2)
project(orange_mapper)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#The version number.
set (orange_mapper_VERSION_MAJOR 0)
set (orange_mapper_VERSION_MINOR 1)
set (orange_mapper_VERSION_PATCH 0)

set(CMAKE_CXX_STANDARD 11)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/OrangeConfig.h.in"
    "${PROJECT_BINARY_DIR}/OrangeConfig.h"
    )


include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_executable(orange_mapper ${PROJECT_SOURCE_DIR}/src/orange_mapper.cpp)
target_link_libraries(orange_mapper bioparser)
