cmake_minimum_required (VERSION 3.2)
project(orange_mapper)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

#The version number.
set (orange_mapper_VERSION_MAJOR 0)
set (orange_mapper_VERSION_MINOR 1)
set (orange_mapper_VERSION_PATCH 0)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/OrangeConfig.h.in"
    "${PROJECT_BINARY_DIR}/OrangeConfig.h"
    )


include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
include_directories("${PROJECT_SOURCE_DIR}/vendor/bioparser/include/bioparser")

add_subdirectory(src)