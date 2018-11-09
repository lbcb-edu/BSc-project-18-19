cmake_minimum_required(VERSION 3.2)
project(brown_mapper)

set(brown_mapper_VERSION_MAJOR 1)
set(brown_mapper_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/brown_mapper.hpp.in"
  "${PROJECT_BINARY_DIR}/brown_mapper.hpp"
)

add_library(brown_alignment STATIC src/brown_alignment.cc)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_executable(brown_mapper src/brown_mapper.cc)
target_link_libraries(brown_mapper bioparser brown_alignment)

#install(TARGETS brown_mapper DESTINATION bin)
#install(FILES "${PROJECT_BINARY_DIR}/brown_mapper.hpp" DESTINATION include)
