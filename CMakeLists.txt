cmake_minimum_required(VERSION 2.8)
project(blue_mapper)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

include_directories("${PROJECT_BINARY_DIR}")

add_executable(blue_mapper src/blue_mapper.cpp)

target_link_libraries(blue_mapper bioparser)
