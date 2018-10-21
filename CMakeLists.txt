cmake_minimum_required(VERSION 2.8)
project(blue_mapper)

set (CMAKE_CXX_STANDARD 11)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_executable(blue_mapper blue_mapper.cpp)

target_link_libraries(blue_mapper bioparser)
