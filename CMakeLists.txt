cmake_minimum_required(VERSION 2.8)
project(blue_mapper)

add_subdirectory(/home/sara/bioinf/blue_mapper/bioparser EXCLUDE_FROM_ALL)

add_executable(blue_mapper main.cpp)

target_link_libraries(blue_mapper bioparser)
