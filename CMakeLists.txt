cmake_minimum_required (VERSION 2.6)
project (Mapper)

add_subdirectory(bioparser EXCLUDE_FROM_ALL)

add_executable (white_mapper white_mapper.cxx)
target_link_libraries (white_mapper bioparser)
