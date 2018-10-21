cmake_minimum_required (VERSION 3.5)

project (lbcb-mapper)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++14")
set (source_dir "${PROJECT_SOURCE_DIR}/src/")

file (GLOB source_files "${source_dir}/*.cpp")

add_executable (pink_mapper ${source_files})

add_subdirectory(bioparser EXCLUDE_FROM_ALL)
target_link_libraries(pink_mapper bioparser)
