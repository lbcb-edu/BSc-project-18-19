cmake_minimum_required (VERSION 2.6)
project (white_mapper)

set (white_mapper_VERSION_MAJOR 1)
set (white_mapper_VERSION_MINOR 0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

configure_file(
	"${PROJECT_SOURCE_DIR}/src/white_mapper.h.in"
	"${PROJECT_BINARY_DIR}/white_mapper.h"
)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory (vendor/bioparser EXCLUDE_FROM_ALL)

add_executable (white_mapper src/white_mapper.cxx)
target_link_libraries (white_mapper bioparser)
