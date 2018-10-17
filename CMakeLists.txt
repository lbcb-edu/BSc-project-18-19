cmake_minimum_required(VERSION 3.2)
project(brown_mapper)

set(brown_mapper_VERSION_MAJOR 1)
set(brown_mapper_VERSION_MINOR 0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

configure_file(
  "${PROJECT_SOURCE_DIR}/brown_mapper.hpp.in"
  "${PROJECT_BINARY_DIR}/brown_mapper.hpp"
)

configure_file(
  "${PROJECT_SOURCE_DIR}/bioparser/include/bioparser/bioparser.hpp"
  "${PROJECT_BINARY_DIR}/bioparser.hpp"
  COPY_ONLY
)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(bioparser EXCLUDE_FROM_ALL)

add_executable(brown_mapper brown_mapper.cc)
target_link_libraries(brown_mapper bioparser)

#install(TARGETS brown_mapper DESTINATION bin)
#install(FILES "${PROJECT_BINARY_DIR}/brown_mapper.hpp" DESTINATION include)
