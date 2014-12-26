cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(line-arrangement-example)

add_executable(line-arrangement-example line_arrangement.cpp)
target_link_libraries(line-arrangement-example visualization ${GMP_LIBRARIES})
