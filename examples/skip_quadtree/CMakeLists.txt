cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(skip-quadtree-example)

add_executable(skip-quadtree-example skip_quadtree.cpp)
target_link_libraries(skip-quadtree-example visualization ${GMP_LIBRARIES})
