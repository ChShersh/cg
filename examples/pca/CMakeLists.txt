cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(pca)

add_executable(pca pca.cpp)
target_link_libraries(pca visualization ${GMP_LIBRARIES})
