cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(svm-linear-classifier)

add_executable(svm-linear-classifier svm.cpp)
target_link_libraries(svm-linear-classifier visualization ${GMP_LIBRARIES})
