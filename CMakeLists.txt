cmake_minimum_required(VERSION 2.6)
project(mgr)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
add_compile_options(-std=c++11)

find_package(Eigen3 REQUIRED)
# find_package(MATIO REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${MATIO_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/include")

add_library(mgr++
    SHARED
    src/SubsetSum.cpp
    src/Backward.cpp
    )

target_link_libraries(mgr++ dw)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

add_subdirectory(test)