cmake_minimum_required(VERSION 2.8.0)
project(phelm CXX C)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Common)
include(CTest)
include(Testing)
find_package(OpenMP)
find_package(Math)
find_package(Fortran)
find_package(CUDA)
find_package(BISON)
find_package(NETCDF)
include(CommonDefinitions)

if (EXISTS "${CMAKE_CURRENT_BINARY_DIR}/settings.cmake")
	include(${CMAKE_CURRENT_BINARY_DIR}/settings.cmake)
endif()

add_subdirectory(src)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(examples)

