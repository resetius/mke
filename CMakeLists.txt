cmake_minimum_required(VERSION 2.6.3)
project(mke CXX C)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(Common)
include(CTest)
include(Testing)
include(FindOpenMP)
include(FindMath)
include(FindFortran)
include(FindSelfUmfPack)
include(FindUmfPack)
include(CommonDefinitions)

add_subdirectory(contrib)
add_subdirectory(src)
add_subdirectory(auxiliary)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(examples)

