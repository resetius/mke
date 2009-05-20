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

if (MSVC)
	add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

add_definitions(${OpenMP_C_FLAGS})

if (CYGWIN)
add_definitions(-DGMRES) # comment it if you want umfpack
endif (CYGWIN)
add_definitions(-DSPARSE)

add_subdirectory(contrib)

if (MSVC)
	if (NOT UMFPACK_SOURCE_FOUND)
		add_definitions(-DGMRES)
	endif (NOT UMFPACK_SOURCE_FOUND)
endif (MSVC)

add_subdirectory(src)
add_subdirectory(auxiliary)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(examples)

