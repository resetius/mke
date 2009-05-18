cmake_minimum_required(VERSION 2.6.3)
project(mke CXX C)

include(CTest)
enable_testing()
IF(BUILD_TESTING)
   SET(BUILDNAME "${BUILDNAME}" CACHE STRING "MKE")
   MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(FindOpenMP)

if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

if (NOT WIN32)
if (NOT CYGWIN)
	set(ICC ".*icpc.*")
	if (CMAKE_CXX_COMPILER MATCHES ${ICC})
		set(ICC_LIBRARIES "./lib-icc") # CHANGE IT!
		set(ICC_STOCK_FORTRAN_LIBRARIES "/opt/intel/Compiler/11.0/081/lib/intel64/") #CHANGE IT!
		set(FLIB "ifcoremt")
		link_directories(${ICC_LIBRARIES} ${ICC_STOCK_FORTRAN_LIBRARIES})
		include_directories(/usr/include/suitesparse/)
		set(UMFPACK umfpack amd blas ${FLIB})
		#add_definitions(-DGMRES)
	else (CMAKE_CXX_COMPILER MATCHES ${ICC})
		set(FLIB "gfortran")
		include_directories(/usr/include/suitesparse/)
		set(UMFPACK umfpack amd)
	endif (CMAKE_CXX_COMPILER MATCHES ${ICC})
endif (NOT CYGWIN)
set(MATH m)
else (NOT WIN32)
set(MATH "")
endif (NOT WIN32)

if (MSVC)
	add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

add_definitions(${OpenMP_C_FLAGS})
if (CYGWIN)
add_definitions(-DGMRES) # comment it if you want umfpack
endif (CYGWIN)
add_definitions(-DSPARSE)

add_subdirectory(contrib)
message("${UMFPACK_SOURCE_FOUND}")

if (MSVC)

set(UMFPACK_SOURCE_FOUND 1)

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/UMFPACK")
	set(UMFPACK_SOURCE_FOUND 0)
	message("${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/UMFPACK not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/UMFPACK")

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack//AMD")
	set(UMFPACK_SOURCE_FOUND 0)
	message("${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack//AMD not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack//AMD")

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack//UFconfig")
	set(UMFPACK_SOURCE_FOUND 0)
	message("${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack//UFconfig not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack//UFconfig")

if (UMFPACK_SOURCE_FOUND)
include_directories(./contrib/umfpack/UMFPACK/Include)
include_directories(./contrib/umfpack/AMD/Include)
include_directories(./contrib/umfpack/UFconfig)
set(UMFPACK umfpack)
else (UMFPACK_SOURCE_FOUND)
add_definitions(-DGMRES)
endif (UMFPACK_SOURCE_FOUND)

endif (MSVC)

add_subdirectory(src)
add_subdirectory(auxiliary)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(examples)


