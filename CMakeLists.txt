cmake_minimum_required(VERSION 2.6.3)
project(mke CXX C)

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
	add_definitions(-DGMRES)
endif (MSVC)

add_definitions(${OpenMP_C_FLAGS})
if (CYGWIN)
add_definitions(-DGMRES) # comment it if you want umfpack
endif (CYGWIN)
add_definitions(-DSPARSE)

add_subdirectory(src)
add_subdirectory(utils)
add_subdirectory(tests)


