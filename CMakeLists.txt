project(mke CXX C)

cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

if (NOT WIN32)
if (NOT CYGWIN)
	set(ICC ".*icpc.*")
	if (CMAKE_CXX_COMPILER MATCHES ${ICC})
		set(OMP_CFLAGS "-openmp")
		set(OMP_LFLAGS "-openmp")
		set(OMP_LIBS   "")
		set(ICC_LIBRARIES "./lib-icc") # CHANGE IT!
		set(ICC_STOCK_FORTRAN_LIBRARIES "/opt/intel/Compiler/11.0/081/lib/intel64/") #CHANGE IT!
		set(FLIB "ifcoremt")
		link_directories(${ICC_LIBRARIES} ${ICC_STOCK_FORTRAN_LIBRARIES})
		set(UMFPACK umfpack amd blas ${FLIB})
	else (CMAKE_CXX_COMPILER MATCHES ${ICC})
		set(OMP_CFLAGS "-fopenmp")
		set(OMP_LFLAGS "-fopenmp")
		set(OMP_LIBS   "")
		set(FLIB "gfortran")
		#set(UMFPACK umfpack amd)
	endif (CMAKE_CXX_COMPILER MATCHES ${ICC})
endif (NOT CYGWIN)
set(MATH m)
else (NOT WIN32)
set(MATH "")
endif (NOT WIN32)

if (MSVC)
	#set(OMP_CFLAGS "/openmp")
	add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

add_definitions(${OMP_CFLAGS})
add_definitions(-DGMRES) # comment it if you want umfpack
add_definitions(-DSPARSE)

add_subdirectory(src)
add_subdirectory(utils)
add_subdirectory(tests)


