project(mke CXX C)

cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

if (NOT WIN32)
if (NOT CYGWIN)
	set(ICC ".*icpc.*")
	if (CMAKE_CXX_COMPILER MATCHES ${ICC})
		set(OMP "-openmp")
		set(ICC_LIBRARIES "./lib-icc") # CHANGE IT!
		set(ICC_STOCK_FORTRAN_LIBRARIES "/opt/intel/Compiler/11.0/081/lib/intel64/") #CHANGE IT!
		set(FLIB "ifcoremt")
		link_directories(${ICC_LIBRARIES} ${ICC_STOCK_FORTRAN_LIBRARIES})
		set(UMFPACK umfpack amd blas ${FLIB})
	else (CMAKE_CXX_COMPILER MATCHES ${ICC})
		set(OMP "-fopenmp")
		set(FLIB "gfortran")
		#set(UMFPACK umfpack amd)
	endif (CMAKE_CXX_COMPILER MATCHES ${ICC})
endif (NOT CYGWIN)
set(MATH m)
else (NOT WIN32)
set(MATH "")
endif (NOT WIN32)

add_definitions(${OMP})
add_definitions(-DGMRES) # comment it if you want umfpack
add_definitions(-DSPARSE)

include_directories(/usr/include/suitesparse/)

add_library(mke mke.cpp util.cpp polynom.cpp laplace.cpp 
	solver.cpp gmres.c quadrature.c sphere_laplace.cpp)
add_executable(test_laplace test_laplace.cpp)
add_executable(test_slaplace test_slaplace.cpp)
add_executable(sphere icosahedron.cpp)
add_executable(rectangle rectangle.cpp)
add_executable(test_polynom test_polynom.cpp)
add_executable(test_integrate test_integrate.cpp)
add_executable(test_gmres test_gmres.cpp)
add_executable(test_quadrature test_quadrature.c)
add_executable(test_funcgen test_funcgen.cpp)

target_link_libraries(test_laplace mke ${UMFPACK})
target_link_libraries(test_slaplace mke ${UMFPACK})
target_link_libraries(test_polynom mke ${UMFPACK})
target_link_libraries(test_integrate mke ${UMFPACK})
target_link_libraries(test_gmres mke)
target_link_libraries(test_quadrature mke ${MATH})
target_link_libraries(test_funcgen mke ${MATH})

if (DEFINED OMP)
set_target_properties(test_laplace PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_slaplace PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_polynom PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_integrate PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(rectangle PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(sphere PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_gmres PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_quadrature PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_funcgen PROPERTIES LINK_FLAGS ${OMP})
endif (DEFINED OMP)


