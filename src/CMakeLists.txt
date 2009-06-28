project(core CXX C)

include_directories(/usr/include/suitesparse/)
set(LIBRARY_OUTPUT_PATH "../lib")

add_library(phelm 
	phelm.cpp 
	util.cpp 
	polynom.cpp 
	solver.cpp 
	ver.cpp
	gmres.c 
	quadrature.c
	phelm.h
	phelm_generators.h
	phelm_private.h
	gmres.h
	polynom.h
	solver.h
	util.h
	)

