project(core CXX C)

include_directories(/usr/include/suitesparse/)
set(LIBRARY_OUTPUT_PATH "../lib")

set(SOURCE 
	phelm.cpp 
	polynom.cpp 
	ver.cpp
	util.cpp
	gmres.c 
	quadrature.c
	phelm.h
	generators.h
	private.h
	solver_private.h
	gmres.h
	polynom.h
	solver.h
	util.h)


add_library(phelm 
	${SOURCE}
	linal.cpp 
	)

if (CUDA_SDK)
add_library(phelm_cu
	${SOURCE}
	linal_cuda.cpp 
	)

set_target_properties(phelm_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_SDK}/include")
endif (CUDA_SDK)


