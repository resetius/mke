project(core CXX C)

include_directories(/usr/include/suitesparse/)
set(LIBRARY_OUTPUT_PATH "../lib")

set(CUDA_SDK "/opt/cuda")

set(SOURCE 
	phelm.cpp 
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
	util.h)


add_library(phelm 
	${SOURCE}
	util.cpp 
	)

if (CUDA_SDK)
add_library(phelm_cu
	${SOURCE}
	util_cu.cpp 
	)

set_target_properties(phelm_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_SDK}/include")
endif (CUDA_SDK)


