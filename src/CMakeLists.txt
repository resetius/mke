project(core CXX C)

include_directories(/usr/include/suitesparse/)
set(LIBRARY_OUTPUT_PATH "../lib")

set(SOURCE 
	phelm.cpp 
	polynom.cpp 
	ver.cpp
	util.cpp
	linal_base.cpp
	quadrature.c
	phelm.h
	base.h
	linal.h
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
	proj.cpp
	)

if (CUDA_FOUND)
cuda_add_library(phelm_cu
	${SOURCE}
	linal_cuda.cpp 
	linal_cuda.h
	proj_cuda.cpp
	alloc_cu.h
	shmem.h
	texture.h
	linal.cu
	phelm.cu
	)
set_target_properties(phelm_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_TOOLKIT_ROOT_DIR}/include")
endif (CUDA_FOUND)


