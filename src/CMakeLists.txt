project(core CXX C)

include_directories(/usr/include/suitesparse/)
include_directories(.)

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
	impl/generators_impl.h
	impl/solver_impl.h
	impl/solver_umfpack.h
	gmres.h
	gmres.cpp
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
	gpu/linal_cuda.cpp 
	gpu/linal_cuda.h
	gpu/proj_cuda.cpp
	gpu/alloc_cu.h
	gpu/shmem.h
	gpu/texture.h
	gpu/linal.cu
	gpu/phelm.cu
	gpu/reduction.h
	gpu/gmres.cu
	)
set_target_properties(phelm_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_TOOLKIT_ROOT_DIR}/include")
endif (CUDA_FOUND)


