project(core CXX C)

include_directories(/usr/include/suitesparse/)
include_directories(.)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

set(SOURCE 
	phelm.cpp 
	polynom.cpp 
	ver.cpp
	util.cpp
	quadrature.c
	phelm.h
	base.h
	generators.h
	impl/generators_impl.h
	gmres.h
	gmres.cpp
	polynom.h
	util.h
	norm.h
	
	laplace.cpp
	laplace.h
	impl/laplace_impl.h

	slaplace.cpp
	slaplace.h
	impl/slaplace_impl.h

	jacobian.cpp
	jacobian.h

	sjacobian.cpp
	sjacobian.h
	)


add_library(phelm 
	${SOURCE}
	proj.cpp
	)

if (CUDA_FOUND)
cuda_add_library(phelm_cu
	${SOURCE}
	gpu/proj_cuda.cpp
	gpu/alloc_cu.h
	gpu/shmem.h
	gpu/texture.h
	gpu/phelm.cu
	)
set_target_properties(phelm_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_TOOLKIT_ROOT_DIR}/include")
endif (CUDA_FOUND)
