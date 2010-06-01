project(core CXX C)

include_directories(${LINAL_INCLUDE})
include_directories(.)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

set(SOURCE 
	phelm.cpp 
	polynom.cpp 
	ver.cpp
	util.cpp
	quadrature.c
	phelm.h
	generators.h
	impl/generators_impl.h
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
	gpu/shmem.h
	gpu/texture.h
	gpu/phelm.cu
	)
set_target_properties(phelm_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_TOOLKIT_ROOT_DIR}/include")
endif (CUDA_FOUND)

