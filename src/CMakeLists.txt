project(core CXX C)

include_directories(${LINAL_INCLUDE})
include_directories(.)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
add_definitions(${LINAL_DEFINES})

set(SOURCE 
	phelm.cpp 
	polynom.cpp 
	ver.cpp
	util.cpp
	quadrature.c
	phelm.h
	generators.h
	generators.cpp
	impl/generators_impl.h
	polynom.h
	util.h
	norm.h
	mesh.cpp
	mesh.h
	func.cpp
	func.h
	point.cpp
	point.h
	
	deriv.h
	impl/deriv_impl.h

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
target_compile_features(phelm PUBLIC cxx_range_for)

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

