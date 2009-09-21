project(auxiliary CXX C)

include_directories(/usr/include/suitesparse/ ../src)
set(LIBRARY_OUTPUT_PATH "../lib")

set(SOURCE
	norm.h        # Norm and dot-product helpers

	baroclin.cpp  # The two-dimensional baroclinic atmosphere equations  on a sphere
	baroclin.h  

	impl/barvortex_impl.h # The Barotropic vorticity equation
	barvortex.h  

	sjacobian.cpp # Spheric Jacobian
	sjacobian.h

	jacobian.cpp  # Flat Jacobian
	jacobian.h

	laplace.cpp   # Flat Laplace equations

	laplace.h
	impl/laplace_impl.h

	slaplace.cpp  # Spheric Laplace equations 
	slaplace.h
	impl/slaplace_impl.h

	chafe.cpp     # Flat Chafe-Infante equation
	chafe.h
	impl/chafe_impl.h

	schafe.cpp    # Spheric Chafe-Infante equation
	schafe.h
	impl/schafe_impl.h
	)

add_library(phelm_aux 
	${SOURCE}
	)

if (CUDA_FOUND)
cuda_add_library(phelm_aux_cu
	${SOURCE}
	)

set_target_properties(phelm_aux_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_TOOLKIT_ROOT_DIR}/include")
endif (CUDA_FOUND)

