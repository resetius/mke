project(auxiliary CXX C)

include_directories(/usr/include/suitesparse/ ../src)
set(LIBRARY_OUTPUT_PATH "../lib")

set(SOURCE 
	baroclin.cpp  
	barvortex.cpp  
	jacobian.cpp  
	laplace.cpp  
	sphere_laplace.cpp
	baroclin.h  
	barvortex.h  
	jacobian.h  
	laplace.h)

add_library(phelm_aux 
	${SOURCE}
	)

add_library(phelm_aux_cu
	${SOURCE}
	)

set_target_properties(phelm_aux_cu PROPERTIES COMPILE_FLAGS "-DGPGPU -I${CUDA_SDK}/include")

