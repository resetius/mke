project(auxiliary CXX C)

include_directories(/usr/include/suitesparse/ ../src)
set(LIBRARY_OUTPUT_PATH "../lib")

add_library(phelm_aux 
	baroclin.cpp  
	barvortex.cpp  
	jacobian.cpp  
	laplace.cpp  
	sphere_laplace.cpp
	baroclin.h  
	barvortex.h  
	jacobian.h  
	laplace.h
	)


