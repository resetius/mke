project(auxiliary CXX C)

include_directories(/usr/include/suitesparse/ ../src)
set(LIBRARY_OUTPUT_PATH "../lib")

add_library(phelm_aux 
	baroclin.cpp  
	barvortex.cpp  
	jacobian.cpp
	sjacobian.cpp  
	laplace.cpp  
	sphere_laplace.cpp
	baroclin.h  
	barvortex.h  
	sjacobian.h
	jacobian.h	
	laplace.h
	)


