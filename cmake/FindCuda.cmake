
find_path(CUDA_SDK "include/cuda.h" 
	HINTS /opt/cuda /cuda /usr/local/cuda c:/CUDA)

message(STATUS "CUDA SDK: ${CUDA_SDK}")
if (CUDA_SDK)
message(STATUS "CUDA SDK Found, building phelm_cu library and examples")
endif (CUDA_SDK)

