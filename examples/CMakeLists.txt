project(phelm_examples CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
include_directories(../src)
link_directories(../lib)

create_test_sourcelist(test_chafe_source test_chafe_driver.cpp test_chafe.cpp test_schafe.cpp)
create_test_sourcelist(test_laplace_source test_laplace_driver.cpp test_laplace.cpp test_slaplace.cpp)

set(chafe_source
	# Spheric Chafe-Infante equation
	schafe.cpp
	schafe.h
	schafe_impl.h
	
	# Flat Chafe-Infante equation
	chafe.cpp
	chafe.h
	chafe_impl.h)

add_executable(test_laplace ${test_laplace_source})
add_executable(test_chafe ${test_chafe_source} ${chafe_source})

add_executable(test_system_laplace test_system_laplace.cpp)
add_executable(test_jacobian test_jacobian.cpp)
add_executable(test_barvortex test_barvortex.cpp barvortex.h barvortex_impl.h)
add_executable(test_baroclin test_baroclin.cpp baroclin.cpp baroclin.h)

target_link_libraries(test_laplace phelm ${UMFPACK} ${SUPERLU})
target_link_libraries(test_chafe phelm ${UMFPACK} ${SUPERLU})

target_link_libraries(test_system_laplace phelm ${UMFPACK} ${SUPERLU})
target_link_libraries(test_jacobian phelm ${UMFPACK} ${SUPERLU})
target_link_libraries(test_barvortex phelm ${UMFPACK} ${SUPERLU})
target_link_libraries(test_baroclin phelm ${UMFPACK} ${SUPERLU})

if (OPENMP_FOUND)
set_target_properties(test_laplace PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_chafe PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

set_target_properties(test_system_laplace PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_jacobian PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_barvortex PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_baroclin PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (OPENMP_FOUND)

if (CUDA_FOUND)
	cuda_add_executable(test_laplace_cu ${test_laplace_source})
	target_link_libraries(test_laplace_cu phelm_cu ${UMFPACK})
	cuda_add_cublas_to_target(test_laplace_cu)
	set_target_properties(test_laplace_cu PROPERTIES COMPILE_FLAGS "-DGPGPU")

	cuda_add_executable(test_chafe_cu ${test_chafe_source} ${chafe_source})
	target_link_libraries(test_chafe_cu phelm_cu ${UMFPACK})
	cuda_add_cublas_to_target(test_chafe_cu)
	set_target_properties(test_chafe_cu PROPERTIES COMPILE_FLAGS "-DGPGPU")

	if (OPENMP_FOUND)
		set_target_properties(test_laplace_cu PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
		set_target_properties(test_chafe_cu PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	endif (OPENMP_FOUND)
endif (CUDA_FOUND)

