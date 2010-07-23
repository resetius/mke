project(phelm_examples CXX C)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

include_directories(../src)
include_directories(${LINAL_INCLUDE})

add_definitions(${LINAL_DEFINES})

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

add_executable(mke_test_laplace ${test_laplace_source})
add_executable(mke_test_chafe ${test_chafe_source} ${chafe_source})

add_executable(mke_test_system_laplace test_system_laplace.cpp)
add_executable(mke_test_jacobian test_jacobian.cpp)
add_executable(mke_test_barvortex test_barvortex.cpp barvortex.h barvortex_impl.h)
add_executable(mke_test_baroclin test_baroclin.cpp baroclin.cpp baroclin.h)

target_link_libraries(mke_test_laplace phelm ${LINAL_LIBS})
target_link_libraries(mke_test_chafe phelm ${LINAL_LIBS})

target_link_libraries(mke_test_system_laplace phelm ${LINAL_LIBS})
target_link_libraries(mke_test_jacobian phelm ${LINAL_LIBS})
target_link_libraries(mke_test_barvortex phelm ${LINAL_LIBS})
target_link_libraries(mke_test_baroclin phelm ${LINAL_LIBS})

if (OPENMP_FOUND)
set_target_properties(mke_test_laplace PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(mke_test_chafe PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

set_target_properties(mke_test_system_laplace PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(mke_test_jacobian PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(mke_test_barvortex PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(mke_test_baroclin PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (OPENMP_FOUND)

if (CUDA_FOUND)
	cuda_add_executable(mke_test_laplace_cu ${test_laplace_source})
	target_link_libraries(mke_test_laplace_cu phelm_cu ${LINAL_CUDA_LIBS})
	cuda_add_cublas_to_target(mke_test_laplace_cu)
	set_target_properties(mke_test_laplace_cu PROPERTIES COMPILE_FLAGS "-DGPGPU")

	cuda_add_executable(mke_test_chafe_cu ${test_chafe_source} ${chafe_source})
	target_link_libraries(mke_test_chafe_cu phelm_cu ${LINAL_CUDA_LIBS})
	cuda_add_cublas_to_target(mke_test_chafe_cu)
	set_target_properties(mke_test_chafe_cu PROPERTIES COMPILE_FLAGS "-DGPGPU")

	if (OPENMP_FOUND)
		set_target_properties(mke_test_laplace_cu PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
		set_target_properties(mke_test_chafe_cu PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	endif (OPENMP_FOUND)
endif (CUDA_FOUND)

