project(phelm CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
include_directories(../src ../auxiliary)
link_directories(../lib)

create_test_sourcelist(small_tests test_small.c test_polynom.cpp test_integrate.cpp test_quadrature.c test_funcgen.cpp test_diff.cpp)

add_executable(test_small ${small_tests})
add_executable(test_solver test_solver.cpp)

target_link_libraries(test_small phelm phelm_aux ${MATH} ${UMFPACK} ${SUPERLU})

target_link_libraries(test_solver phelm phelm_aux ${UMFPACK} ${SUPERLU})

if (OPENMP_FOUND)
set_target_properties(test_small PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

set_target_properties(test_solver PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (OPENMP_FOUND)

if (WIN32)
set (LUA ${EXECUTABLE_OUTPUT_PATH}/lua)
else (WIN32)
set (LUA lua)
endif (WIN32)

add_test(test_polynom ${EXECUTABLE_OUTPUT_PATH}/test_small test_polynom)
add_test(test_quadrature ${EXECUTABLE_OUTPUT_PATH}/test_small test_quadrature)
add_test(test_integrate ${EXECUTABLE_OUTPUT_PATH}/test_small test_integrate)
add_test(test_laplace ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH})
add_test(test_system_laplace ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_system_laplace.lua ${EXECUTABLE_OUTPUT_PATH})
add_test(test_slaplace_full ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} full)
add_test(test_slaplace_half ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} half)
add_test(test_slaplace_test ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} test)
add_test(test_slaplace_test2 ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} test2)
add_test(test_barvortex_full ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_barvortex.lua ${EXECUTABLE_OUTPUT_PATH} full)
add_test(test_barvortex_half ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_barvortex.lua ${EXECUTABLE_OUTPUT_PATH} half)
add_test(test_jacobian_full ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} full)
add_test(test_jacobian_half ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} half)
add_test(test_jacobian_test ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} test)
add_test(test_jacobian_test2 ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} test2)

if (CUDA_FOUND)
	cuda_add_executable(test_gpu_reduction test_gpu_reduction.cpp)
	target_link_libraries(test_gpu_reduction phelm_cu)
	cuda_add_cublas_to_target(test_gpu_reduction)
	set_target_properties(test_gpu_reduction PROPERTIES COMPILE_FLAGS "-DGPGPU")

	if (OPENMP_FOUND)
		set_target_properties(test_gpu_reduction PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	endif (OPENMP_FOUND)

	add_test(test_gpu_reduction ${EXECUTABLE_OUTPUT_PATH}/test_gpu_reduction)
	add_test(test_laplace_cu ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH} test_laplace_cu)
endif (CUDA_FOUND)

