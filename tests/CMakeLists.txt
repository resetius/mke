project(phelm CXX C)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
include_directories(../src)
include_directories(../utils)
include_directories(${LINAL_INCLUDE})

create_test_sourcelist(small_tests 
	test_small.c test_polynom.cpp 
	test_integrate.cpp test_quadrature.c 
	test_funcgen.cpp test_diff.cpp 
	test_triag_rotate.cpp
	test_func_op.cpp
	)

add_executable(mke_test_small ${small_tests})

target_link_libraries(mke_test_small phelm ${MATH} ${LINAL_LIBS} mesh_builder)

# TODO: check lua version (>= 5.3)
find_program(LUA lua)

add_test(mke_func_op ${EXECUTABLE_OUTPUT_PATH}/mke_test_small test_func_op)
add_test(mke_test_polynom ${EXECUTABLE_OUTPUT_PATH}/mke_test_small test_triag_rotate)
add_test(mke_test_polynom ${EXECUTABLE_OUTPUT_PATH}/mke_test_small test_polynom)
add_test(mke_test_quadrature ${EXECUTABLE_OUTPUT_PATH}/mke_test_small test_quadrature)
add_test(mke_test_integrate ${EXECUTABLE_OUTPUT_PATH}/mke_test_small test_integrate)
add_test(mke_test_laplace ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH})
add_test(mke_test_system_laplace ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_system_laplace.lua ${EXECUTABLE_OUTPUT_PATH})
add_test(mke_test_slaplace_full ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} full)
add_test(mke_test_slaplace_half ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} half)
add_test(mke_test_slaplace_test ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} test)
add_test(mke_test_slaplace_test2 ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_slaplace.lua ${EXECUTABLE_OUTPUT_PATH} test2)
add_test(mke_test_barvortex_full ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_barvortex.lua ${EXECUTABLE_OUTPUT_PATH} full)
add_test(mke_test_barvortex_half ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_barvortex.lua ${EXECUTABLE_OUTPUT_PATH} half)
add_test(mke_test_jacobian_full ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} full)
add_test(mke_test_jacobian_half ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} half)
add_test(mke_test_jacobian_test ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} test)
add_test(mke_test_jacobian_test2 ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_jacobian.lua ${EXECUTABLE_OUTPUT_PATH} test2)

if (CUDA_FOUND)
	add_test(mke_test_laplace_cu ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH} mke_test_laplace_cu)
endif (CUDA_FOUND)

