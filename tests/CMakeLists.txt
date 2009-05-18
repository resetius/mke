project(mke CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
include_directories(../src ../auxiliary)
link_directories(../lib)

add_executable(test_polynom test_polynom.cpp)
add_executable(test_integrate test_integrate.cpp)
add_executable(test_gmres test_gmres.cpp)
add_executable(test_quadrature test_quadrature.c)
add_executable(test_funcgen test_funcgen.cpp)

target_link_libraries(test_polynom mke mke_aux ${UMFPACK})
target_link_libraries(test_integrate mke mke_aux ${UMFPACK})
target_link_libraries(test_gmres mke mke_aux)
target_link_libraries(test_quadrature mke mke_aux ${MATH})
target_link_libraries(test_funcgen mke mke_aux ${MATH} ${UMFPACK})

if (OPENMP_FOUND)
set_target_properties(test_polynom PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_integrate PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_gmres PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_quadrature PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(test_funcgen PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (OPENMP_FOUND)

if (WIN32)
set (LUA ${EXECUTABLE_OUTPUT_PATH}/lua)
else (WIN32)
set (LUA lua)
endif (WIN32)

add_test(test_polynom ${EXECUTABLE_OUTPUT_PATH}/test_polynom)
add_test(test_quadrature ${EXECUTABLE_OUTPUT_PATH}/test_quadrature)
add_test(test_integrate ${EXECUTABLE_OUTPUT_PATH}/test_integrate)
add_test(test_laplace ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH})
add_test(test_slaplace_full ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH} full)
add_test(test_slaplace_half ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH} half)
add_test(test_slaplace_test ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH} test)
add_test(test_slaplace_test2 ${LUA} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/test_laplace.lua ${EXECUTABLE_OUTPUT_PATH} test2)


