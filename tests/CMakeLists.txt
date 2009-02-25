project(mke CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
include_directories("../src")
link_directories("../lib")

add_executable(test_laplace test_laplace.cpp laplace.cpp)
add_executable(test_slaplace test_slaplace.cpp sphere_laplace.cpp)
add_executable(test_schafe test_schafe.cpp sphere_laplace.cpp)
add_executable(test_polynom test_polynom.cpp)
add_executable(test_integrate test_integrate.cpp)
add_executable(test_gmres test_gmres.cpp)
add_executable(test_quadrature test_quadrature.c)
add_executable(test_funcgen test_funcgen.cpp)

target_link_libraries(test_laplace mke ${UMFPACK})
target_link_libraries(test_slaplace mke ${UMFPACK})
target_link_libraries(test_schafe mke ${UMFPACK})
target_link_libraries(test_polynom mke ${UMFPACK})
target_link_libraries(test_integrate mke ${UMFPACK})
target_link_libraries(test_gmres mke)
target_link_libraries(test_quadrature mke ${MATH})
target_link_libraries(test_funcgen mke ${MATH})

if (DEFINED OMP)
set_target_properties(test_laplace PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_slaplace PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_schafe PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_polynom PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_integrate PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_gmres PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_quadrature PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(test_funcgen PROPERTIES LINK_FLAGS ${OMP})
endif (DEFINED OMP)

