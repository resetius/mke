project(utils CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
add_executable(sphere icosahedron.cpp)
add_executable(rectangle rectangle.cpp)

if (DEFINED OMP_LFLAGS)
set_target_properties(sphere PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(rectangle PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (DEFINED OMP_LFLAGS)

