project(utils CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
add_executable(sphere icosahedron.cpp)
add_executable(rectangle rectangle.cpp)

if (OPENMP_FOUND)
set_target_properties(sphere PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(rectangle PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (OPENMP_FOUND)

