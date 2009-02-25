project(utils CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
add_executable(sphere icosahedron.cpp)
add_executable(rectangle rectangle.cpp)


if (DEFINED OMP)
set_target_properties(rectangle PROPERTIES LINK_FLAGS ${OMP})
set_target_properties(sphere PROPERTIES LINK_FLAGS ${OMP})
endif (DEFINED OMP)


