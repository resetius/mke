project(utils CXX C)

set(EXECUTABLE_OUTPUT_PATH "../bin")
add_executable(sphere icosahedron.cpp)
add_executable(rectangle rectangle.cpp)
add_executable(strm2mesh strm2mesh.cpp)

include_directories("../src")

target_link_libraries(strm2mesh mke)

if (OPENMP_FOUND)
set_target_properties(sphere PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(rectangle PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(strm2mesh PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (OPENMP_FOUND)
