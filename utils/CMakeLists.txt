project(utils CXX C)

if (NETCDF_FOUND AND BISON_FOUND)
	include_directories(${NETCDF_INCLUDE_DIR})
	set (BUILD_NETCDF 1)
endif (NETCDF_FOUND AND BISON_FOUND)


set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

include_directories(".")
include_directories(${LINAL_INCLUDE})
include_directories("../src")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(mesh_builder mesh_builder.cpp)
add_executable(mke_mesh_sphere icosahedron.cpp)
add_executable(mke_mesh_rectangle rectangle.cpp)
add_executable(mke_topo30_mesh srtm2mesh.cpp)
add_executable(mke_dump_operator_matrix dump_operator_matrix.cpp)

target_compile_features(mesh_builder PUBLIC cxx_range_for)

target_link_libraries(mke_mesh_sphere mesh_builder)
target_link_libraries(mke_mesh_rectangle mesh_builder)
target_link_libraries(mke_dump_operator_matrix phelm ${LINAL_LIBS})

if (BUILD_NETCDF)
	set (netcdf_source netcdf_reader.cpp)
	bison_target(netcdf_parser netcdf_cmd.ypp ${CMAKE_CURRENT_BINARY_DIR}/netcdf_cmd.cpp)
	add_executable(netcdf_reader ${BISON_netcdf_parser_OUTPUTS} ${netcdf_source})
	target_link_libraries(netcdf_reader ${NETCDF_LIBRARIES})
endif (BUILD_NETCDF)
