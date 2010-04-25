project(utils CXX C)

if (NETCDF_FOUND AND BISON_FOUND)
	include_directories(${NETCDF_INCLUDE_DIR})
	set (BUILD_NETCDF 1)
endif (NETCDF_FOUND AND BISON_FOUND)


set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

include_directories(".")
add_library(mesh_builder mesh_builder.cpp)
add_executable(mke_mesh_sphere icosahedron.cpp)
add_executable(mke_mesh_rectangle rectangle.cpp)
add_executable(mke_topo30_mesh srtm2mesh.cpp)

target_link_libraries(mke_mesh_sphere mesh_builder)
target_link_libraries(mke_mesh_rectangle mesh_builder)

if (BUILD_NETCDF)
	set (netcdf_source netcdf_reader.cpp)
	bison_target(netcdf_parser netcdf_cmd.ypp ${CMAKE_CURRENT_BINARY_DIR}/netcdf_cmd.cpp)
	add_executable(netcdf_reader ${netcdf_source} ${BISON_netcdf_parser_OUTPUTS})
	target_link_libraries(netcdf_reader ${NETCDF_LIBRARIES})
endif (BUILD_NETCDF)

if (OPENMP_FOUND)
set_target_properties(mke_mesh_sphere PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(mke_mesh_rectangle PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(mke_topo30_mesh PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

if (BUILD_NETCDF)
set_target_properties(netcdf_reader PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (BUILD_NETCDF)
endif (OPENMP_FOUND)
