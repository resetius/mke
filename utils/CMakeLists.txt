project(utils CXX C)

if (NETCDF_FOUND AND BISON_FOUND)
	include_directories(${NETCDF_INCLUDE_DIR})
	set (BUILD_NETCDF 1)
endif (NETCDF_FOUND AND BISON_FOUND)

set(EXECUTABLE_OUTPUT_PATH "../bin")
set(LIBRARY_OUTPUT_PATH "../lib")

add_library(mesh_builder mesh_builder.cpp)
add_executable(sphere icosahedron.cpp)
add_executable(rectangle rectangle.cpp)
add_executable(srtm2mesh srtm2mesh.cpp)
add_executable(analizer analizer.cpp)

target_link_libraries(sphere mesh_builder)
target_link_libraries(rectangle mesh_builder)

if (BUILD_NETCDF)
	set (netcdf_source netcdf_reader.cpp)
	bison_target(netcdf_parser netcdf_cmd.ypp ${CMAKE_CURRENT_BINARY_DIR}/netcdf_cmd.cpp)
	add_executable(netcdf_reader ${netcdf_source} ${BISON_netcdf_parser_OUTPUTS})
	target_link_libraries(netcdf_reader netcdf_c++ netcdf)
endif (BUILD_NETCDF)

if (OPENMP_FOUND)
set_target_properties(sphere PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(rectangle PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(srtm2mesh PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(analizer PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})

if (BUILD_NETCDF)
set_target_properties(netcdf_reader PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif (BUILD_NETCDF)
endif (OPENMP_FOUND)
