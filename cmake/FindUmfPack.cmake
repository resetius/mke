if (UMFPACK_SOURCE_FOUND)
	include_directories(${PHELM_SUBDIR}contrib/umfpack/UMFPACK/Include)
	include_directories(${PHELM_SUBDIR}contrib/umfpack/AMD/Include)
	include_directories(${PHELM_SUBDIR}contrib/umfpack/UFconfig)
	set(UMFPACK umfpack)
else (UMFPACK_SOURCE_FOUND)
	# UMFPACK_LIBRARY_PATH, system libs
	set(UMFPACK_LIBRARY_PATH /usr/lib;/usr/lib64)
	if (ICC_FOUND)
		include_directories(/usr/include/suitesparse/)
		set(UMFPACK umfpack amd blas ${FLIB})
		find_library(UMFPACK_LIB umfpack HINTS UMFPACK_LIBRARY_PATH)
		find_library(AMD HINTS amd UMFPACK_LIBRARY_PATH)
	else (ICC_FOUND)
		if (NOT WIN32)
			include_directories(/usr/include/suitesparse/)
			set(UMFPACK umfpack amd)
			find_library(UMFPACK_LIB umfpack HINTS UMFPACK_LIBRARY_PATH)
			find_library(AMD amd HINTS UMFPACK_LIBRARY_PATH)
		endif (NOT WIN32)
	endif (ICC_FOUND)
	if (UMFPACK_LIB AND AMD)
		message(STATUS "UMFPACK:${UMFPACK_LIB},${AMD}")
	else (UMFPACK_LIB AND AMD)
		set(UMFPACK "")
	endif (UMFPACK_LIB AND AMD)
endif (UMFPACK_SOURCE_FOUND)

if (UMFPACK)
message(STATUS "UMFPACK library: ${UMFPACK}")
else (UMFPACK)
add_definitions(-DGMRES)
endif (UMFPACK)

