set(UMFPACK_SOURCE_FOUND 1)
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/UMFPACK")
        set(UMFPACK_SOURCE_FOUND 0)
	message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/UMFPACK not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/UMFPACK")

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/AMD")
        set(UMFPACK_SOURCE_FOUND 0)
	message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/AMD not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/AMD")

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/UFconfig")
        set(UMFPACK_SOURCE_FOUND 0)
	message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/UFConfig not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/umfpack/UFconfig")

if (UMFPACK_SOURCE_FOUND)
	message(STATUS "UMFPACK sources found, using it!")
else (UMFPACK_SOURCE_FOUND)
	message(STATUS "UMFPACK sources not found")
endif (UMFPACK_SOURCE_FOUND)

