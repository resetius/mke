set(UMFPACK_SOURCE_FOUND 1)

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/UMFPACK")
        set(UMFPACK_SOURCE_FOUND 0)
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/UMFPACK")

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/AMD")
        set(UMFPACK_SOURCE_FOUND 0)
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/AMD")

if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/UFconfig")
        set(UMFPACK_SOURCE_FOUND 0)
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/umfpack/UFconfig")

if (UMFPACK_SOURCE_FOUND)
	message(STATUS "UMFPACK sources found, using it!")
else (UMFPACK_SOURCE_FOUND)
	message(STATUS "UMFPACK sources not found")
endif (UMFPACK_SOURCE_FOUND)

