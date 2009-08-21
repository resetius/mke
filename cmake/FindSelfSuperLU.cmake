set(SUPERLU_SOURCE_FOUND 1)
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/superlu/superlu_timer.c")
	set(SUPERLU_SOURCE_FOUND 0)
	message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/superlu/superlu_timer.c not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/superlu/superlu_timer.c")

if (SUPERLU_SOURCE_FOUND)
	set(SUPERLU superlu blas)
	include_directories(${PHELM_SUBDIR}contrib/superlu)
	add_definitions(-DSUPERLU)
	message(STATUS "SuperLU source found!")
endif (SUPERLU_SOURCE_FOUND)

set(SUPERLU_MT_SOURCE_FOUND 1)
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/superlu_mt/superlu_timer.c")
	set(SUPERLU_MT_SOURCE_FOUND 0)
	message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/superlu_mt/superlu_timer.c not found")
endif (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${PHELM_SUBDIR}/contrib/superlu_mt/superlu_timer.c")

if (SUPERLU_MT_SOURCE_FOUND)
	set(SUPERLU superlu blas) # fixme
	include_directories(${PHELM_SUBDIR}contrib/superlu_mt)
	add_definitions(-DSUPERLU -DSUPERLU_MT)
	message(STATUS "SuperLU_MT source found!")
endif (SUPERLU_MT_SOURCE_FOUND)

if (SUPERLU)
message(STATUS "SuperLU: ${SUPERLU}")
endif (SUPERLU)
