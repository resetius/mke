enable_testing()

IF(BUILD_TESTING)
	SET(BUILDNAME "${BUILDNAME}" CACHE STRING "MKE")
	MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)
