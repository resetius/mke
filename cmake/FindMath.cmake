if (MSVC)
	set(MATH "")
else (MSVC)
	set(MATH "m")
endif (MSVC)

message(STATUS "Math library: ${MATH}")

