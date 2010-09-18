if (NOT _FINDMATH_)

if (MSVC)
	set(MATH "")
else (MSVC)
	set(MATH "m")
endif (MSVC)

message(STATUS "Math library: ${MATH}")

set (_FINDMATH_ TRUE)
endif ()
