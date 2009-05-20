set (ICC_PATTERN ".*icpc.*")

if (NOT WIN32)
	if (NOT CYGWIN)
		if (CMAKE_CXX_COMPILER MATCHES ${ICC_PATTERN})
			set(ICC_LIB_PATH "")
			set(ICC_COMPILE "1")
			set(ICC_PATH ${CMAKE_CXX_COMPILER})
			string(REPLACE "/bin/"  "/lib/" ICC_LIB_PATH ${ICC_PATH})
			string(REPLACE "/icpc/" "/"     ICC_LIB_PATH ${ICC_LIB_PATH})

			link_directories(${ICC_LIB_PATH})
			set(FLIB "ifcoremt")
		else (CMAKE_CXX_COMPILER MATCHES ${ICC_PATTERN})
			set(FLIB "gfortran")
		endif (CMAKE_CXX_COMPILER MATCHES ${ICC_PATTERN})
	endif (NOT CYGWIN)
endif (NOT WIN32)

message(STATUS "Fortran library: ${FLIB}")

