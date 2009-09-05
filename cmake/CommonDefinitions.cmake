if (MSVC)
        add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

#set(OpenMP_C_FLAGS "${OpenMP_C_FLAGS} -m32")
add_definitions(${OpenMP_C_FLAGS})
add_definitions(-DSPARSE)

