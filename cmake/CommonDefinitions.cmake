if (MSVC)
        add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

add_definitions(${OpenMP_C_FLAGS})
add_definitions(-DSPARSE)

