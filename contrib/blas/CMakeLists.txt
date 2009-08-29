project(blas)


if (WIN32)

set(LIBRARY_OUTPUT_PATH "../../lib")

add_definitions(-DNO_BLAS_WRAP)
file(GLOB SRC2 BLAS/SRC/*.c)
file(GLOB SRC3 BLAS/WRAP/*.c)


include_directories(.)

#add_library(blas ${SRC2} ${SRC3})

endif(WIN32)
