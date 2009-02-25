project(core CXX C)

include_directories(/usr/include/suitesparse/)
set(LIBRARY_OUTPUT_PATH "../lib")

add_library(mke mke.cpp util.cpp polynom.cpp 
	solver.cpp gmres.c quadrature.c)


