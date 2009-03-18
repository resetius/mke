#!/bin/sh
rm -f *.bc
for i in `find -name "*.cpp"` ; do echo $i ; llvm-g++ -O3 -I src -emit-llvm $i -DGMRES -DSPARSE -fopenmp -c -o `basename $i .cpp`.bc ; done
for i in `find -name "*.c"` ; do echo $i ; llvm-gcc -O3 -I src -emit-llvm $i -DGMRES -DSPARSE -fopenmp -c -o `basename $i .c`.bc ; done
rm -f *.s 
for i in `find -name "*.bc"` ; do echo $i ; opt -O3 -f $i -o `basename $i .bc`-opt.bc ; done
for i in `find -name "*-opt.bc"` ; do echo $i ; llc $i -o `basename $i -opt.bc`.s ; done

g++ -fopenmp rectangle.s -o rectangle
g++ -fopenmp icosahedron.s -o icosahedron
g++ -fopenmp mke.s util.s polynom.s solver.s gmres.s quadrature.s test_laplace.s laplace.s -o test_laplace
llvm-link mke.bc util.bc polynom.bc solver.bc gmres.bc quadrature.bc test_laplace.bc laplace.bc -f -o test_laplace.bc

