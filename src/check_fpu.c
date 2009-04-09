#include <stdio.h>

int main()
{
	unsigned m;
	asm ("fstcw %0" : : "m" (*&m));
	asm ("fwait");
	printf("%x\n", (int)m);
	m ^= 0x7f;
	printf("%x\n", (int)m);
	asm ("fldcw %0" : : "m" (*&m));
	asm ("fwait");
	m = 0;
	asm ("fstcw %0" : : "m" (*&m));
	asm ("fwait");



	printf("%x\n", (int)m);

	asm ("stmxcsr %0" : : "m" (*&m));
	printf("%x\n", (int)m);
	m = 0;
	asm ("ldmxcsr %0" : : "m" (*&m));

	printf("%lf\n", 0.0 / 0.0);
	printf("%lf\n", 1e100*1e100*100000000000000.0 * 1000000000000000000 * 10000000000000000000.00 * 1000000000000000000.00 * 1e50);
}

