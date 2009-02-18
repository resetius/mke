#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"

double f1(double x, void * unused)
{
	return x * x;
}

double f2(double x, void * unused)
{
	return x / cos(x);
}

int main()
{
	double i;

	i = gauss_kronrod15(0, 1, f1, 0);
	printf("int 0 1 x * x = %.16le\n", i);

	i = gauss_kronrod15(-1, 1, f2, 0);
	printf("int -1 1 x / cos(x) = %.16le\n", i);
	return 0;
}

