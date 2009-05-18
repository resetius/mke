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

void check(double a, double b)
{
	if (fabs(a - b) > 1e-14) {
		fprintf(stderr, "%.16lf != %.16lf\n", a, b);
		exit(-1);
	}
}

int main()
{
	double i;

	i = gauss_kronrod15(0, 1, f1, 0);
	fprintf(stderr, "int 0 1 x * x = %.16le\n", i);
	check(i, 0.3333333333333333333);


	i = gauss_kronrod15(-1, 1, f2, 0);
	fprintf(stderr, "int -1 1 x / cos(x) = %.16le\n", i);
	check(i, 0);
	return 0;
}

