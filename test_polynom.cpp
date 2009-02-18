#include <stdio.h>
#include <stdlib.h>

#include "polynom.h"

using namespace std;

void test_primitives()
{
	fprintf(stderr, "primitives:\n");
	P2X.print();
	P2Y.print();
	fprintf(stderr, "done\n");
}

void test_primitives_operatins()
{
	fprintf(stderr, "primitives operations:\n");
	Polynom p1 = (P2X - 1) * (2 - 5) - (P2Y - 11) * (-9 - 5);
	Polynom p2 = (P2X - 2) * (2 - 10) - (P2Y - (-11)) * (10 + 5);
	p1 /= 15.0;
	fprintf(stderr, "((x - 1) * (2 - 5) - (y - 11) * (-9 - 5)) / 15.0:\n");
	p1.print();
	fprintf(stderr, "((x - 2) * (2 - 10) - (y + 11) * (10 + 5)):\n");
	p2.print();
	fprintf(stderr, "p1 * p2:\n");
	(p1 * p2).print();
	fprintf(stderr, "done\n");
}

int main()
{
	test_primitives();
	test_primitives_operatins();
	return 0;
}

