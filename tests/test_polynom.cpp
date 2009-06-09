#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "polynom.h"

using namespace std;
using namespace MKE;

void check(const string & str1, const string & str2)
{
	if (str1 != str2) {
		fprintf(stderr, "%s != %s\n", str1.c_str(), str2.c_str());
		exit(-1);
	}
}

void test_primitives()
{
	fprintf(stderr, "primitives:\n");
	check(P2X.print(), "+1.0 x^1 ");
	check(P2Y.print(), "+1.0 y^1 ");
	fprintf(stderr, "done\n");
}

void test_primitives_operatins()
{
	fprintf(stderr, "primitives operations:\n");
	Polynom p1 = (P2X - 1) * (2 - 5) - (P2Y - 11) * (-9 - 5);
	Polynom p2 = (P2X - 2) * (2 - 10) - (P2Y - (-11)) * (10 + 5);
	fprintf(stderr, "((x - 1) * (2 - 5):\n");
	check(((P2X - 1) * (2 - 5)).print(), "+3.0 -3.0 x^1 ");
	fprintf(stderr, "(y - 11) * (-9 - 5)):\n");
	check(((P2Y - 11) * (-9 - 5)).print(), "+154.0 -14.0 y^1 ");
	fprintf(stderr, "((x - 1) * (2 - 5) - (y - 11) * (-9 - 5)):\n");
	check(p1.print(), "-151.0 +14.0 y^1 -3.0 x^1 " );
	p1 /= 15.0;
	fprintf(stderr, "((x - 1) * (2 - 5) - (y - 11) * (-9 - 5)) / 15.0:\n");
	check(p1.print(), "-10.1 +0.9 y^1 -0.2 x^1 ");
	fprintf(stderr, "((x - 2) * (2 - 10) - (y + 11) * (10 + 5)):\n");
	check(p2.print(), "-149.0 -15.0 y^1 -8.0 x^1 ");
	fprintf(stderr, "p1 * p2:\n");
	check((p1 * p2).print(), "+1499.9 +11.9 y^1 -14.0 y^2 +110.3 x^1 -4.5 x^1 y^1 +1.6 x^2 ");
	fprintf(stderr, "diff(p1, x):\n");
	check(diff(p1, 0).print(), "-0.2 ");
	fprintf(stderr, "done\n");
}

int main()
{
	test_primitives();
	test_primitives_operatins();
	return 0;
}

