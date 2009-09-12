/* -*- charset: utf-8 -*- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "phelm.h"
#include "util.h"

using namespace std;
using namespace phelm;

static double ans(double x, double y, double t)
{
	//return sin(x + t) * cos(y + t);
	return x * (double)rand() / (double)RAND_MAX +
		y * (double)rand() / (double)RAND_MAX;
}

template < typename T >
void init_func(Mesh & mesh, vector < double > & F, double t, T f)
{
	F.resize(mesh.ps.size());
	for (size_t i = 0; i < mesh.ps.size(); ++i)
	{
		F[i] = f(mesh.ps[i].x(), mesh.ps[i].y(), t);
	}
}

static double x(double u, double v)
{
	return cos(u) * cos(v);
}

static double y(double u, double v)
{
	return cos(u) * sin(v);
}

static double z(double u, double v)
{
	return sin(u);
}

static void usage(const char * name)
{
	fprintf(stderr, "usage: %s [mesh.txt|-]\n", name);
	exit(1);
}

extern "C" int test_funcgen(int argc, char * argv[])
{
	Mesh mesh;
	vector < double > F;

	if (argc > 1) {
		FILE * f = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "rb");
		if (!f) {
			usage(argv[0]);
		}
		mesh.load(f);
		fclose(f);
	} else {
		usage(argv[0]);
	}

	for (int i = 1; i < 10000; ++i) {
		init_func(mesh, F, (double)i * 0.01, ans);
		print_function(stdout, &F[0], mesh, x, y, z);
	}
	return 0;
}
