/* -*- charset: utf-8 -*- */
 
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#include "phelm.h"
#include "jacobian.h"
#include "util.h"

using namespace std;
using namespace phelm;

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [-f|--file mesh.txt|-] [--verbose|-v number]\n", name);
	exit(1);
}

double f1 (double x, double y)
{
	return sin (x) * sin (y) * cos(x);
}

double f2 (double x, double y)
{
	return 1000.0 * cos (x) * cos (y) * sin(x);
}

double f3 (double x, double y)
{
	return -(-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) ) / cos (x);
}

double an (double x, double y)
{
	return (-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) ) / cos (x);
	//return (-sin(x)*cos(y)*sin(x)*cos(y)) / cos(x);
	//return (-cos(x)*sin(y)*cos(x)*sin(y))/cos(x);

	//return cos(x) * sin(y) / cos(x);
	//return sin(x) * cos(y) / cos(x);
}

double test_jacobian_f1 (double x, double y)
{
	return sin (x) * sin (y);
}

double test_jacobian_f2 (double x, double y)
{
	return cos (x) * cos (y);
}

double test_jacobian_an (double x, double y)
{
	return (-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) );
}

void test_jacobian (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int rs = (int)m.inner.size();
	int os = (int)m.outer.size();

	Jacobian j (m);
	vector < double > F1 (sz);
	vector < double > F2 (sz);
	vector < double > ans1 (sz);
	vector < double > rans1 (sz);
	vector < double > bnd (os);

	proj (&F1[0], m, test_jacobian_f1);
	proj (&F2[0], m, test_jacobian_f2);
	proj (&rans1[0], m, test_jacobian_an);
	proj_bnd (&bnd[0], m, test_jacobian_an);

	j.calc1 (&ans1[0], &F1[0], &F2[0], &bnd[0]);

	fprintf (stdout, "jacobian err=%.2le\n",
	         dist (&ans1[0], &rans1[0], m, sphere_scalar_cb, (void*)0) );

	//vector < double > p1(m.inner.size());
	//u2p(&p1[0], &rans1[0], m);
	//vector_print(&p1[0], p1.size());
	//u2p(&p1[0], &ans1[0], m);
	//vector_print(&p1[0], p1.size());

	//vector_print(&rans1[0], rans1.size());
	//vector_print(&ans1[0], ans1.size());
}

static double x (double u, double v)
{
	return cos (u) * cos (v);
}

static double y (double u, double v)
{
	return cos (u) * sin (v);
}

static double z (double u, double v)
{
	return sin (u);
}

int main (int argc, char *argv[])
{
	Mesh mesh;
	string task;
	int verbose = 0;

	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
		{
			usage(argv[0]);
		} else if (!strcmp(argv[i], "--file") || !strcmp(argv[i], "-f")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			FILE * f = (strcmp(argv[i + 1], "-") == 0) ? stdin : fopen(argv[i + 1], "rb");

			if (!f) {
				usage(argv[0]);
			}
			if (!mesh.load(f)) {
				usage(argv[0]);
			}

			fclose(f);
		} else if (!strcmp(argv[i], "--verbose") || !strcmp(argv[i], "-v")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}
			verbose = atoi(argv[i + 1]);
		}
	}

	if (mesh.ps.empty()) {
		usage(argv[0]);
	}

#if defined(WIN32)
	set_fpe_except();
#endif

	mesh.info();

	test_jacobian(mesh);
	return 0;
}

