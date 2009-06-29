/* -*- charset: utf-8 -*- */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vector>

#include "phelm.h"
#include "util.h"
#include "laplace.h"

using namespace std;
using namespace phelm;

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [mesh.txt|-]\n", name);
	exit(1);
}

double ans(double x, double y, double t)
{
	//return exp(t) * sin(y) * sin(2.0 * x);
	return x*sin(y+t)*ipow(cos(x),4);
}

double bnd(double x, double y, double t)
{
	return ans(x, y, t);
}

double nr2(double * a, double * b, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt(sum);
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

int main(int argc, char *argv[])
{
	Mesh mesh;
	int i, steps = 10000;
	double tau   = 0.01;
	double mu    = 1.0;
	double sigma = +70;

	if (argc > 1) {
		FILE * f = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "rb");
		if (!f) {
			usage(argv[0]);
		}
		if (!mesh.load(f)) {
			usage(argv[0]);
		}
		fclose(f);
	} else {
		usage(argv[0]);
	}

	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();

	vector < double > U(sz);
	vector < double > B(mesh.outer.size());
	vector < double > Ans(sz);
	vector < double > P(rs);

	proj(&U[0], mesh, ans, 0.0);

//	print_function(stdout, &U[0], mesh, x, y, z);
//	fflush(stdout);

	SphereChafe schafe(mesh, tau, sigma, mu);

	for (i = 0; i < steps; ++i) {
		proj_bnd(&B[0], mesh, bnd, tau * (i + 1));
		schafe.solve(&U[0], &U[0], &B[0], tau * (i));

		// check
		{
			proj(&Ans[0], mesh, ans, tau * (i + 1));
			fprintf(stderr, "time %lf/ norm %le\n", tau * (i + 1), 
				dist(&U[0], &Ans[0], mesh, sphere_scalar_cb, (void*)0));
//			vector_print(&U[0], U.size());
//			vector_print(&Ans[0], U.size());

			//u2p(&P[0], &U[0], mesh);
			//vector_print(&P[0], P.size());
			//u2p(&P[0], &Ans[0], mesh);
			//vector_print(&P[0], P.size());
		}

		fprintf(stderr, "OUTPUT !!! \n");
		print_function(stdout, &U[0], mesh, x, y, z);

//		print_function(stdout, &F[0], mesh, x, y, z);
//		print_function(stdout, &Ans[0], mesh, x, y, z);
//		fflush(stdout);
	}

	return 0;
}

