/* -*- charset: utf-8 -*- */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vector>

#include "util.h"
#include "schafe.h"
#include "norm.h"

using namespace std;
using namespace phelm;

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [-f|--file mesh.txt|-] [-d|--double] [-t|--threads number]\n", name);
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

template < typename T >
void calc_schafe(Mesh & mesh)
{
	int i, steps = 100;
	double tau   = 0.01;
	double mu    = 1.0;
	double sigma = +70;


	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();

	ArrayHost < T > U(sz);
	ArrayHost < T > B(mesh.outer.size());
	ArrayHost < T > Ans(sz);
	ArrayHost < T > P(rs);

	ArrayDevice < T > cU(sz);
	ArrayDevice < T > cB(mesh.outer.size());
	ArrayDevice < T > cAns(sz);
	ArrayDevice < T > cP(rs);

	SphereNorm < T > nr(mesh);

	proj(&U[0], mesh, ans, 0.0);
	vec_copy_from_host(&cU[0],    &U[0], sz);

//	print_function(stdout, &U[0], mesh, x, y, z);
//	fflush(stdout);

	SphereChafe < T > schafe(mesh, tau, sigma, mu);

	for (i = 0; i < steps; ++i) {
		Timer t;

		proj_bnd(&B[0], mesh, bnd, tau * (i + 1));
		vec_copy_from_host(&cB[0], &B[0], (int)B.size());
		schafe.solve(&cU[0], &cU[0], &cB[0], tau * (i));

		// check
		{
			proj(&Ans[0], mesh, ans, tau * (i + 1));
			vec_copy_from_host(&cAns[0], &Ans[0], (int)Ans.size());
			fprintf(stderr, "time %lf/ norm %le: %lf\n", tau * (i + 1), 
				(double)nr.dist(&cU[0], &cAns[0]), t.elapsed());
//			vector_print(&U[0], U.size());
//			vector_print(&Ans[0], U.size());

			//u2p(&P[0], &U[0], mesh);
			//vector_print(&P[0], P.size());
			//u2p(&P[0], &Ans[0], mesh);
			//vector_print(&P[0], P.size());
		}

//		fprintf(stderr, "OUTPUT !!! \n");
//		vec_copy_from_device(&cU[0], &U[0], (int)U.size());
//		print_function(stdout, &U[0], mesh, x, y, z);

//		print_function(stdout, &F[0], mesh, x, y, z);
//		print_function(stdout, &Ans[0], mesh, x, y, z);
//		fflush(stdout);
	}

}

int main(int argc, char *argv[])
{
	Mesh mesh;
	bool use_double = false;

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
		} else if (!strcmp(argv[i], "--threads") || !strcmp(argv[i], "-t")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			int threads = atoi(argv[i + 1]);
			set_num_threads(threads);
		} else if (!strcmp(argv[i], "--double") || !strcmp(argv[i], "-d")) {
			use_double = true;
		}
	}

	if (mesh.ps.empty()) {
		usage(argv[0]);
	}

	mesh.info();

	phelm_init();

	Timer t;
	try {
		if (check_device_supports_double() && use_double)
		{
			fprintf(stderr, "using double\n");
			calc_schafe < float > (mesh);
		}
		else
		{
			fprintf(stderr, "using float\n");
			calc_schafe < float > (mesh);
		}
	} catch (const std::exception & e) {
		fprintf(stderr, "exception: %s\n", e.what());
	}
	fprintf(stderr, "elapsed: %lf\n", t.elapsed());

	phelm_shutdown();

	return 0;
}

