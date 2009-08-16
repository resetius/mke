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
	fprintf(stderr, "usage: %s [-f|--file mesh.txt|-] [-d|--double] [-t|--threads number]\n", name);
	exit(1);
}

double rp(double x, double y)
{
	//return -6.0 * sin(y) * sin(2.0 * x);
	
	double sx = sin(x);
	double cx = cos(x);
	double sy = sin(y);
	double cy = cos(y);

	double ans = 0.0;

	ans += (sx * sx * (cx * cy - 1) * cy);
	ans += - (cx * cy - 1) * cy;
	ans += -2 * sx * (cx * cy - sx) * (-sx * cy - cx);
	ans += -2 * (cx * cy - sx) * cy;
	if (fabs(ans) > 1e-15) {
		ans /= cx;
	}

	ans += sx * sx * cy * cy;
	ans += -cx * (cx * cy - 1) * cy;
	ans += sy * sy;
	ans += 2 * ipow(-sx * cy - cx, 2);
	ans += 2 * (cx * cy - sx) * (-cx * cy + sx);
	ans += 2 * sy * sy;
	ans += -8 * (sx - 3.0) * sx + 4.0 * cx * cx;
	return ans;
}

double ans(double x, double y)
{
	//return sin(y) * sin(2.0 * x);

	double sx = sin(x);
	double cx = cos(x);
	//double sy = sin(y);
	double cy = cos(y);

	return 0.5 * ipow(cx * cy - 1.0, 2) +
		ipow(cx * cy - sx, 2) +
		2.0 * ipow(sx - 3, 2);
}

double bnd(double x, double y)
{
	return ans(x, y);
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
void test_invert(const Mesh & mesh)
{
	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();
	int os = (int)mesh.outer.size();

	ArrayHost < T > F(sz);
	ArrayHost < T > B(std::max(os, 1));
	ArrayHost < T > Ans(sz);
	ArrayHost < T > rans(sz);

	ArrayDevice < T > cF(sz);
	ArrayDevice < T > cB(std::max(os, 1));
	ArrayDevice < T > cAns(sz);
	ArrayDevice < T > crans(sz);

	SphereNorm < T > nr(mesh);

	proj(&F[0], mesh, rp);
	proj(&rans[0], mesh, ans);
	proj_bnd(&B[0], mesh, bnd);

	vec_copy_from_host(&cF[0],    &F[0], sz);
	vec_copy_from_host(&cB[0],    &B[0], (int)B.size());
	vec_copy_from_host(&crans[0], &rans[0], sz);

	Timer t;
	SphereLaplace < T > l(mesh);
	fprintf(stderr, "l->%lf\n", t.elapsed());
	l.solve(&cAns[0], &cF[0], &cB[0]);

	fprintf(stdout, "L1: invert  err=%.2le\n", (double)nr.dist(&cAns[0], &crans[0]));

	{
		FILE * f = fopen("invert_answer.txt", "wb");
		vec_copy_from_device(&Ans[0], &cAns[0], sz);
		print_function(f, &Ans[0], mesh, x, y, z);
		fclose(f);

		fprintf(stderr, "invert answer saved to 'invert_answer.txt'\n");

		f = fopen("answer.txt", "wb");
		print_function(f, &rans[0], mesh, x, y, z);
		fclose(f);

		f = fopen("rp.txt", "wb");
		print_function(f, &F[0], mesh, x, y, z);
		fclose(f);
	}
}

template < typename T >
void test_laplace(const Mesh & mesh)
{
	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();
	int os = (int)mesh.outer.size();

	ArrayHost < T > U(sz);
	ArrayHost < T > LU(sz);
	ArrayHost < T > LU1(sz);
	ArrayHost < T > B1(std::max(os, 1));
	ArrayHost < T > B2(std::max(os, 1));

	ArrayDevice < T > cU(sz);
	ArrayDevice < T > cLU(sz);
	ArrayDevice < T > cLU1(sz);
	ArrayDevice < T > cB1(std::max(os, 1));
	ArrayDevice < T > cB2(std::max(os, 1));

	SphereNorm < T > nr(mesh);

	proj(&U[0], mesh, ans);
	proj(&LU[0], mesh, rp);
	proj_bnd(&B1[0], mesh, rp);
	proj_bnd(&B2[0], mesh, ans);

	vec_copy_from_host(&cU[0], &U[0], sz);
	vec_copy_from_host(&cLU[0], &LU[0], sz);
	vec_copy_from_host(&cB1[0], &B1[0], (int)B1.size());
	vec_copy_from_host(&cB2[0], &B2[0], (int)B2.size());

	SphereLaplace < T > l(mesh);
	l.calc1(&cLU1[0], &cU[0], &cB1[0]);

	fprintf(stdout, "L2: laplace err=%.2le\n", (double)nr.dist(&cLU[0], &cLU1[0]));
	{
		FILE * f = fopen("slu_real.txt", "w");
		print_function(f, &LU[0], mesh, x, y, z);
		fclose(f);
		f = fopen("slu_calc.txt", "w");
		vec_copy_from_device(&LU1[0], &LU1[0], sz);
		print_function(f, &LU1[0], mesh, x, y, z);
		fclose(f);
	}

	l.solve(&cLU[0], &cLU1[0], &cB2[0]);

	fprintf(stdout, "L3: laplace err=%.2le\n", (double)nr.dist(&cU[0], &cLU[0]));
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


	mesh.info();
	Timer t;
	try {
		if (check_device_supports_double() && use_double)
		{
			fprintf(stderr, "using double\n");
			test_invert < double > (mesh);
			test_laplace < double > (mesh);
		}
		else
		{
			fprintf(stderr, "using float\n");
			test_invert < float > (mesh);
			test_laplace < float > (mesh);
		}
	} catch (const std::exception & e) {
		fprintf(stderr, "exception: %s\n", e.what());
	}

	fprintf(stderr, "elapsed: %lf\n", t.elapsed());

	phelm_shutdown();

	return 0;
}

