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

double rp(double x, double y)
{
	//return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
	return 2.0 * y * y - 2.0 * y + 2.0 * x * x - 2.0 * x;
}

double ans(double x, double y)
{
	//return sin(M_PI * x) * sin(M_PI * y);// + 1.0;
	return x * (x - 1) * y * (y - 1);// + 1.0;
}

double bnd(double x, double y)
{
	return ans(x, y);
	//return 0.0;
	//return 1.0;
}

double nr2(double * a, double * b, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt(sum);
}

void test_invert(Mesh & mesh)
{
	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();
	int os = (int)mesh.outer.size();

	vec F(sz);
	vec B(os);
	vec Ans(sz);
	vec rans(sz);

	proj(&F[0], mesh, rp);
	proj(&rans[0], mesh, ans);
	proj_bnd(&B[0], mesh, bnd);

	Timer t;
	Laplace l(mesh);
	fprintf(stderr, "l -> %lf\n", t.elapsed());
	l.solve(&Ans[0], &F[0], &B[0]);
	{
		FILE * f = fopen("lu_1_real.txt", "w");
		print_function(f, &rans[0], mesh);
		fclose(f);
		f = fopen("lu_1_calc.txt", "w");
		print_function(f, &Ans[0], mesh);
		fclose(f);
	}

	fprintf(stdout, "L1: invert  err=%.2le\n", 
		dist(&Ans[0], &rans[0], mesh));
}

void test_laplace(Mesh & mesh)
{
	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();
	int os = (int)mesh.outer.size();

	vec U(sz);
	vec LU(sz);
	vec LU1(sz);
	vec B1(os);
	vec B2(os);
	vec P(rs);
	vec P1(rs);

	fprintf(stderr, "prepare ... \n");

	proj(&U[0], mesh, ans);
	proj(&LU[0], mesh, rp);
	proj_bnd(&B1[0], mesh, rp);
	proj_bnd(&B2[0], mesh, ans);

	fprintf(stderr, "start ... \n");

	Laplace l(mesh);
	l.calc1(&LU1[0], &U[0], &B1[0]);
	//l.calc1(&LU1[0], &U[0], &B2[0]);

	fprintf(stdout, "L2: laplace err=%.2le\n", dist(&LU[0], &LU1[0], mesh));
	{
		FILE * f = fopen("lu_real.txt", "w");
		print_inner_function (f, &LU[0], mesh);
		fclose(f);
		f = fopen("lu_calc.txt", "w");
		print_inner_function (f, &LU1[0], mesh);
		fclose(f);

		f = fopen("lu_diff.txt", "w");
		vec tmp(LU.size());
		vec_diff(&tmp[0], &LU[0], &LU1[0], (int)LU.size());
		print_inner_function (f, &tmp[0], mesh);
		fclose(f);
	}

	l.solve(&LU[0], &LU1[0], &B2[0]);
	
//	vector_print(&U[0], U.size());
//	vector_print(&LU[0], LU.size());
	
	fprintf(stdout, "L3: laplace err=%.2le\n", dist(&U[0], &LU[0], mesh));
}

//#include <omp.h>

int main(int argc, char *argv[])
{
	Mesh mesh;

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

	//int nprocs = 1;
	//omp_set_num_threads(nprocs);

	mesh.info();
	//test_invert(mesh);
	//getchar();

	phelm_init();
	test_laplace(mesh);
	phelm_shutdown();
	return 0;
}

