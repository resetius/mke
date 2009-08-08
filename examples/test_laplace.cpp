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

template < typename T >
void test_invert(Mesh & mesh)
{
	typedef Array < T, Allocator < T > > vector;

	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();
	int os = (int)mesh.outer.size();

	std::vector < T > F(sz);
	std::vector < T > B(os);
	std::vector < T > Ans(sz);
	std::vector < T > rans(sz);

	vector cF(sz);
	vector cB(os);
	vector cAns(sz);
	vector crans(sz);

	proj(&F[0], mesh, rp);
	proj(&rans[0], mesh, ans);
	proj_bnd(&B[0], mesh, bnd);

	vec_copy_from_host(&cF[0], &F[0], sz);
	vec_copy_from_host(&crans[0], &rans[0], sz);
	vec_copy_from_host(&cB[0], &B[0], os);

	Timer t;
	Laplace < T > l(mesh);
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

template < typename T > 
void test_laplace(Mesh & mesh)
{
	typedef Array < T, Allocator < T > > vector;
	int sz = (int)mesh.ps.size();
	int rs = (int)mesh.inner.size();
	int os = (int)mesh.outer.size();
	FlatNorm < T > nr(mesh);

	std::vector < T > U(sz);
	std::vector < T > LU(sz);
	std::vector < T > LU1(sz);
	std::vector < T > B1(os);
	std::vector < T > B2(os);

	vector cU(sz);
	vector cLU(sz);
	vector cLU1(sz);
	vector cB1(os);
	vector cB2(os);

	fprintf(stderr, "prepare ... \n");

	proj(&U[0], mesh, ans);
	proj(&LU[0], mesh, rp);
	proj_bnd(&B1[0], mesh, rp);
	proj_bnd(&B2[0], mesh, ans);

	vec_copy_from_host(&cU[0], &U[0], sz);
	vec_copy_from_host(&cLU[0], &LU[0], sz);
	vec_copy_from_host(&cB1[0], &B1[0], os);
	vec_copy_from_host(&cB2[0], &B2[0], os);

	fprintf(stderr, "start ... \n");

	Laplace < T > l(mesh);
	l.calc1(&cLU1[0], &cU[0], &cB1[0]);

	vec_copy_from_device(&LU1[0], &cLU1[0], sz);

	fprintf(stdout, "L2: laplace err=%.2le\n", 
		nr.dist(&LU[0], &LU1[0]));
	{
		FILE * f = fopen("lu_real.txt", "w");
		print_inner_function (f, &LU[0], mesh);
		fclose(f);
		f = fopen("lu_calc.txt", "w");
		print_inner_function (f, &LU1[0], mesh);
		fclose(f);

		f = fopen("lu_diff.txt", "w");
		vector tmp1(sz);
		std::vector < T > tmp2(sz);

		vec_diff(&tmp1[0], &cLU[0], &cLU1[0], sz);
		vec_copy_from_device(&tmp2[0], &tmp1[0], sz);

		print_inner_function (f, &tmp2[0], mesh);
		fclose(f);
	}

	l.solve(&cLU[0], &cLU1[0], &cB2[0]);
	vec_copy_from_device(&LU[0], &cLU[0], sz);
	
//	vector_print(&U[0], U.size());
//	vector_print(&LU[0], LU.size());
	
	fprintf(stdout, "L3: laplace err=%.2le\n", nr.dist(&U[0], &LU[0]));
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


	phelm_init();

	//test_invert < float > (mesh);
	//getchar();

	test_laplace < float > (mesh);
	//test_laplace < double > (mesh);
	
	phelm_shutdown();
	return 0;
}
