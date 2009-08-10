/* -*- charset: utf-8 -*- */

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "gmres.h"
#include "util.h"

using namespace phelm;
using namespace std;

bool cmp(double a, double b)
{
	return (fabs(a - b) < 1e-15);
}

bool cmp(float a, float b)
{
	return (fabsf(a - b) < 1e-6f);
}

template < typename T >
bool test_gmres()
{
	int i, j = 0;
	int n  = 300000;
//	int n  = 30;
	int nz = n + n - 1 + n - 1;

	Array < int, Allocator < int > > cAp((n + 1));
	Array < int, Allocator < int > > cAi(nz);
	Array < T, Allocator < T > > cAx(nz);
	Array < T, Allocator < T > > cb(n);
	Array < T, Allocator < T > > cx(n);

	vector < int > Ap((n + 1));
	vector < int > Ai(nz);
	vector < T > Ax(nz);
	vector < T > x(n);
	vector < T > b(n);

	Sparse_t < T > A;

	/**
	 * -4  1  0  ....  0
	 *  1 -4  1  ....  0
	 *  0  1 -4  1 ..  0
	 *  0  0  1 -4 ..  0
	 * ............
	 * ......    1 -4  1
	 * .........   1  -4
	 */

	/* матрицу записываем по строкам! */
	Ap[0] = 0;
	for (i = 0; i < n; ++i) {
		if (i == 0) {
			/* 2 1 0 .... 0 */
			Ap[i + 1] = Ap[i] + 2;

			Ax[j    ] = -4;
			Ax[j + 1] =  1;

			Ai[j    ] = i;
			Ai[j + 1] = i + 1;

			j += 2;
		} else if (i < n - 1) {
			/* 1 2 1 .... 0 */
			Ap[i + 1] = Ap[i] + 3;

			Ax[j    ] =  1;
			Ax[j + 1] = -4;
			Ax[j + 2] =  1;

			Ai[j    ] = i - 1;
			Ai[j + 1] = i;
			Ai[j + 2] = i + 1;

			j += 3;
		} else {
			/* 0 .... 1 2 1 */
			Ap[i + 1] = Ap[i] + 2;

			Ax[j    ] =  1;
			Ax[j + 1] = -4;

			Ai[j    ] = i - 1;
			Ai[j + 1] = i;

			j += 2;
		}

		b[i] = 1;
	}

	vec_copy_from_host(&cAp[0], &Ap[0], (int)Ap.size());
	vec_copy_from_host(&cAi[0], &Ai[0], (int)Ai.size());
	vec_copy_from_host(&cAx[0], &Ax[0], (int)Ax.size());
	vec_copy_from_host(&cb[0], &b[0], (int)b.size());

	A.Ax = &cAx[0];
	A.Ap = &cAp[0];
	A.Ai = &cAi[0];

	for (int k = 0; k < 10; ++k) {
		gmres(&cx[0], &A, &cb[0], 
			MatMultiplier < Sparse_t < T > > :: mult_vector_l, n, 10, 100);

//		vec_copy_from_device(&x[0], &cx[0], n);
//		vec_print(&x[0], n);		
	}
//	vec_copy_from_device(&x[0], &cx[0], n);
//	vec_print(&x[0], n);

	return true;
}

template < typename T >
bool test_matvect()
{
	int i, j = 0;
	int n  = 500000;
//	int n  = 500;
	int nz = n + n - 1 + n - 1;

	Array < int, Allocator < int > > cAp((n + 1));
	Array < int, Allocator < int > > cAi(nz);
	Array < T, Allocator < T > > cAx(nz);
	Array < T, Allocator < T > > cb(n);
	Array < T, Allocator < T > > cx(n);

	vector < int > Ap((n + 1));
	vector < int > Ai(nz);
	vector < T > Ax(nz);
	vector < T > b(n);
	vector < T > x(n);
	Sparse_t < T > A;

	/**
	 * 2 1 0 .... 0
	 * 1 2 1 .... 0
	 * 0 1 2 1 .. 0
	 * 0 0 1 2 .. 0
	 * ............
	 * ...... 1 2 1
	 * .........1 2
	 */

	/* матрицу записываем по строкам! */
	Ap[0] = 0;
	for (i = 0; i < n; ++i) {
		if (i == 0) {
			/* 2 1 0 .... 0 */
			Ap[i + 1] = Ap[i] + 2;

			Ax[j    ] = 2;
			Ax[j + 1] = 1;

			Ai[j    ] = i;
			Ai[j + 1] = i + 1;

			j += 2;
		} else if (i < n - 1) {
			/* 1 2 1 .... 0 */
			Ap[i + 1] = Ap[i] + 3;

			Ax[j    ] = 1;
			Ax[j + 1] = 2;
			Ax[j + 2] = 1;

			Ai[j    ] = i - 1;
			Ai[j + 1] = i;
			Ai[j + 2] = i + 1;

			j += 3;
		} else {
			/* 0 .... 1 2 1 */
			Ap[i + 1] = Ap[i] + 2;

			Ax[j    ] = 2;
			Ax[j + 1] = 1;

			Ai[j    ] = i - 1;
			Ai[j + 1] = i;

			j += 2;
		}

		b[i] = 1;
	}

	vec_copy_from_host(&cAp[0], &Ap[0], (int)Ap.size());
	vec_copy_from_host(&cAi[0], &Ai[0], (int)Ai.size());
	vec_copy_from_host(&cAx[0], &Ax[0], (int)Ax.size());
	vec_copy_from_host(&cb[0], &b[0], (int)b.size());

	A.Ax = &cAx[0];
	A.Ap = &cAp[0];
	A.Ai = &cAi[0];

	for (int k = 0; k < 1000; ++k) {
		sparse_mult_vector_l(&cx[0], &A, &cb[0], n);
//		vec_copy_from_device(&x[0], &cx[0], n);
//		vec_print(&x[0], n);
	}

	vec_copy_from_device(&x[0], &cx[0], n);

	if (!cmp(x[0], (T)3) || !cmp(x[n - 1], (T)3)) {
		fprintf(stderr, "fail on edge: x[0]=%lf != 3, x[n - 1]=%lf != 3\n", 
				(double)x[0], (double)x[n - 1]);
		return false;
	}

	for (int i = 1; i < n - 1; ++i) {
		if (!cmp(x[i], (T)4)) {
			fprintf(stderr, " x[i]=%lf != 4 \n", (double)x[i]);
			return false;
		}
	}
	return true;
}

template < typename T >
bool test_sum()
{
	T s = 0.0;
	int n = 10000;
	Array < T, Allocator < T > > ca(n);
	Array < T, Allocator < T > > cb(n);
	Array < T, Allocator < T > > cr(n);

	Array < T, std::allocator < T > > a(n);
	Array < T, std::allocator < T > > b(n);
	Array < T, std::allocator < T > > r(n);

	for (int i = 0; i < n; ++i) {
		a[i] = (T)1.5;
		b[i] = (T)2.5;
	}

	vec_copy_from_host(&ca[0], &a[0], n);
	vec_copy_from_host(&cb[0], &b[0], n);

	vec_sum2(&cr[0], &ca[0], &cb[0], -1.0, n);
	s = vec_scalar2(&ca[0], &cb[0], n);

	vec_copy_from_device(&r[0], &cr[0], n);
	fprintf(stderr, "%lf\n", (double)s);
	//vec_print(&r[0], n);
	
	for (int i = 0; i < n; ++i) {
		if (!cmp(r[i], (T)-1)) {
			fprintf(stderr, " r[i]=%lf != -1\n", (double)r[i]);
			return false;
		}
	}
	return true;
}

int main(int argc, char * argv[])
{
	bool result = true;
	try {
		phelm_init();

		int has_double = check_device_supports_double();
		fprintf(stderr, "has double: %d\n", has_double);

		Timer t;
		if (has_double) {
			fprintf(stderr, "testing double:\n");

			t.restart(); result &= test_sum < float > ();
			fprintf(stderr, "test_sum < float > (): %lf, %d\n", t.elapsed(), (int)result);
			t.restart(); result &= test_gmres < float > ();
			fprintf(stderr, "test_gmres < float > (): %lf, %d\n", t.elapsed(), (int)result);
			t.restart(); result &= test_matvect < float > ();
			fprintf(stderr, "test_matvect < float > (): %lf, %d\n", t.elapsed(), (int)result);

			t.restart(); result &= test_sum < double > ();
			fprintf(stderr, "test_sum < double > (): %lf, %d\n", t.elapsed(), (int)result);
			t.restart(); result &= test_gmres < double > ();
			fprintf(stderr, "test_gmres < double > (): %lf, %d\n", t.elapsed(), (int)result);
			t.restart(); result &= test_matvect < double > ();
			fprintf(stderr, "test_matvect < double > (): %lf, %d\n", t.elapsed(), (int)result);

		} else {
			t.restart(); result &= test_sum < float > ();
			fprintf(stderr, "test_sum < float > (): %lf, %d\n", t.elapsed(), (int)result);
			t.restart(); result &= test_gmres < float > ();
			fprintf(stderr, "test_gmres < float > (): %lf, %d\n", t.elapsed(), (int)result);
			t.restart(); result &= test_matvect < float > ();
			fprintf(stderr, "test_matvect < float > (): %lf, %d\n", t.elapsed(), (int)result);
		}
		fprintf(stderr, "elapsed: %lf\n", t.elapsed());

		phelm_shutdown();
	} catch (const std::exception & e) {
		fprintf(stderr, "exception: %s\n", e.what());
	}
	return (int)(!result);
}

