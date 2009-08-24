/* -*- charset: utf-8 -*- */

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "gmres.h"
#include "util.h"
#include "solver.h"

using namespace phelm;
using namespace std;

bool check(float val)
{
	return fabs(val) < 1e-5;
}

bool check(double val)
{
	return fabs(val) < 1e-12;
}

template < typename T >
bool test_solver()
{
	int i, j = 0;
	int n  = 320000;
//	int n  = 500000;
//	int n  = 50000;
//	int n = 20;

	/**
	 * 4 2 0 .... 0
	 * 1 4 2 .... 0
	 * 0 1 4 2 .. 0
	 * 0 0 1 4 .. 0
	 * ............
	 * ...... 1 4 2
	 * .........1 4
	 */

	bool ret = true;

	vector < T > b(n);
	vector < T > x1(n);
	vector < T > x2(n);
	vector < T > x3(n);
	vector < T > v(n);

	/* матрицу записываем по строкам! */
	SparseSolver  < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > > M1(n);
#ifdef UMFPACK
	UmfPackSolver < T, StoreCSR < T , Allocator > > M2(n);
#endif

#ifdef SUPERLU
	SuperLUSolver < T, StoreCSR < T , Allocator > > M3(n);
#endif

	for (i = 0; i < n; ++i) {
		if (i == 0) {
			/* 4 2 0 .... 0 */
			M1.add(i, i,     4);
			M1.add(i, i + 1, 2);
#ifdef UMFPACK
			M2.add(i, i,     4);
			M2.add(i, i + 1, 2);
#endif
#ifdef SUPERLU
			M3.add(i, i,     4);
			M3.add(i, i + 1, 2);
#endif
		} else if (i < n - 1) {
			/* 1 4 2 .... 0 */
			M1.add(i, i - 1, 1);
			M1.add(i, i,     4);
			M1.add(i, i + 1, 2);
#ifdef UMFPACK
			M2.add(i, i - 1, 1);
			M2.add(i, i,     4);
			M2.add(i, i + 1, 2);
#endif
#ifdef SUPERLU
			M3.add(i, i - 1, 1);
			M3.add(i, i,     4);
			M3.add(i, i + 1, 2);
#endif
		} else {
			/* 0 .... 1 4 */
			M1.add(i, i - 1, 1);
			M1.add(i, i,     4);
#ifdef UMFPACK
			M2.add(i, i - 1, 1);
			M2.add(i, i,     4);
#endif
#ifdef SUPERLU
			M3.add(i, i - 1, 1);
			M3.add(i, i,     4);
#endif
		}

		b[i] = 1;
	}

	Timer t;

	t.restart();
	for (int k = 0; k < 1; ++k) {
		M1.solve(&x1[0], &b[0]);
	}
	fprintf(stderr, "gmres solve: %lf\n", t.elapsed());
#ifdef UMFPACK
	t.restart();
	for (int k = 0; k < 1; ++k) {
		M2.solve(&x2[0], &b[0]);
	}
	fprintf(stderr, "umfpack solve: %lf\n", t.elapsed());
#endif
#ifdef SUPERLU
	t.restart();
	for (int k = 0; k < 1; ++k) {
		M3.solve(&x3[0], &b[0]);
	}
	fprintf(stderr, "superlu solve: %lf\n", t.elapsed());
#endif
#if defined(UMFPACK) && defined(SUPERLU)
	T nr;
	vec_diff(&v[0], &x2[0], &x3[0], (int)x1.size());
	nr = vec_norm2(&v[0], (int)x2.size());
	fprintf(stderr, "%.16le\n", (double)nr);
	ret &= check(nr);
#endif
	vec_diff(&v[0], &x1[0], &x2[0], (int)x1.size());
	nr = vec_norm2(&v[0], (int)x2.size());
	fprintf(stderr, "%.16le\n", (double)nr);
//	ret &= check(nr);

	return ret;
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

			t.restart(); result &= test_solver < double > ();
			fprintf(stderr, "test_solver < double > (): %lf, %d\n", t.elapsed(), (int)result);

			fprintf(stderr, "testing float:\n");

			t.restart(); result &= test_solver < float > ();
			fprintf(stderr, "test_solver < float > (): %lf, %d\n", t.elapsed(), (int)result);

		} else {
			;
		}
		fprintf(stderr, "elapsed: %lf\n", t.elapsed());

		phelm_shutdown();
	} catch (const std::exception & e) {
		fprintf(stderr, "exception: %s\n", e.what());
	}
	return (int)(!result);
}

