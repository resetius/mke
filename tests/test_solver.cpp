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
	return fabs(val) < 1e-6;
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
	int nz = n + n - 1 + n - 1;

	vector < int > Ap((n + 1));
	vector < int > Ai(nz);
	vector < T > Ax(nz);
	vector < T > b(n);
	vector < T > x(n);
	vector < T > y(n);

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

	Timer t;

	t.restart();
	UmfPackMatrix < T > lu2(&Ap[0], &Ai[0], &Ax[0], n, nz);
	for (int k = 0; k < 1000; ++k) {
		lu2.solve(&x[0], &b[0]);
	}
	fprintf(stderr, "umfpack solve: %lf\n", t.elapsed());

	t.restart();
	SuperLUMatrix < T > lu1(&Ap[0], &Ai[0], &Ax[0], n, nz);
	for (int k = 0; k < 1000; ++k) {
		lu1.solve(&y[0], &b[0]);
	}
	fprintf(stderr, "superlu solve: %lf\n", t.elapsed());

	vec_diff(&y[0], &y[0], &x[0], (int)x.size());
	T nr = vec_norm2(&y[0], (int)y.size());
	fprintf(stderr, "%lf\n", (double)nr);
	return check(nr);
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

