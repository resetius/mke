/* -*- charset: utf-8 -*- */

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "gmres.h"
#include "util.h"

using namespace phelm;
using namespace std;

template < typename T >
void test_gmres()
{
	int i, j = 0;
	int n  = 300000;
	int nz = n + n - 1 + n - 1;
	vector < int > Ap((n + 1));
	vector < int > Ai(nz);
	vector < T, phelm_allocator < T > > Ax(nz);
	vector < T, phelm_allocator < T > > b(n);
	vector < T, phelm_allocator < T > > x(n);
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

	A.Ax = &Ax[0];
	A.Ap = &Ap[0];
	A.Ai = &Ai[0];

	for (int k = 0; k < 10; ++k) {
		gmres(&x[0], &A, &b[0], 
			MatMultiplier < Sparse_t < T > > :: mult_vector_l, n, 10, 100);
	}
}

template < typename T >
void test_matvect()
{
	int i, j = 0;
	int n  = 5000000;
	int nz = n + n - 1 + n - 1;
	vector < int > Ap((n + 1));
	vector < int > Ai(nz);
	vector < T, phelm_allocator < T > > Ax(nz);
	vector < T, phelm_allocator < T > > b(n);
	vector < T, phelm_allocator < T > > x(n);
	Sparse A;

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

	A.Ax = &Ax[0];
	A.Ap = &Ap[0];
	A.Ai = &Ai[0];

	for (int k = 0; k < 1000; ++k) {
		sparse_mult_vector_l(&x[0], &A, &b[0], n);
	}
}

int main(int argc, char * argv[])
{
	test_gmres < double > ();
	test_gmres < float > ();
	//test_matvect();
	return 0;
}
