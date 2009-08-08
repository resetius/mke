/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (Алексей Озерицкий)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Redistributions in any form must be accompanied by information on
 *    how to obtain complete source code for the Phelm software and any
 *    accompanying software that uses the Phelm software.  The source code
 *    must either be included in the distribution or be available for no
 *    more than the cost of distribution plus a nominal fee, and must be
 *    freely redistributable under reasonable conditions.  For an
 *    executable file, complete source code means the source code for all
 *    modules it contains.  It does not include source code for modules or
 *    files that typically accompany the major components of the operating
 *    system on which the executable file runs.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "util.h"
#include "ver.h"

VERSION("$Id$");

namespace phelm {
/**
 * Gauss
 * Copyright (c) 2009 Andrey Kornev, Andrey Ivanchikov, Alexey Ozeritsky
 * from manipro 
 */
static void gauss_reverse (const double *A, const double *b, double *x, int n, int diagUnit)
{
	int j, k;

	for (k = n - 1; k >= 0; k--) {
		x[k] = b[k];
		for (j = k + 1; j < n; j++) {
			x[k] = x[k] - x[j] * A[k*n+j];
		}
		if (!diagUnit) {
			x[k] = x[k] / A[k*n+k];
		}
	}
}

int gauss (double *A, double *b, double *x, int n)
{
	int i, j, k;
	double p;
	int imax;
	double Eps = 1.e-15;

	for (k = 0; k < n; k++) {
		imax = k;

		for (i = k + 1; i < n; i++) {
			if (fabs(A[i*n+k]) > fabs(A[imax*n+k])) imax = i;
		}

		for (j = k; j < n; j++) {
			p = A[imax*n+j];
			A[imax*n+j] = A[k*n+j];
			A[k*n+j] = p;
		}
		p = b[imax];
		b[imax] = b[k];
		b[k] = p;

		p = A[k*n+k];

		if (fabs(p) < Eps) {
			printf("Warning in %s %s : Near-null zero element\n", __FILE__, __FUNCTION__);
			return -1;
		}

		for (j = k; j < n; j++) {
			A[k*n+j] = A[k*n+j] / p;
		}
		b[k] = b[k] / p;

		for (i = k + 1; i < n; i++) {
			p = A[i*n+k];
			for (j = k; j < n; j++) {
				A[i*n+j] = A[i*n+j] - A[k*n+j] * p;
			}
			b[i] = b[i] - b[k] * p;
		}
	}

	gauss_reverse(A, b, x, n, true);

	return 0;
}

template < typename T >
void mat_print_(const T * A, int n)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%9.2le ", A[i * n + j]);
		}
		printf("\n");
	}
}

void mat_print(const float * A, int n)
{
	mat_print_(A, n);
}

void mat_print(const double * A, int n)
{
	mat_print_(A, n);
}

template < typename T >
void vec_print_(const T * A, int n)
{
	for (int i = 0; i < n; ++i) {
		printf("%9.2le ", A[i]);
	}
	printf("\n");
}

void vec_print(const double * A, int n)
{
	vec_print_(A, n);
}

void vec_print(const float * A, int n)
{
	vec_print_(A, n);
}

template < typename T >
void mat_mult_vector_(T * r, const T * A, const T * x, int n)
{
	int i, j;
	for (i = 0; i < n; ++i) {
		r[i] = 0.0;
		for (j = 0; j < n; ++j) {
			r[i] += A[i * n + j] * x[j];
		}
	}
}

void mat_mult_vector(double * r, const double * A, const double * x, int n)
{
	mat_mult_vector_(r, A, x, n);
}

void mat_mult_vector(float * r, const float * A, const float * x, int n)
{
	mat_mult_vector_(r, A, x, n);
}

template < typename T, typename Sparse >
void sparse_mult_vector_r_(T * r, const Sparse * A, const T * x, int n)
{
	int i0, i, j;
	const T * p = A->Ax;

	memset(r, 0, n * sizeof(T));
	for (j = 0; j < n; ++j) {
		for (i0 = A->Ap[j]; i0 < A->Ap[j + 1]; ++i0, ++p) {
			i = A->Ai[i0];
			r[i] += *p * x[j];
		}
	}
}

void sparse_mult_vector_r(double * r, const Sparse * A, const double * x, int n)
{
	sparse_mult_vector_r_(r, A, x, n);
}

void sparse_mult_vector_r(float * r, const Sparsef * A, const float * x, int n)
{
	sparse_mult_vector_r_(r, A, x, n);
}

template < typename T >
void vec_mult_(T * r, const T * a, const T *b, int n)
{
	int i;
//#pragma omp parallel for
	for (i = 0; i < n; ++i) {
		r[i] = a[i] * b[i];
	}
}

void vec_mult(double * r, const double * a, const double *b, int n)
{
	vec_mult_(r, a, b, n);
}

void vec_mult(float * r, const float * a, const float *b, int n)
{
	vec_mult_(r, a, b, n);
}

template < typename Sparse >
void sparse_print_(const Sparse * A, int n, FILE * f)
{
	int i, i0, j, k, i_old;
	const typename Sparse::data_type * p = A->Ax;
	for (j = 0; j < n; ++j) {
		i_old = -1;
		for (i0 = A->Ap[j]; i0 < A->Ap[j + 1]; ++i0, ++p) {
			i = A->Ai[i0];
			for (k = i_old; k < i - 1; ++k) {
				fprintf(f, "%8.3lf ", 0.0);
			}
			fprintf(f, "%8.3lf ", (double)*p);
			i_old = i;
		}

		for (k = i_old + 1; k < n; ++k) {
			fprintf(f, "%8.3lf ", 0.0);
		}
		fprintf(f, "\n");
	}
}

void sparse_print(const Sparse * A, int n, FILE * f)
{
	sparse_print_ (A, n, f);
}

void sparse_print(const Sparsef * A, int n, FILE * f)
{
	sparse_print_ (A, n, f);
}

} /* namespace */