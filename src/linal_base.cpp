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

VERSION ("$Id$");

namespace phelm
{

void set_num_threads (int threads)
{
#ifdef _OPENMP
	if (threads)
	{
		omp_set_num_threads (threads);
	}
#endif
}

/**
 * Gauss
 * Copyright (c) 2009 Andrey Kornev, Andrey Ivanchikov, Alexey Ozeritsky
 * from manipro
 */
static void gauss_reverse (const double *A, const double *b, double *x, int n, int diagUnit)
{
	int j, k;

	for (k = n - 1; k >= 0; k--)
	{
		x[k] = b[k];
		for (j = k + 1; j < n; j++)
		{
			x[k] = x[k] - x[j] * A[k*n+j];
		}
		if (!diagUnit)
		{
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

	for (k = 0; k < n; k++)
	{
		imax = k;

		for (i = k + 1; i < n; i++)
		{
			if (fabs (A[i*n+k]) > fabs (A[imax*n+k]) ) imax = i;
		}

		for (j = k; j < n; j++)
		{
			p = A[imax*n+j];
			A[imax*n+j] = A[k*n+j];
			A[k*n+j] = p;
		}
		p = b[imax];
		b[imax] = b[k];
		b[k] = p;

		p = A[k*n+k];

		if (fabs (p) < Eps)
		{
			printf ("Warning in %s %s : Near-null zero element\n", __FILE__, __FUNCTION__);
			return -1;
		}

		for (j = k; j < n; j++)
		{
			A[k*n+j] = A[k*n+j] / p;
		}
		b[k] = b[k] / p;

		for (i = k + 1; i < n; i++)
		{
			p = A[i*n+k];
			for (j = k; j < n; j++)
			{
				A[i*n+j] = A[i*n+j] - A[k*n+j] * p;
			}
			b[i] = b[i] - b[k] * p;
		}
	}

	gauss_reverse (A, b, x, n, true);

	return 0;
}

template < typename T >
void mat_print_ (const T * A, int n)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			printf ("%9.2le ", A[i * n + j]);
		}
		printf ("\n");
	}
}

void mat_print (const float * A, int n)
{
	mat_print_ (A, n);
}

void mat_print (const double * A, int n)
{
	mat_print_ (A, n);
}

double vec_norm2 (const double * v, int n)
{
	return sqrt (vec_scalar2 (v, v, n) );
}

float vec_norm2 (const float * v, int n)
{
	float f = vec_scalar2 (v, v, n);
	return sqrtf (f);
}

template < typename T >
void vec_print_ (const T * A, int n)
{
	for (int i = 0; i < n; ++i)
	{
		printf ("%9.2le ", A[i]);
	}
	printf ("\n");
}

void vec_print (const double * A, int n)
{
	vec_print_ (A, n);
}

void vec_print (const float * A, int n)
{
	vec_print_ (A, n);
}

template < typename T >
void mat_mult_vector_ (T * r, const T * A, const T * x, int n)
{
#pragma omp parallel
	{
		int block_dim = 400; //cache size = (block_dim * block_dim * 8)
		int blocks = (n + block_dim - 1) / block_dim;

#pragma omp for
		for (int i = 0; i < n; ++i)
		{
			r[i] = 0;
		}

		for (int l = 0; l < blocks; ++l )
		{
			int fl = n * l;
			fl /= blocks;
			int ll = n * (l + 1);
			ll = ll / blocks - 1;

			for (int m = 0; m < blocks; ++m)
			{
				int fm = n * m;
				fm /= blocks;
				int lm = n * (m + 1);
				lm = lm / blocks - 1;

				// blocks:
				// R[fl] += A[fl, fm] * X[fl]

#pragma omp for
				for (int i = fl; i <= ll; ++i)
				{
					const T * ax = &A[i * n + fm];
					const T * xx = &x[fm];

					T s = 0.0;
					for (int j = fm; j <= lm; ++j)
					{
						s += *ax++ * *xx++;
					}
					r[i] += s;
				}
			}
		}
	}
}

void mat_mult_vector (double * r, const double * A, const double * x, int n)
{
	mat_mult_vector_ (r, A, x, n);
}

void mat_mult_vector (float * r, const float * A, const float * x, int n)
{
	mat_mult_vector_ (r, A, x, n);
}

template < typename T >
void mat_mult_vector_stupid_ (T * r, const T * A, const T * x, int n)
{
#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		T s = 0.0;
		for (int j = 0; j < n; ++j)
		{
			s += A[i * n + j] * x[j];
		}
		r[i] = s;
	}
}

void mat_mult_vector_stupid (double * r, const double * A, const double * x, int n)
{
	mat_mult_vector_stupid_ (r, A, x, n);
}

void mat_mult_vector_stupid (float * r, const float * A, const float * x, int n)
{
	mat_mult_vector_stupid_ (r, A, x, n);
}

template < typename T >
void csr_print_ (const int * Ap, const int * Ai,
                 const T * Ax, int n, FILE * f)
{
	int i, i0, j, k, i_old;
	const T * p = Ax;
	for (j = 0; j < n; ++j)
	{
		i_old = -1;
		for (i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++p)
		{
			i = Ai[i0];
			for (k = i_old; k < i - 1; ++k)
			{
				fprintf (f, "%8.3lf ", 0.0);
			}
			fprintf (f, "%8.3lf ", (double) *p);
			i_old = i;
		}

		for (k = i_old + 1; k < n; ++k)
		{
			fprintf (f, "%8.3lf ", 0.0);
		}
		fprintf (f, "\n");
	}
}

void sparse_print (const int * Ap, const int * Ai,
                   const double * Ax, int n, FILE * f)
{
	csr_print_ (Ap, Ai, Ax, n, f);
}

void sparse_print (const int * Ap, const int * Ai,
                   const float * Ax, int n, FILE * f)
{
	csr_print_ (Ap, Ai, Ax, n, f);
}

} /* namespace */
