/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky
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

#include "linal_cuda.h"
#include "ver.h"

VERSION("$Id$");

void vector_splay (int n, int threads_min, int threads_max, 
	int grid_width, int *ctas, 
	int *elems_per_cta, int * threads_per_cta)
{
	if (n < threads_min) {
		*ctas              = 1;
		*elems_per_cta     = n;
		*threads_per_cta   = threads_min;
	} else if (n < (grid_width * threads_min)) {
		*ctas              = ((n + threads_min - 1) / threads_min);
		*threads_per_cta   = threads_min;
		*elems_per_cta     = *threads_per_cta;
	} else if (n < (grid_width * threads_max)) {
		int grp;
		*ctas              = grid_width;
		grp                = ((n + threads_min - 1) / threads_min);
		*threads_per_cta   = (((grp + grid_width -1) / grid_width) * threads_min);
		*elems_per_cta     = *threads_per_cta;
	} else {
		int grp;
		*ctas              = grid_width;
		*threads_per_cta   = threads_max;
		grp                = ((n + threads_min - 1) / threads_min);
		grp                = ((grp + grid_width - 1) / grid_width);
		*elems_per_cta     = grp * threads_min;
	}
}

namespace phelm {

template < typename T >
__global__ void sparse_mult_vector_l_(T * r, 
	const int * Ap, 
	const int * Ai, 
	const T * Ax,
	const T * x, 
	int n)
{
	int threads = gridDim.x  * blockDim.x;
	//int j       = blockDim.x * blockIdx.x + threadIdx.x;

	int j;
	//int num_per_block  = n / gridDim.x + (n % gridDim.x > 0);
	//int num   = num_per_block / blockDim.x + (num_per_block % blockDim.x > 0);
	int start = blockDim.x * blockIdx.x + threadIdx.x;
	//int fin   = start + num;

	for (j = start; j < n; j += threads) {
		const T *p = &Ax[Ap[j]];
		int i0;
		r[j] = 0;

		for (i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++p) {
			int i = Ai[i0];
			r[j] += *p * x[i];
		}
	}
}

__host__ void sparse_mult_vector_ld(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	sparse_mult_vector_l_ <<< ctas, threads >>> (r, Ap, Ai, Ax, x, n);
}

__host__ void sparse_mult_vector_lf(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n)
{
	int blocks, threads, elems;
	vector_splay (n, 32, 128, 80, &blocks, &elems, &threads);
	//threads = 4; blocks  = 4;//n / threads + ((n % threads) == 0?0:1);
	//fprintf(stderr, "%d, %d: %d\n", n, blocks, threads);
	sparse_mult_vector_l_ <<< blocks, threads >>> (r, Ap, Ai, Ax, x, n);
}

/* r = k1 * a + k2 * b */
template < typename T >
__global__ void vec_sum1_(T * r, const T * a, const T *b, T k1, T k2, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < n; i += threads) {
		r[i] = k1 * a[i] + k2 * b[i];
	}
}

__host__ void vec_sum1(float * r, const float * a, const float *b, float k1, float k2, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_sum1_ <<< ctas, threads >>> (r, a, b, k1, k2, n);
}

__host__ void vec_sum1(double * r, const double * a, const double *b, double k1, double k2, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_sum1_ <<< ctas, threads >>> (r, a, b, k1, k2, n);
}

/* r = a + k2 * b */
template < typename T >
__global__ void vec_sum2_(T * r, const T * a, const T *b, T k2, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < n; i += threads) {
		r[i] = a[i] + k2 * b[i];
	}
}

__host__ void vec_sum2(float * r, const float * a, const float *b, float k2, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_sum2_ <<< ctas, threads >>> (r, a, b, k2, n);
}

__host__ void vec_sum2(double * r, const double * a, const double *b, double k2, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_sum2_ <<< ctas, threads >>> (r, a, b, k2, n);
}

/* r = a + b */
template < typename T >
__global__ void vec_sum_(T * r, const T * a, const T * b, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;;
	for (;i < n; i += threads) {
		r[i] = a[i] + b[i];
	}
}

__host__ void vec_sum(float * r, const float * a, const float *b, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_sum_ <<< ctas, threads >>> (r, a, b, n);
}

__host__ void vec_sum(double * r, const double * a, const double *b, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_sum_ <<< ctas, threads >>> (r, a, b, n);
}

/* r = a * b */
template < typename T >
__global__ void vec_mult_(T * r, const T * a, const T * b, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;;
	for (;i < n; i += threads) {
		r[i] = a[i] * b[i];
	}
}

__host__ void vec_mult(float * r, const float * a, const float *b, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_mult_ <<< ctas, threads >>> (r, a, b, n);
}

__host__ void vec_mult(double * r, const double * a, const double *b, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_mult_ <<< ctas, threads >>> (r, a, b, n);
}

/* r = a - b*/
template < typename T >
__global__ void vec_diff_(T * r, const T * a, const T *b, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < n; i += threads) {
		r[i] = a[i] - b[i];
	}
}

__host__ void vec_diff(float * r, const float * a, const float *b, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_diff_ <<< ctas, threads >>> (r, a, b, n);
}

__host__ void vec_diff(double * r, const double * a, const double *b,  int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_diff_ <<< ctas, threads >>> (r, a, b, n);
}

/* r = b * k*/
template < typename T >
__global__ void vec_mult_scalar_(T * r, const T * b, T k, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < n; i += threads) {
		r[i] = k * b[i];
	}
}

__host__ void vec_mult_scalar(float * r, const float * b, float k, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_mult_scalar_ <<< ctas, threads >>> (r, b, k, n);
}

__host__ void vec_mult_scalar(double * r, const double * b, double k, int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	vec_mult_scalar_ <<< ctas, threads >>> (r, b, k, n);
}

}

