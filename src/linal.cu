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

#include "linal_cuda.h"

static void vector_splay (int n, int threads_min, int threads_max, 
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

__global__ void _sparse_mult_vector_ld(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n)
{
	int j;

	for (j = 0; j < n; ++j) {
		const double *p = &Ax[Ap[j]];
		int i0;
		r[j] = 0;

		for (i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++p) {
			int i = Ai[i0];
			r[j] += *p * x[i];
		}
	}
}

__host__ void API sparse_mult_vector_ld(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	_sparse_mult_vector_ld <<< ctas, threads >>> (r, Ap, Ai, Ax, x, n);
}

__global__ void _sparse_mult_vector_lf(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n)
{
	int j;
	int tid     = threadIdx.x;
	int threads = gridDim.x  * blockDim.x;
	int start   = blockDim.x * blockIdx.x;

	for (j = start + tid; j < n; j += threads) {
		const float *p = &Ax[Ap[j]];
		int i0;
		r[j] = 0;

		for (i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++p) {
			int i = Ai[i0];
			r[j] += *p * x[i];
		}
	}
}

__host__ void API sparse_mult_vector_lf(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n)
{
	int ctas, threads, elems;
	vector_splay (n, 32, 128, 80, &ctas, &elems, &threads);
	_sparse_mult_vector_lf <<< ctas, threads >>> (r, Ap, Ai, Ax, x, n);
}
