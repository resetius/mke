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

#include <stdio.h>

#ifndef WIN32
#include <stdint.h>
#endif

#include "linal_cuda.h"
#include "shmem.h"
#include "texture.h"
#include "ver.h"
#include "reduction.h"

VERSION("$Id$");

void vector_splay (int n, int threads_min, int threads_max, 
	int grid_width, int *blocks, 
	int *elems_per_block, int * threads_per_block)
{
	if (n < threads_min) {
		*blocks            = 1;
		*elems_per_block   = n;
		*threads_per_block = threads_min;
	} else if (n < (grid_width * threads_min)) {
		*blocks            = ((n + threads_min - 1) / threads_min);
		*threads_per_block = threads_min;
		*elems_per_block   = *threads_per_block;
	} else if (n < (grid_width * threads_max)) {
		int grp;
		*blocks            = grid_width;
		grp                = ((n + threads_min - 1) / threads_min);
		*threads_per_block = (((grp + grid_width -1) / grid_width) * threads_min);
		*elems_per_block   = *threads_per_block;
	} else {
		int grp;
		*blocks            = grid_width;
		*threads_per_block = threads_max;
		grp                = ((n + threads_min - 1) / threads_min);
		grp                = ((grp + grid_width - 1) / grid_width);
		*elems_per_block   = grp * threads_min;
	}
}

register_texture(float, texX1);
register_texture(float, texAX);
register_texture(int, texAP);
register_texture(int, texAI);

register_texture(float, texA);
register_texture(float, texB);

namespace phelm {

template < typename T, typename APR, typename AIR, typename XR, typename AXR >
__global__ void sparse_mult_vector_csr_(T * r, 
	APR Ap, 
	AIR Ai, 
	AXR Ax,
	XR x, 
	int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i, i0, to, j;
	int start = blockDim.x * blockIdx.x + threadIdx.x;

	for (j = start; j < n; j += threads) {
		i0 = Ap.get(j);    
		to = Ap.get(j + 1);

		T rj = (T)0.0;

		for (; i0 < to; ++i0) {
			i   = Ai.get(i0);
			rj += Ax.get(i0) * x.get(i);
		}

		r[j] = rj;
	}
}

__host__ void csr_mult_vector_r(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n,
	int nz)
{
	SPLAY(n);
	simple_reader < double > XR(x);
	simple_reader < double > AXR(Ax);
	simple_reader < int > AIR(Ap);
	simple_reader < int > APR(Ap);

	sparse_mult_vector_csr_ <<< blocks, threads >>> (r, APR, AIR, AXR, XR, n);
}

__host__ void csr_mult_vector_r(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n,
	int nz)
{
	SPLAY(n);

	bool useTexture;

	useTexture = ((n + 1 < MAX_1DBUF_SIZE) && (nz < MAX_1DBUF_SIZE));

	if (n < 1000) /* experimental bound */
	{
		useTexture = false;
	}

	if (useTexture) {
		texture_reader(texX1) XR(x, n);
		texture_reader(texAX) AXR(Ax, nz);
		texture_reader(texAI) AIR(Ai, nz);
		texture_reader(texAP) APR(Ap, n + 1);

		sparse_mult_vector_csr_ <<< blocks, threads >>> (r, APR, AIR, AXR, XR, n);
	} else {
		simple_reader < float > XR(x);
		simple_reader < float > AXR(Ax);
		simple_reader < int > AIR(Ai);
		simple_reader < int > APR(Ap);

		sparse_mult_vector_csr_ <<< blocks, threads >>> (r, APR, AIR, AXR, XR, n);
	}
}

template < typename T, typename AIR, typename AXR, typename XR >
__global__ void ell_mult(
			   T * r, 
			   AIR Ai, 
			   AXR Ax,
			   XR x, 
			   int n,
			   int cols, 
			   int stride)
{
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	if (row < n) {
		T sum = 0;

		for (int i0 = 0; i0 < cols; i0++){
			const T A_ij = Ax.get(stride * i0 + row);

			if (A_ij != 0) {
				const int col = Ai.get(stride * i0 + row);
				sum += A_ij * x.get(col);
			}
		}
	    r[row] = sum;
	}
}

__host__ void 
ell_mult_vector_r(float * r, const int * Ai, const float * Ax, 
	const float * x, int n, int cols, int stride)
{
	SPLAY2(n);

	texture_reader(texX1) XR(x, n);
	texture_reader(texAX) AXR(Ax, cols * stride);
	texture_reader(texAI) AIR(Ai, cols * stride);
	
	ell_mult<<<blocks, threads>>>(r, AIR, AXR, XR, n, cols, stride);
}

__host__ void 
ell_mult_vector_r(double * r, const int * Ai, const double * Ax, 
	const double * x, int n, int cols, int stride)
{
	SPLAY2(n);

	simple_reader < double > XR(x);
	simple_reader < double > AXR(Ax);
	simple_reader < int > AIR(Ai);
	
	ell_mult<<<blocks, threads>>>(r, AIR, AXR, XR, n, cols, stride);
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
	SPLAY(n);
	vec_sum1_ <<< blocks, threads >>> (r, a, b, k1, k2, n);
}

__host__ void vec_sum1(double * r, const double * a, const double *b, double k1, double k2, int n)
{
	SPLAY(n);
	vec_sum1_ <<< blocks, threads >>> (r, a, b, k1, k2, n);
}

/* r = a + k2 * b */
template < typename T, typename AR, typename BR >
__global__ void vec_sum2_(T * r, AR a, BR b, T k2, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < n; i += threads) {
		r[i] = a.get(i) + k2 * b.get(i);//a[i] + k2 * b[i];
	}
}

__host__ void vec_sum2(float * r, const float * a, const float *b, float k2, int n)
{
	SPLAY2(n);

	bool useTexture;

	useTexture = (n < MAX_1DBUF_SIZE);

	if ((n < 10000) ||
		((!(((uintptr_t) a) % WORD_ALIGN)) && 
		(!(((uintptr_t) b) % WORD_ALIGN)))) 
	{
		useTexture = false;
	}

	if (useTexture) {
		texture_reader(texA) AR(a, n);
		texture_reader(texB) BR(b, n);

		vec_sum2_ <<< blocks, threads >>> (r, AR, BR, k2, n);
	} else {
		simple_reader < float > AR(a);
		simple_reader < float > BR(b);

		vec_sum2_ <<< blocks, threads >>> (r, AR, BR, k2, n);
	}
}

__host__ void vec_sum2(double * r, const double * a, const double *b, double k2, int n)
{
	SPLAY2(n);

	simple_reader < double > AR(a);
	simple_reader < double > BR(b);

	vec_sum2_ <<< blocks, threads >>> (r, AR, BR, k2, n);
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
	SPLAY(n);
	vec_sum_ <<< blocks, threads >>> (r, a, b, n);
}

__host__ void vec_sum(double * r, const double * a, const double *b, int n)
{
	SPLAY(n);
	vec_sum_ <<< blocks, threads >>> (r, a, b, n);
}

/* r = a * b */
template < typename T >
__global__ void vec_mult_(T * r, const T * a, const T * b, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < n; i += threads) {
		r[i] = a[i] * b[i];
	}
}

__host__ void vec_mult(float * r, const float * a, const float *b, int n)
{
	SPLAY(n);
	vec_mult_ <<< blocks, threads >>> (r, a, b, n);
}

__host__ void vec_mult(double * r, const double * a, const double *b, int n)
{
	SPLAY(n);
	vec_mult_ <<< blocks, threads >>> (r, a, b, n);
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
	SPLAY(n);
	vec_diff_ <<< blocks, threads >>> (r, a, b, n);
}

__host__ void vec_diff(double * r, const double * a, const double *b,  int n)
{
	SPLAY(n);
	vec_diff_ <<< blocks, threads >>> (r, a, b, n);
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
	SPLAY2(n);
	vec_mult_scalar_ <<< blocks, threads >>> (r, b, k, n);
}

__host__ void vec_mult_scalar(double * r, const double * b, double k, int n)
{
	SPLAY2(n);
	vec_mult_scalar_ <<< blocks, threads >>> (r, b, k, n);
}

/*
template < typename T >
__global__ void reduction_(T * out, unsigned N, unsigned BlockStride)
{
	unsigned int i      = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int Stride = 2 * BlockStride;
	unsigned int j      = blockDim.x;
	while (j > 0)
	{
		if (Stride * i< N)
			out[Stride*i] += out[Stride*i+ (Stride>>1)];
		Stride <<= 1;
		j >>= 1;
		__syncthreads();
	}
}
*/

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

#include <malloc.h>
#ifdef WIN32
#define alloca _alloca
#endif

template < typename T, typename AR, typename BR >
struct Multiplier
{
	AR a_;
	BR b_;

	Multiplier(AR a, BR b): a_(a), b_(b) {}
	__device__ T get(int i) { 
		return a_.get(i) * b_.get(i);
	}
};

template < typename T >
__host__ T vec_scalar2_(const T * a, const T * b, int n)
{
	int maxThreads = 256;
	int maxBlocks  = 64;

	//int threads = maxThreads;
	//int blocks  = maxBlocks;

	int threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
	int blocks  = (n + (threads * 2 - 1)) / (threads * 2);
	blocks  = min(maxBlocks, blocks);

	T * v2 = 0;
	T answer = (T)0.0;;

	cudaMalloc((void**)&v2, blocks * sizeof(T));

	{
		texture_reader(texA) AR(a, n);
		texture_reader(texB) BR(b, n);
		Multiplier < T, texture_reader(texA), texture_reader(texB) > m(AR, BR);

		//simple_reader < T > AR(a);
		//simple_reader < T > BR(b);
		//Multiplier < T, simple_reader < T >, simple_reader < T > > m(AR, BR);

		reduce6 (threads, blocks, m, v2, n);
	}

	int N = blocks;
	int final_threshold = 1;

	texture_reader(texAX) VR(v2, blocks);
	//simple_reader < T > VR(v2);

	while (N > final_threshold) {
		threads = (N < maxThreads*2) ? nextPow2((N + 1)/ 2) : maxThreads;
		blocks  = (N + (threads * 2 - 1)) / (threads * 2);

		reduce5(threads, blocks, VR, v2, N);

		N = (N + (threads*2-1)) / (threads*2);
	}

	if (final_threshold > 1) {
		T * final = (T*)alloca(N * sizeof(T));
		cudaMemcpy(final, v2, N * sizeof(T), cudaMemcpyDeviceToHost);
		for (int i = 0; i < N; ++i) 
		{
			answer += final[i];
		}
	} else {
		cudaMemcpy(&answer, v2, N * sizeof(T), cudaMemcpyDeviceToHost);
	}

	cudaFree(v2);

	return answer;
}

/*
__host__ double vec_scalar2(const double * a, const double * b, int n)
{
	return vec_scalar2_(a, b, n);
}
*/

__host__ float vec_scalar2(const float * a, const float * b, int n)
{
	return vec_scalar2_(a, b, n);
}

}
