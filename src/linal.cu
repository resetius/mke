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
#include "shmem.h"
#include "ver.h"

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

texture < float > texX;
texture < float > texAx;
texture < int > texAp;
texture < int > texAi;

texture < float > texA;
texture < float > texB;

namespace phelm {

template < typename T >
__global__ void sparse_mult_vector_l_(T * r, 
	const int * Ap, 
	const int * Ai, 
	const T * Ax,
	const T * x, 
	size_t off_x,
	size_t off_ap,
	size_t off_ai,
	size_t off_ax,
	int n)
{
	int threads = gridDim.x  * blockDim.x;
	int j;
	int start = blockDim.x * blockIdx.x + threadIdx.x;

	for (j = start; j < n; j += threads) {
		int i0 = tex1Dfetch(texAp, off_ap + j);
		int to = tex1Dfetch(texAp, off_ap + j + 1);
		
		//const T *p = &Ax[i0];
		T rj = (T)0.0;

		//for (i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++p) {
		for (; i0 < to; ++i0/*, ++p*/) {
			int i = Ai[i0];
			//int i = tex1Dfetch(texAi, off_ai + i0);
			//rj   += *p * tex1Dfetch(texX, off_x + i);//x[i];
			rj   += tex1Dfetch(texAx, off_ax + i0) * tex1Dfetch(texX, off_x + i);
		}

		r[j] = rj;
	}
}

__host__ void sparse_mult_vector(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n,
	int nz)
{
	SPLAY(n);
	sparse_mult_vector_l_ <<< blocks, threads >>> (r, Ap, Ai, Ax, x, 1, 1, 1, 1, n);
}

__host__ void sparse_mult_vector(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n,
	int nz)
{
	SPLAY(n);

	size_t x_off;
	size_t Ap_off;
	size_t Ai_off;
	size_t Ax_off;

	if (cudaBindTexture (&x_off,  texX, x, n * sizeof(float)) != cudaSuccess) {
		exit(1);
	}
	if (cudaBindTexture (&Ax_off,  texAx, Ax, nz * sizeof(float)) != cudaSuccess) {
		exit(1);
	}
	if (cudaBindTexture (&Ap_off, texAp, Ap, (n + 1) * sizeof(int)) != cudaSuccess) {
		exit(1);
	}
	if (cudaBindTexture (&Ai_off, texAi, Ai, n * sizeof(int)) != cudaSuccess) {
		exit(1);
	}

	Ax_off /= sizeof(float);
	x_off /= sizeof(float);
	Ap_off /= sizeof(int);
	Ai_off /= sizeof(int);

	sparse_mult_vector_l_ <<< blocks, threads >>> (r, Ap, Ai, Ax, x, x_off, Ap_off, Ai_off, Ax_off, n);
	cudaUnbindTexture (texX);
	cudaUnbindTexture (texAx);
	cudaUnbindTexture (texAp);
	cudaUnbindTexture (texAi);
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
template < typename T >
__global__ void vec_sum2_(T * r, const T * a, const T *b, T k2, 
	size_t off_a, size_t off_b, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < n; i += threads) {
		r[i] = tex1Dfetch(texA, off_a + i) + k2 * tex1Dfetch(texB, off_b + i);//a[i] + k2 * b[i];
	}
}

__host__ void vec_sum2(float * r, const float * a, const float *b, float k2, int n)
{
	int blocks, threads, elems;
	vector_splay (n, 32, 128, 80, &blocks, &elems, &threads);
	size_t a_off, b_off;

	if (cudaBindTexture (&a_off,  texA, a, n * sizeof(float)) != cudaSuccess) {
		exit(1);
	}

	if (cudaBindTexture (&b_off,  texB, b, n * sizeof(float)) != cudaSuccess) {
		exit(1);
	}

	vec_sum2_ <<< blocks, threads >>> (r, a, b, k2, a_off/sizeof(float), b_off/sizeof(float), n);

	cudaUnbindTexture(texA);
	cudaUnbindTexture(texB);
}

__host__ void vec_sum2(double * r, const double * a, const double *b, double k2, int n)
{
	SPLAY(n);
	vec_sum2_ <<< blocks, threads >>> (r, a, b, k2, 1, 1, n);
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
	int blocks, threads, elems;
	vector_splay (n, 32, 128, 80, &blocks, &elems, &threads);
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
	SPLAY(n);
	vec_mult_scalar_ <<< blocks, threads >>> (r, b, k, n);
}

__host__ void vec_mult_scalar(double * r, const double * b, double k, int n)
{
	SPLAY(n);
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

template <class T>
__global__ void
reduce (T *g_idata, T *g_odata, unsigned int n)
{
#if 0
    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata[tid] = (i < n) ? g_idata[i] : 0;
    
    __syncthreads();

    // do reduction in shared mem
	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
	if (tid < s) {
		sdata[tid] += sdata[tid + s];
	}
	__syncthreads();
	}

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
#endif

	SharedMemory < T > shmem;
	T * sdata = shmem.getPointer();

	unsigned int tid = threadIdx.x;
	int blockSize = blockDim.x;
	unsigned int i = blockIdx.x*(blockSize*2) + tid;
	unsigned int gridSize = blockSize*2*gridDim.x;
	sdata[tid] = 0;
	
	while (i < n) { sdata[tid] += g_idata[i] + g_idata[i+blockSize]; i += gridSize; }
	__syncthreads();
	if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
	if (tid < 32) {
		if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
		if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
		if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
		if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
		if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
		if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
	}
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

unsigned int nextPow2( unsigned int x ) 
{
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

template < typename T >
__host__ T vec_scalar2_(const T * a, const T * b, int n)
{
	int threads = 128;
	int blocks  = (n + threads - 1) / threads;

	int maxThreads = 128;
	int maxBlocks  = blocks;

	T * v1;
	T * v2;
	T * final;
	T answer = (T)0.0;;

	cudaMalloc((void**)&v1, n * sizeof(T));
	cudaMalloc((void**)&v2, n * sizeof(T));

	vec_mult(v1, a, b, n);

	int smemsize = threads * sizeof(float);

	reduce <<< blocks, threads, smemsize >>> (v1, v2, n);

	int N = blocks;
	int final_threshold = 1000;

	while (N > final_threshold) {
        threads = (N < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
        blocks  = (N + (threads * 2 - 1)) / (threads * 2);
		blocks  = min(maxBlocks, blocks);

		reduce <<< blocks, threads, smemsize >>> (v2, v2, N);

		N = (N + (threads*2-1)) / (threads*2);
	}

	final = (T*)malloc(N * sizeof(T));
	cudaMemcpy(final, v2, N * sizeof(T), cudaMemcpyDeviceToHost);
	for (int i = 0; i < N; ++i) 
	{
		answer += final[i];
	}

	cudaFree(v1);
	cudaFree(v2);
	free(final);
	return answer;
}

/*
__host__ double vec_scalar2(const double * a, const double * b, int n)
{
	return vec_scalar2_(a, b, n);
}

__host__ float vec_scalar2(const float * a, const float * b, int n)
{
	return vec_scalar2_(a, b, n);
}
*/

}
