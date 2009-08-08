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

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <cublas.h>
#include <cuda_runtime_api.h>

#include "linal.h"
#include "ver.h"

VERSION("$Id$");

namespace phelm {

/**
 * r = k1 * a + k2 * b
 */
void vec_sum1(double * r, const double * a, const double *b, double k1, double k2, int n)
{
	cublasDcopy(n, a, 1, r, 1);
	cublasDscal(n, k1, r, 1);
	cublasDaxpy(n, k2, b, 1, r, 1);
}

void vec_sum1(float * r, const float * a, const float *b, float k1, float k2, int n)
{
	cublasScopy(n, a, 1, r, 1);
	cublasSscal(n, k1, r, 1);
	cublasSaxpy(n, k2, b, 1, r, 1);
}

void vec_sum(double * r, const double * a, const double *b, int n)
{
	cublasDcopy(n, a, 1, r, 1);
	cublasDaxpy(n, 1.0, b, 1, r, 1);
}

void vec_sum(float * r, const float * a, const float *b, int n)
{
	cublasScopy(n, a, 1, r, 1);
	cublasSaxpy(n, 1.0f, b, 1, r, 1);
}

void vec_sum2(double * r, const double * a, const double *b, double k2, int n)
{
	cublasDcopy(n, b, 1, r, 1);
	cublasDaxpy(n, k2, a, 1, r, 1);
}

void vec_sum2(float * r, const float * a, const float *b, float k2, int n)
{
	cublasScopy(n, b, 1, r, 1);
	cublasSaxpy(n, k2, a, 1, r, 1);
}

/**
 * a = b * k
 */
void vec_mult_scalar(double * a, const double * b, double k, int n)
{
	cublasDcopy(n, b, 1, a, 1);
	cublasDscal(n, k, a, 1);
}

void vec_mult_scalar(float * a, const float * b, float k, int n)
{
	cublasScopy(n, b, 1, a, 1);
	cublasSscal(n, k, a, 1);
}

void vec_copy(double * a, const double * b, int n)
{
	cublasDcopy(n, b, 1, a, 1);
}

void vec_copy(float * a, const float * b, int n)
{
	cublasScopy(n, b, 1, a, 1);
}

void vec_copy_from_host(double * a, const double * b, int n)
{
	cublasSetVector(n, sizeof(double), b, 1, a, 1);
}

void vec_copy_from_host(float * a, const float * b, int n)
{
	cublasSetVector(n, sizeof(float), b, 1, a, 1);
}

void vec_copy_from_host(int * a, const int * b, int n)
{
	cublasSetVector(n, sizeof(int), b, 1, a, 1);
}

void vec_copy_from_device(double * a, const double * b, int n)
{
	cublasGetVector(n, sizeof(double), b, 1, a, 1);
}

void vec_copy_from_device(float * a, const float * b, int n)
{
	cublasGetVector(n, sizeof(float), b, 1, a, 1);
}

void vec_copy_from_device(int * a, const int * b, int n)
{
	cublasGetVector(n, sizeof(int), b, 1, a, 1);
}

/**
 * r = a - b
 */
void vec_diff(double * r, const double * a, const double * b, int n)
{
	cublasDcopy(n, a, 1, r, 1);
	cublasDaxpy(n, -1.0, b, 1, r, 1);
}

void vec_diff(float * r, const float * a, const float * b, int n)
{
	cublasScopy(n, a, 1, r, 1);
	cublasSaxpy(n, -1.0f, b, 1, r, 1);
}

double vec_norm2(const double * v, int n)
{
	return cublasDnrm2(n, v, 1);
}

float vec_norm2(const float * v, int n)
{
	return cublasSnrm2(n, v, 1);
}

double vec_scalar2(const double * a, const double * b, int n)
{
	return cublasDdot(n, a, 1, b, 1);
}

float vec_scalar2(const float * a, const float * b, int n)
{
	return cublasSdot(n, a, 1, b, 1);
}

extern "C" {
#ifdef WIN32
#define API __stdcall
#else
#define API
#endif

void API sparse_mult_vector_ld(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n);

void API sparse_mult_vector_lf(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n);
}

void sparse_mult_vector_l(double * r, const Sparse * A, const double * x, int n)
{
	sparse_mult_vector_ld(r, A->Ap, A->Ai, A->Ax, x, n);
}

void sparse_mult_vector_l(float * r, const Sparsef * A, const float * x, int n)
{
	sparse_mult_vector_lf(r, A->Ap, A->Ai, A->Ax, x, n);
}

int check_device_supports_double()
{
	double pi1 = M_PI;
	double pi2 = M_PI;
	cudaSetDoubleForDevice(&pi1);
	return pi1 == pi2;
}

void phelm_init()
{
	if (cublasInit() != CUBLAS_STATUS_SUCCESS)
	{
		fprintf(stderr, "failed to initialize library\n");
		exit(1);
	} else {
		fprintf(stderr, "library initialized\n");
	}
}

void phelm_shutdown()
{
	cublasShutdown();
}

} /* namespace phelm */
