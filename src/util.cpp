/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (јлексей ќзерицкий)
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
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
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

extern "C" {

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

#include "cnk_6.h"

//це из ен по ка
inline double CNK(int n,int k)
{
	assert(n < 6 && k < 6);
	return Cnk[n][k];
}

inline double 
ssum(double x, int k, int n, double k1, double b1, 
			  double k2, double b2)
{
	double sum1 = 0, sum2 = 0, sum;
	int i;
	if (k2 && b2) {
		for(i = 0; i <= n + 1; i++)
		{
			sum1 += CNK(n + 1, i) * ipow(k2, n + 1 - i)
				* ipow(x, n + k + 2 - i) * ipow(b2, i)
				/ (double)(n + k + 2 - i);
		}
	} else if (k2 && !b2) {
		sum1 = ipow(k2, n + 1) * ipow(x, k + n + 2) / (double)(k + n + 2);
	} else if (!k2 && b2) {
		// ?
		sum1 = ipow(b2, n + 1) * ipow(x, k + 1) / (double)(k + 1);
	} else { //if (!k2 && !b2) {
		sum1 = 0;
	}

	if (k1 && b1) {
		for (i = 0; i <= n + 1; i++)
		{
			sum2 += CNK(n + 1, i) * ipow(k1, n + 1 - i)
				* ipow(x, n + k + 2 - i)
				* ipow(b1, i) / (double)(n + k + 2 - i);
		}
	} else if (k1 && !b1) {
		sum2 = ipow(k1, n + 1) * ipow(x, k + n + 2) / (double)(k + n + 2);
	} else if (!k1 && b1) {
		sum2 = ipow(b1, n + 1) * ipow(x, k + 1) / (double)(k + 1);
	} else { //if (!k1 && !b1) {
		sum2 = 0;
	}

	sum = sum1 - sum2;
	return sum;
}

// интеграл от x^k * y^n по трапеции
// y = k1 x + b1 y = k2 x + b2 
// пр€мые ограничивающие трапецию
// x мен€етс€ от x1 до x3 
double trapezoid_integral(int k, int n, 
						  double k1, double b1, 
						  double k2, double b2, 
						  double x1, double x3)
{
	// надо вычислить интеграл от x1 до x3 от :
	// 1 / (n + 1) x^{k+1} ((k2 x + b2)^{n+1}-(k1 x + b1)^{n+1})
	double retval=(ssum(x3, k, n, k1, b1, k2, b2) 
		- ssum(x1, k, n, k1, b1, k2, b2)) / (double)(n + 1);
	return retval;
}

struct f_cos_data
{
	int k;
	int n;
	double k1;
	double b1;
	double k2;
	double b2;

	f_cos_data(int k_, int n_, double k1_, double b1_,
			double k2_, double b2_):
		k(k_), n(n_), k1(k1_), b1(b1_), k2(k2_), b2(b2_)
	{
	}
};

/* x^k y^{n+1}/{n+1} cos x | */
static double f_cos(double x, f_cos_data * d)
{
	double pt1 = (ipow((d->k2 * x + d->b2), d->n + 1) 
			- ipow((d->k1 * x + d->b1), d->n + 1)) / (double)(d->n + 1);
	return ipow(x, d->k) * pt1 * cos(x);
}

/* x^k y^{n+1}/{n+1} sin x | */
static double f_sin(double x, f_cos_data * d)
{
	double pt1 = (ipow((d->k2 * x + d->b2), d->n + 1)
			- ipow((d->k1 * x + d->b1), d->n + 1)) / (double)(d->n + 1);
	return ipow(x, d->k) * pt1 * sin(x);
}

/* x^k y^{n+1}/{n+1}/cos x | */
static double f_1_cos(double x, f_cos_data * d)
{
	double pt1 = (ipow((d->k2 * x + d->b2), d->n + 1)
			- ipow((d->k1 * x + d->b1), d->n + 1)) / (double)(d->n + 1);
	return ipow(x, d->k) * pt1 / cos(x);
}

/**
 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n cos x) dx dy
 */
double trapezoid_integral_cos(int k, int n,
		double k1, double b1,
		double k2, double b2,
		double x1, double x3)
{
	f_cos_data data(k, n, k1, b1, k2, b2);
	return gauss_kronrod15(x1, x3, (fx_t)f_cos, &data);
}

/**
 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n sin x) dx dy
 */
double trapezoid_integral_sin(int k, int n,
	double k1, double b1,
	double k2, double b2,
	double x1, double x3)
{
	f_cos_data data(k, n, k1, b1, k2, b2);
	return gauss_kronrod15(x1, x3, (fx_t)f_sin, &data);
}

/**
 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n/cos x) dx dy
 */
double trapezoid_integral_1_cos(int k, int n,
		double k1, double b1,
		double k2, double b2,
		double x1, double x3)
{
	f_cos_data data(k, n, k1, b1, k2, b2);
	return gauss_kronrod15(x1, x3, (fx_t)f_1_cos, &data);
}

void mat_print(const double * A, int n)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%9.2le ", A[i * n + j]);
		}
		printf("\n");
	}
}

void vec_print(const double * A, int n)
{
	for (int i = 0; i < n; ++i) {
		printf("%9.2le ", A[i]);
	}
	printf("\n");
}

void mat_mult_vector(double * r, const double * A, const double * x, int n)
{
	int i, j;
	for (i = 0; i < n; ++i) {
		r[i] = 0.0;
		for (j = 0; j < n; ++j) {
			r[i] += A[i * n + j] * x[j];
		}
	}
}

void sparse_mult_vector_l(double * r, const struct Sparse * A, const double * x, int n)
{
	int j;

#pragma omp parallel for
	for (j = 0; j < n; ++j) {
		double *p = &A->Ax[A->Ap[j]];
		int i0;
		r[j] = 0;

		for (i0 = A->Ap[j]; i0 < A->Ap[j + 1]; ++i0, ++p) {
			int i = A->Ai[i0];
			r[j] += *p * x[i];
		}
	}
}

void sparse_mult_vector_r(double * r, const struct Sparse * A, const double * x, int n)
{
	int i0, i, j;
	const double * p = A->Ax;

	memset(r, 0, n * sizeof(double));
	for (j = 0; j < n; ++j) {
		for (i0 = A->Ap[j]; i0 < A->Ap[j + 1]; ++i0, ++p) {
			i = A->Ai[i0];
			r[i] += *p * x[j];
		}
	}
}

/**
 * r = k1 * a + k2 * b
 */
void vec_sum1(double * r, const double * a, const double *b, double k1, double k2, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i) {
		r[i] = k1 * a[i] + k2 * b[i];
	}
}

void vec_sum(double * r, const double * a, const double *b, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i) {
		r[i] = a[i] + b[i];
	}
}

void vec_mult(double * r, const double * a, const double *b, int n)
{
	int i;
//#pragma omp parallel for
	for (i = 0; i < n; ++i) {
		r[i] = a[i] * b[i];
	}
}

/**
 * a = b * k
 */
void vec_mult_scalar(double * a, const double * b, double k, int n)
{
	int i;
//#pragma omp parallel for
	for (i = 0; i < n; ++i) {
		a[i] = b[i] * k;
	}
}

/**
 * r = a - b
 */
void vec_diff(double * r, const double * a, const double * b, int n)
{
	int i;
//#pragma omp parallel for
	for (i = 0; i < n; ++i) {
		r[i] = a[i] - b[i];
	}
}

void sparse_print(const struct Sparse * A, int n, FILE * f)
{
	int i, i0, j, k, i_old;
	const double * p = A->Ax;
	for (j = 0; j < n; ++j) {
		i_old = -1;
		for (i0 = A->Ap[j]; i0 < A->Ap[j + 1]; ++i0, ++p) {
			i = A->Ai[i0];
			for (k = i_old; k < i - 1; ++k) {
				fprintf(f, "%8.3lf ", 0.0);
			}
			fprintf(f, "%8.3lf ", *p);
			i_old = i;
		}

		for (k = i_old + 1; k < n; ++k) {
			fprintf(f, "%8.3lf ", 0.0);
		}
		fprintf(f, "\n");
	}
}


#ifdef WIN32
#include <windows.h>
void set_fpe_except()
{
	int cw = _controlfp(0, 0);
	cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);
	_controlfp(cw, MCW_EM);
}
#else


/*
	For reference, the layout of the MXCSR register:
	FZ:RC:RC:PM:UM:OM:ZM:DM:IM:Rsvd:PE:UE:OE:ZE:DE:IE
	15 14 13 12 11 10  9  8  7   6   5  4  3  2  1  0

	And the layout of the 387 FPU control word register:
	Rsvd:Rsvd:Rsvd:X:RC:RC:PC:PC:Rsvd:Rsvd:PM:UM:OM:ZM:DM:IM
	 15   14   13 12 11 10  9  8   7    6   5  4  3  2  1  0

	Where:
		Rsvd - Reserved
		FZ   - Flush to Zero
		RC   - Rounding Control
		PM   - Precision Mask
		UM   - Underflow Mask
		OM   - Overflow Mask
		ZM   - Zerodivide Mask
		DM   - Denormal Mask
		IM   - Invalid Mask
		PE   - Precision Exception
		UE   - Underflow Exception
		OE   - Overflow Exception
		ZE   - Zerodivide Exception
		DE   - Denormal Exception
		IE   - Invalid Exception
		X    - Infinity control (unused on 387 and higher)
		PC   - Precision Control

	Source: Intel Architecture Software Development Manual, Volume 1, Basic Architecture
*/

void set_fpe_except()
{
	//http://www.website.masmforum.com/tutorials/fptute/fpuchap1.htm
	unsigned m;
	asm ("fstcw %0" : : "m" (*&m));
	asm ("fwait");
	m ^= 0x1f; //turn on all exceptions!
	m |= (1 << 8);
	m |= (1 << 9);
	m ^= (1 << 8);
	m ^= (1 << 9);
	m |= (1 << 9); // double precision
	asm ("fldcw %0" : : "m" (*&m));
	asm ("fwait");

	m = 1 << 12;
	asm ("ldmxcsr %0" : : "m" (*&m));
}
#endif

#include <time.h>
 
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

#ifdef _MSC_VER
#include <windows.h>

struct timezone 
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
	FILETIME ft;
	unsigned __int64 tmpres = 0;
	static int tzflag;

	if (NULL != tv)
	{
		GetSystemTimeAsFileTime(&ft);

		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;

		/*converting file time to unix epoch*/
		tmpres /= 10;  /*convert into microseconds*/
		tmpres -= DELTA_EPOCH_IN_MICROSECS; 
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	}

	if (NULL != tz)
	{
		if (!tzflag)
		{
			_tzset();
			tzflag++;
		}
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	}

	return 0;
}
#else
#include <sys/time.h>
#endif

double get_full_time()
{
	struct timeval tv;
	gettimeofday(&tv, 0);
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}

} /* extern "C" */

