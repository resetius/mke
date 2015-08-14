/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2015 Alexey Ozeritsky
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "util.h"
#include "ver.h"
#include "timer.h"
#include "pow.h"

VERSION ("$Id$");

extern "C"
{
	double atan3(double y, double x) {
		if (fabs(x) < 1e-15) {
			if (y > 0) {
				return M_PI_2;
			}
			else if (y < 0) {
				return M_PI + M_PI_2;
			}
		}

		if (x > 0 && y >= 0) {
			return atan(y / x);
		}
		if (x < 0 && y >= 0) {
			return M_PI + atan(y / x);
		}
		if (x < 0 && y < 0) {
			return M_PI + atan(y / x);
		}
		if (x > 0 && y < 0) {
			return M_PI + M_PI + atan(y / x);
		}
		//abort();
		return 0.0;
	}

	void write_header (int argc, char ** argv, const char * mes)
	{
		fprintf (stderr, "#%s\n", mes);
		fprintf (stderr, "#cmd:");
		for (int i = 0; i < argc; ++i)
		{
			fprintf (stderr, "%s ", argv[i]);
		}
		fprintf (stderr, "\n");
	}

#include "cnk_6.h"

//це из ен по ка
	inline double CNK (int n, int k)
	{
		assert (n < 10 && k < 10);
		return Cnk[n][k];
	}

	inline double
	ssum (double x, int k, int n, double k1, double b1,
	      double k2, double b2)
	{
		double sum1 = 0, sum2 = 0, sum;
		int i;
		if (k2 && b2)
		{
			for (i = 0; i <= n + 1; i++)
			{
				sum1 += CNK (n + 1, i) * ipow (k2, n + 1 - i)
				        * ipow (x, n + k + 2 - i) * ipow (b2, i)
				        / (double) (n + k + 2 - i);
			}
		}
		else if (k2 && !b2)
		{
			sum1 = ipow (k2, n + 1) * ipow (x, k + n + 2) / (double) (k + n + 2);
		}
		else if (!k2 && b2)
		{
			// ?
			sum1 = ipow (b2, n + 1) * ipow (x, k + 1) / (double) (k + 1);
		}
		else   //if (!k2 && !b2) {
		{
			sum1 = 0;
		}

		if (k1 && b1)
		{
			for (i = 0; i <= n + 1; i++)
			{
				sum2 += CNK (n + 1, i) * ipow (k1, n + 1 - i)
				        * ipow (x, n + k + 2 - i)
				        * ipow (b1, i) / (double) (n + k + 2 - i);
			}
		}
		else if (k1 && !b1)
		{
			sum2 = ipow (k1, n + 1) * ipow (x, k + n + 2) / (double) (k + n + 2);
		}
		else if (!k1 && b1)
		{
			sum2 = ipow (b1, n + 1) * ipow (x, k + 1) / (double) (k + 1);
		}
		else   //if (!k1 && !b1) {
		{
			sum2 = 0;
		}

		sum = sum1 - sum2;
		return sum;
	}

// интеграл от x^k * y^n по трапеции
// y = k1 x + b1 y = k2 x + b2
// прямые ограничивающие трапецию
// x меняется от x1 до x3
	double trapezoid_integral (int k, int n,
	                           double k1, double b1,
	                           double k2, double b2,
	                           double x1, double x3)
	{
		// надо вычислить интеграл от x1 до x3 от :
		// 1 / (n + 1) x^{k+1} ((k2 x + b2)^{n+1}-(k1 x + b1)^{n+1})
		double retval = (ssum (x3, k, n, k1, b1, k2, b2)
		                 - ssum (x1, k, n, k1, b1, k2, b2) ) / (double) (n + 1);
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

		f_cos_data (int k_, int n_, double k1_, double b1_,
		            double k2_, double b2_) :
				k (k_), n (n_), k1 (k1_), b1 (b1_), k2 (k2_), b2 (b2_)
		{
		}
	};

	/* x^k y^{n+1}/{n+1} cos x | */
	static double f_cos (double x, f_cos_data * d)
	{
		double pt1 = (ipow ( (d->k2 * x + d->b2), d->n + 1)
		              - ipow ( (d->k1 * x + d->b1), d->n + 1) ) / (double) (d->n + 1);
		return ipow (x, d->k) * pt1 * cos (x);
	}

	/* x^k y^{n+1}/{n+1} sin x | */
	static double f_sin (double x, f_cos_data * d)
	{
		double pt1 = (ipow ( (d->k2 * x + d->b2), d->n + 1)
		              - ipow ( (d->k1 * x + d->b1), d->n + 1) ) / (double) (d->n + 1);
		return ipow (x, d->k) * pt1 * sin (x);
	}

	/* x^k y^{n+1}/{n+1}/cos x | */
	static double f_1_cos (double x, f_cos_data * d)
	{
		double pt1 = (ipow ( (d->k2 * x + d->b2), d->n + 1)
		              - ipow ( (d->k1 * x + d->b1), d->n + 1) ) / (double) (d->n + 1);
		return ipow (x, d->k) * pt1 / cos (x);
	}

	/**
	 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n cos x) dx dy
	 */
	double trapezoid_integral_cos (int k, int n,
	                               double k1, double b1,
	                               double k2, double b2,
	                               double x1, double x3)
	{
		f_cos_data data (k, n, k1, b1, k2, b2);
		return gauss_kronrod15 (x1, x3, (fx_t) f_cos, &data);
	}

	/**
	 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n sin x) dx dy
	 */
	double trapezoid_integral_sin (int k, int n,
	                               double k1, double b1,
	                               double k2, double b2,
	                               double x1, double x3)
	{
		f_cos_data data (k, n, k1, b1, k2, b2);
		return gauss_kronrod15 (x1, x3, (fx_t) f_sin, &data);
	}

	/**
	 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n/cos x) dx dy
	 */
	double trapezoid_integral_1_cos (int k, int n,
	                                 double k1, double b1,
	                                 double k2, double b2,
	                                 double x1, double x3)
	{
		f_cos_data data (k, n, k1, b1, k2, b2);
		return gauss_kronrod15 (x1, x3, (fx_t) f_1_cos, &data);
	}

	void burn (double secs)
	{
		double t1 = get_full_time();
#pragma omp parallel sections
		{
#pragma omp section
			for (;;)
			{
				double t2 = get_full_time();
				if (t2 - t1 > secs * 100)
				{
					break;
				}
			}
#pragma omp section
			for (;;)
			{
				double t2 = get_full_time();
				if (t2 - t1 > secs * 100)
				{
					break;
				}
			}
#pragma omp section
			for (;;)
			{
				double t2 = get_full_time();
				if (t2 - t1 > secs * 100)
				{
					break;
				}
			}
#pragma omp section
			for (;;)
			{
				double t2 = get_full_time();
				if (t2 - t1 > secs * 100)
				{
					break;
				}
			}
		}
	}

} /* extern "C" */

