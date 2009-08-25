/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (РђР»РµРєСЃРµР№ РћР·РµСЂРёС†РєРёР№)
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gmres.h"
#include "ver.h"

VERSION("$Id$");

void vec_diff(double * r, const double * a, const double * b, int n);
void vec_mult_scalar(double * a, const double * b, double k, int n);
void vec_sum2(double * r, const double * a, const double *b, double k2, int n);

double vec_norm2(const double *v, int n);
double vec_scalar2(const double * a, const double * b, int n);

double get_full_time();

/**
 * Demmel Algorithm  6.9 p 303
 */

static double
algorithm6_9(double * x, const void * A, const double * b, 
			 Ax_t Ax, double eps, int n, int k)
{
	double * q  = 0;
	double * r  = malloc(n * sizeof(double)); /* b - Ax */ // -> device
	double * ax = malloc(n * sizeof(double)); // -> device
	double * h  = 0;
	double * gamma = 0;

	double * s = 0;
	double * c = 0;

	double gamma_0;
	double beta;
	double ret;

	int i, j;

	int hz = k + 1;

	/* x0 = x */
	/* r0 = b - Ax0 */
	Ax(ax, A, x, n);
	vec_diff(r, b, ax, n);

	gamma_0 = vec_norm2(r, n);

	if (gamma_0 <= eps) {
		ret = gamma_0;
		goto end;
	}

	h = calloc(hz * hz, sizeof(double));     // -> host
	q = malloc(hz * n * sizeof(double));     // -> device
	vec_mult_scalar(q, r, 1.0 / gamma_0, n); 
	gamma = malloc(hz * sizeof(double));     // -> host

	s = malloc(hz * sizeof(double));         // -> host
	c = malloc(hz * sizeof(double));         // -> host

	gamma[0] = gamma_0;

	for (j = 0; j < k; ++j) {
		double nr1, nr2;

		Ax(ax, A, &q[j * n], n);
		nr1 = vec_norm2(ax, n);

		for (i = 0; i <= j; ++i) {
			h[i * hz + j] = vec_scalar2(&q[i * n], ax, n); //-> j x j
			//ax = ax - h[i * hz + j] * q[i * n];
			vec_sum2(ax, ax, &q[i * n], -h[i * hz + j], n);
		}
		
		// h -> (j + 1) x j
		h[(j + 1) * hz + j] = vec_norm2(ax, n);

		// loss of orthogonality detected
		// C. T. Kelley: Iterative Methods for Linear and Nonlinear Equations, SIAM, ISBN 0-89871-352-8
		nr2 = 0.001 * h[(j + 1) * hz + j] + nr1;
		if (fabs(nr2  - nr1) < eps) {
			/*fprintf(stderr, "reortho!\n");*/
			for (i = 0; i <= j; ++i) {
				double hr = vec_scalar2(&q[i * n], ax, n);
				h[i * hz + j] += hr;
				vec_sum2(ax, ax, &q[i * n], -hr, n);
			}
			h[(j + 1) * hz + j] = vec_norm2(ax, n);
		}

		// rotate
		for (i = 0; i <= j - 1; ++i) {
			double x = h[i * hz + j];
			double y = h[(i + 1) * hz + j];

			h[i * hz + j]       = x * c[i + 1] + y * s[i + 1];
			h[(i + 1) * hz + j] = x * s[i + 1] - y * c[i + 1];
		}

		beta = sqrt(h[j * hz + j] * h[j * hz + j] + h[(j + 1) * hz + j] * h[(j + 1) * hz + j]);
		s[j + 1]      = h[(j + 1) * hz + j] / beta;
		c[j + 1]      = h[j * hz + j] / beta;
		h[j * hz + j] = beta;

		gamma[j + 1] = s[j + 1] * gamma[j];
		gamma[j]     = c[j + 1] * gamma[j];

		if (gamma[j + 1]  > eps) {
			vec_mult_scalar(&q[(j + 1) * n], ax, 1.0 / h[(j + 1) * hz + j], n);
		} else {
			goto done;
		}
	}

	--j;

done:
	ret = gamma[j + 1];

	{
		double * y = malloc(hz * sizeof(double));  // -> host
		for (i = j; i >= 0; --i) {
			double sum = 0.0;
			for (k = i + 1; k <= j; ++k) {
				sum += h[i * hz + k] * y[k];
			}
			y[i] = (gamma[i] - sum) / h[i * hz + i];
		}

		for (i = 0; i <= j; ++i) {
			vec_sum2(x, x, &q[i * n], y[i], n);
		}
		free(y);
	}

end:
	free(q);  free(r);
	free(ax); free(h);
	free(gamma);
	free(s); free(c);

	return ret;
}

void gmres(double * x, const void * A, const double * b, 
			 Ax_t Ax, int n, int k_dim, int max_it)
{
	double tol = 1e-10;
	double bn  = vec_norm2(b, n);
	int i;

	/* x0 = b */
	memcpy(x, b, n * sizeof(double));
	fprintf(stderr, "  gmres: ||b|| = %le\n", bn);

	if (bn < tol) {
		return;
	}

	for (i = 0; i < max_it; ++i)
	{
		double t1 = get_full_time();
		double e  = algorithm6_9(x, A, b, Ax, tol * bn, n, k_dim);
		double t2 = get_full_time();
		//double xn = vec_norm2(x, n);
		e /= bn;
		fprintf(stderr, "  gmres: iters = %d, eps = %le, t = %lf\n", i, e, 
			(t2 - t1) / 100.0);
		if (e < tol) {
			return;
		}
	}
}


