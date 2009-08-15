#ifndef GMRES_H
#define GMRES_H
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
 
/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @section DESCRIPTION
 *
 * Generalized minimal residual method (GMRES).
 * @see { http://en.wikipedia.org/wiki/GMRES }
 */

#include <math.h>

#include "base.h"
#include "linal.h"
#include "util.h"

namespace phelm {

/**
 * Demmel Algorithm  6.9 p 303
 */

template < typename T, typename Mat, typename Ax_t >
T algorithm6_9(T * x, const Mat * A, const T * b, 
			 Ax_t Ax, T eps, int n, int k)
{
	ArrayDevice < T > q;        // -> device
	ArrayDevice < T > r(n);     /* b - Ax */ // -> device
	ArrayDevice < T > ax(n);    // -> device

	ArrayHost < T > h;
	ArrayHost < T > gamma;

	ArrayHost < T > s;
	ArrayHost < T > c;

	T gamma_0;
	T beta;
	T ret = 1;

	int i, j;

	int hz = k + 1;

	/* x0 = x */
	/* r0 = b - Ax0 */
	Ax(&ax[0], A, x, n);
	vec_diff(&r[0], &b[0], &ax[0], n);

	gamma_0 = vec_norm2(&r[0], n);

	if (gamma_0 <= eps) {
		ret = gamma_0;
		goto end;
	}

	h.resize(hz * hz);    // -> host
	q.resize(hz * n);     // -> device
	vec_mult_scalar(&q[0], &r[0], (T)1.0 / gamma_0, n); 
	gamma.resize(hz);     // -> host

	s.resize(hz);         // -> host
	c.resize(hz);         // -> host

	gamma[0] = gamma_0;

	for (j = 0; j < k; ++j) {
		T nr1, nr2;

		Ax(&ax[0], A, &q[j * n], n);
		nr1 = vec_norm2(&ax[0], n);

		for (i = 0; i <= j; ++i) {
			h[i * hz + j] = vec_scalar2(&q[i * n], &ax[0], n); //-> j x j
			//ax = ax - h[i * hz + j] * q[i * n];
			vec_sum2(&ax[0], &ax[0], &q[i * n], -h[i * hz + j], n);
		}
		
		// h -> (j + 1) x j
		h[(j + 1) * hz + j] = vec_norm2(&ax[0], n);

		// loss of orthogonality detected
		// C. T. Kelley: Iterative Methods for Linear and Nonlinear Equations, SIAM, ISBN 0-89871-352-8
		nr2 = (T)0.001 * h[(j + 1) * hz + j] + nr1;
		if (fabs(nr2  - nr1) < eps) {
			/*fprintf(stderr, "reortho!\n");*/
			for (i = 0; i <= j; ++i) {
				T hr = vec_scalar2(&q[i * n], &ax[0], n);
				h[i * hz + j] += hr;
				vec_sum2(&ax[0], &ax[0], &q[i * n], -hr, n);
			}
			h[(j + 1) * hz + j] = vec_norm2(&ax[0], n);
		}

		// rotate
		for (i = 0; i <= j - 1; ++i) {
			T x = h[i * hz + j];
			T y = h[(i + 1) * hz + j];

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
			vec_mult_scalar(&q[(j + 1) * n], &ax[0], (T)1.0 / h[(j + 1) * hz + j], n);
		} else {
			goto done;
		}
	}

	--j;

done:
	ret = gamma[j + 1];

	{
		ArrayHost < T > y(hz);  // -> host
		for (i = j; i >= 0; --i) {
			T sum = 0.0;
			for (k = i + 1; k <= j; ++k) {
				sum += h[i * hz + k] * y[k];
			}
			y[i] = (gamma[i] - sum) / h[i * hz + i];
		}

		for (i = 0; i <= j; ++i) {
			vec_sum2(&x[0], &x[0], &q[i * n], y[i], n);
		}
	}

end:
	return ret;
}

/**
 * Solve equation Ax=b.
 *
 * @param x - answer
 * @param A - matrix
 * @param b - right part
 * @param Ax - callback that calculates y=Ax
 * @param n  - dimension of system
 * @param k_dim  - Krylov dimension
 * @param max_it - maximum number of iterations
 */
template < typename T, typename Mat, typename Ax_t >
void gmres(T * x, const Mat * A, const T * b, 
			 Ax_t Ax, int n, int k_dim, int max_it)
{
	T tol = (sizeof(T) >= 8) ? 1e-12 : 1e-5f;
	T bn  = vec_norm2(b, n);
	int i;

	/* x0 = b */
	vec_copy(x, b, n);
	fprintf(stderr, "  gmres: ||b|| = %le\n", (double)bn);

	if (bn < tol) {
		return;
	}

	for (i = 0; i < max_it; ++i)
	{
		Timer t;
		T e  = algorithm6_9(x, A, b, Ax, tol * bn, n, k_dim);
		e /= bn;
		fprintf(stderr, "  gmres: iters = %d, eps = %le, t = %lf\n", i, e, 
			t.elapsed() / 100.0);
		if (e < tol) {
			return;
		}
	}
}

} /* namespace */

#endif /* GMRES_H */
