#ifndef UTIL_H
#define UTIL_H
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

/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @section DESCRIPTION
 * Misc functions.
 */


#ifdef __cplusplus
extern "C"
{
#endif

	/**
	 * @defgroup misc Miscellaneous functions and classes.
	 * @{
	 */

	double atan3(double y, double x);

	/**
	 * Take the integral of \f$x^k y^n\f$ over trapezoid.
	 * y=k1x+b1 y=k2x+b2 - bounding lines.
	 * x belongs to segment [x1, x3].
	 \f[
	 \int_{x_1}^{x_3}\int_{k_1 x+b_1}^{k_2 x+b_2} x^k y ^n dx dy
	 \f]
	 * @param k - power
	 * @param n - power
	 * @param k1 - coefficient
	 * @param b1 - coefficient
	 * @param k2 - coefficient
	 * @param b2 - coefficient
	 * @param x1 - the begging of segment
	 * @param x3 - the end of segment
	 */
	double trapezoid_integral (int k, int n,
	                           double k1, double b1,
	                           double k2, double b2,
	                           double x1, double x3);

	/**
	 * Take the integral of \f$x^k y^n cos(x)\f$ over trapezoid.
	 \f[
	 \int_{x_1}^{x_3}\int_{k_1 x+b_1}^{k_2 x+b_2} x^k y ^n cos(x) dx dy
	 \f]
	 * @param k - power
	 * @param n - power
	 * @param k1 - coefficient
	 * @param b1 - coefficient
	 * @param k2 - coefficient
	 * @param b2 - coefficient
	 * @param x1 - the begging of segment
	 * @param x3 - the end of segment
	 */
	double trapezoid_integral_cos (int k, int n,
	                               double k1, double b1,
	                               double k2, double b2,
	                               double x1, double x3);

	/**
	 * Take the integral of \f$x^k y^n sin(x)\f$ over trapezoid.
	 \f[
	 \int_{x_1}^{x_3}\int_{k_1 x+b_1}^{k_2 x+b_2} x^k y ^n sin(x) dx dy
	 \f]
	 * @param k - power
	 * @param n - power
	 * @param k1 - coefficient
	 * @param b1 - coefficient
	 * @param k2 - coefficient
	 * @param b2 - coefficient
	 * @param x1 - the begging of segment
	 * @param x3 - the end of segment
	 */
	double trapezoid_integral_sin (int k, int n,
	                               double k1, double b1,
	                               double k2, double b2,
	                               double x1, double x3);

	/**
	 * Take the integral of \f$x^k y^n / cos(x)\f$ over trapezoid.
	 \f[
	 \int_{x_1}^{x_3}\int_{k_1 x+b_1}^{k_2 x+b_2} \frac{x^k y ^n}{cos(x)} dx dy
	 \f]
	 * @param k - power
	 * @param n - power
	 * @param k1 - coefficient
	 * @param b1 - coefficient
	 * @param k2 - coefficient
	 * @param b2 - coefficient
	 * @param x1 - the begging of segment
	 * @param x3 - the end of segment
	 */
	double trapezoid_integral_1_cos (int k, int n,
	                                 double k1, double b1,
	                                 double k2, double b2,
	                                 double x1, double x3);

	/**
	 * Function-callback that is passed to gauss_kronrod15.
	 * @param x - function argument
	 * @param data - user data
	 * @return the value of f(x)
	 */
	typedef double (*fx_t) (double x, void * data);

	/**
	 * Gauss Kronrod quadrature formula.
	 * http://en.wikipedia.org/wiki/Gauss-Kronrod_quadrature
	 *
	 * @param a
	 * @param b
	 * @param fm - function
	 * @param data - the user data that are passed to the function
	 * @return the integral of fm over the segment [a, b]
	 */
	double gauss_kronrod15 (double a, double b, fx_t fm, void * data);

	/**
	 * writes header.
	 */
	void write_header (int argc, char ** argv, const char * mes);

	/**
	 * burn loop.
	 * @param secs - the number of seconds
	 */
	void burn (double secs);

	/** @} */

#ifdef __cplusplus
}
#endif

#endif /* UTIL_H */

