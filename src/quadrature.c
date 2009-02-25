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

typedef double (*fx_t)(double x, void * data);

/* w_i */
static double g7w[] = 
{
	0.129484966168870,
	0.129484966168870,

	0.279705391489277,
	0.279705391489277,

	0.381830050505119,
	0.381830050505119,

	0.417959183673469,
};

/*x_i, w_i*/
static double k15n[][2] = 
{
	/* gauss: 0, 1 */
	{+0.949107912342759, 0.063092092629979},
	{-0.949107912342759, 0.063092092629979},

	/* gauss: 2, 3 */
	{+0.741531185599394, 0.140653259715525},
	{-0.741531185599394, 0.140653259715525},

	/* gauss: 4, 5 */
	{+0.405845151377397, 0.190350578064785},
	{-0.405845151377397, 0.190350578064785},

	/* gauss: 6 */
	{ 0.000000000000000, 0.209482141084728},

	/* 7, 8 */
	{+0.991455371120813, 0.022935322010529},
	{-0.991455371120813, 0.022935322010529},

	/* 9, 10 */
	{+0.864864423359769, 0.104790010322250},
	{-0.864864423359769, 0.104790010322250},

	/* 11, 12 */
	{+0.586087235467691, 0.169004726639267},
	{-0.586087235467691, 0.169004726639267},

	/* 13, 14 */
	{+0.207784955007898, 0.204432940075298},
	{-0.207784955007898, 0.204432940075298},
};

//#define PRINT_EPS

double gauss_kronrod15(double a, double b, fx_t f, void * d)
{
#ifdef PRINT_EPS
	double g7  = 0.0;
#endif
	double k15 = 0.0;
	double p1  = 0.5 * (b - a);
	double p2  = 0.5 * (b + a);
	int i;

	for (i = 0; i < 7; ++i) {
		double f1 = f(p1 * k15n[i][0] + p2, d);
		k15 += k15n[i][1] * f1;
#ifdef PRINT_EPS
		g7 += g7w[i] * f1;
#endif
	}

	for (i = 7; i < 15; ++i) {
		k15 += k15n[i][1] * f(p1 * k15n[i][0] + p2, d);
	}

	k15 *= p1;
#ifdef PRINT_EPS
	g7 *= p1;
	fprintf(stderr, "eps=%le\n", pow(200.0 * fabs(k15 - g7), 1.5));
#endif
	return k15;
}

