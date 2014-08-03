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

#include "schafe.h"

#if 0
static double f (double u, double x, double y, double t,
                 double mu, double sigma)
{
	double a = exp (t) * sin (y) * sin (2.0 * x);
	double b = -6.0 * exp (t) * sin (y) * sin (2.0 * x);

	return a - mu * b + sigma * a;

//	return (1.0 + 6.0 * mu + sigma) * u;
//	return -u * u * u;
}
#endif


namespace SphereChafe_Private
{

double f (/*double u,*/ double x, double y, double t,
                        double mu, double sigma)
{
	return ipow (cos (x), 2) *
	       (x*cos (y + t) *ipow (cos (x), 2)
	        + sigma*sin (y + t) *ipow (cos (x), 2) *x + 9*mu*sin (y + t) *sin (x) *cos (x)
	        - 15*mu*sin (y + t) *x + 20*mu*sin (y + t) *x*ipow (cos (x), 2) );
}

double
schafe_integrate_cb ( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk, int z,
                      const Mesh & m,
                      int point_i, int point_j,
                      int, int,
                      SphereChafeConfig * d)
{
	double tau   = d->tau_;
	double mu    = d->mu_;
	double sigma = d->sigma_;

	double pt1, pt2;

	pt1  = integrate_cos (phi_i * phi_j, trk, z);
	pt1 *= 1.0 / tau + sigma * 0.5;

	pt2  =  slaplace (phi_j, phi_i, trk, z);
	pt2 *= -0.5 * mu;

	return pt1 + pt2;
}

}

