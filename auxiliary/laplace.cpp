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

#include <assert.h>
#include <math.h>

#include "phelm.h"
#include "util.h"
#include "solver.h"
#include "laplace.h"
#include "ver.h"

VERSION("$Id$");

using namespace std;
using namespace phelm;

struct laplace_right_part_cb_data
{
	const double * F;
	const double * bnd;
};

double
laplace(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, const Mesh::points_t & ps)
{
	Polynom poly = diff(phi_j, 0) * diff(phi_i, 0)
		+ diff(phi_j, 1) * diff(phi_i, 1);
	return -integrate(poly, trk, ps);
}

static double 
laplace_right_part_cb( const Polynom & phi_i,
                       const Polynom & phi_j,
                       const Triangle & trk, /* номер треугольника */
                       const Mesh & m,
                       int point_i,
                       int point_j,
                       int, int,
                       laplace_right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

	//b = F[point_j] * integrate(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) {         // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += - bnd[j0] * laplace(phi_i, phi_j, trk, m.ps);
	}
	else {
		b += F[point_j] * integrate(phi_i * phi_j, trk, m.ps);
	}

	return b;
}

namespace Laplace_Private {

double
laplace_bnd1_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int, int,
		void * d)
{
	return integrate(phi_i * phi_j, trk, m.ps);
}

double
laplace_bnd2_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int, int,
		void * )
{
	return -laplace(phi_i, phi_j, trk, m.ps);
}

double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i, int point_j,
		int, int,
		void *)
{
	return integrate(phi_i * phi_j, trk, m.ps);
}

double 
laplace_integrate_cb( const Polynom & phi_i,
                      const Polynom & phi_j, 
                      const Triangle & trk, /* номер треугольника */
                      const Mesh & m,
                      int point_i,
                      int point_j,
                      int, int,
                      void * user_data)
{
	return laplace(phi_i, phi_j, trk, m.ps);
}

} /* namespace */

namespace Chafe_Private {

double 
chafe_integrate_cb( const Polynom & phi_i,
                    const Polynom & phi_j, 
                    const Triangle & trk, 
                    const Mesh & m, int point_i, int point_j,
                    int, int,
                    const ChafeConfig * d)
{
	double tau   = d->tau_;
	double mu    = d->mu_;
	double sigma = d->sigma_;

	double pt1, pt2;

	pt1  = integrate(phi_j * phi_i, trk, m.ps);
	pt1 *= 1.0 / tau + sigma * 0.5;

	pt2  = laplace(phi_i, phi_j, trk, m.ps);
	pt2 *= -0.5 * mu;

	return pt1 + pt2;
}

} /* namespace Chafe_Private */

static double lp_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i, int point_j,
		laplace_right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

	b = F[point_j] * laplace(phi_i, phi_j, trk, m.ps);

	//return F[point_j] * laplace(phi_i, phi_j, trk, m);
#if 0
	if (m.ps_flags[point_j] == 1)
	{
		b = F[point_j] * id_cb(phi_i, phi_j, trk, m, point_i, point_j, 0);
	}

	//if (m.ps_flags[point_j] == 1 && d->bnd)
	//{
	//	int j0       = m.p2io[point_j];
	//	b += - d->bnd[j0] * id_cb(phi_i, phi_j, 
	//			trk, m, point_i, point_j, 0);
	//}
#endif
	return b;
}
