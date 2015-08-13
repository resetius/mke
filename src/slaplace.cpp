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

VERSION ("$Id$");

using namespace std;
using namespace phelm;

/*
 * phi    = x - latitude
 * lambda = y - longitude
 */

/*
 * laplace = laplace1 + laplace2
 */

static double integrate_cos_func (double x, double y, Polynom * poly)
{
	return cos (x) * poly->apply (x, y);
}

static double integrate_1_cos_func (double x, double y, Polynom * poly)
{
	return poly->apply (x, y) / cos (x);
}

/*
 * laplace1 = (1 / cos phi d / d phi phi_i, cos phi d / phi phi_j)
 */
static double laplace1 (const Polynom & phi_i, const Polynom & phi_j,
                        const Triangle & trk, int z)
{
	Polynom poly1 = diff (phi_i, 0) * diff (phi_j, 0);
	Polynom poly2 = diff (phi_i, 0) * phi_i * phi_j;

	return -integrate_cos(poly1, trk, z);
//	return -integrate_generic (trk, (fxy_t) integrate_cos_func, &poly1)
//		+ integrate_boundary_y(trk, (fxy_t) integrate_cos_func, &poly2)
//		;
}

/*
 * laplace2 = (1 / cos phi d / d lambda phi_i, 1 / cos phi d / d lambda phi_j)
 */
static double laplace2 (const Polynom & phi_i, const Polynom & phi_j,
                        const Triangle & trk, int z)
{
	Polynom poly1 = diff (phi_i, 1) * diff (phi_j, 1);
	Polynom poly2 = diff (phi_i, 1) * phi_j;

	return -integrate_1_cos(poly1, trk, z);
//	return -integrate_generic (trk, (fxy_t) integrate_1_cos_func, &poly1)
//		+ integrate_boundary_x(trk, (fxy_t) integrate_1_cos_func, &poly2); // << ?
}

double slaplace (const Polynom & phi_i, const Polynom & phi_j,
                 const Triangle & trk, int z)
{
	return (laplace1 (phi_i, phi_j, trk, z) + laplace2 (phi_i, phi_j, trk, z) );
}

struct slaplace_right_part_cb_data
{
	const double * F;
	const double * bnd;
};

#if 0
static double
slaplace_right_part_cb ( const Polynom & phi_i,
                         const Polynom & phi_j,
                         const Triangle & trk,
                         const Mesh & m,
                         int point_i,
                         int point_j,
                         int, int,
                         slaplace_right_part_cb_data * d)
{
	const double * F = d->F;
	double b;

	b = F[point_j] * integrate_cos (phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1 && d->bnd)   // на границе
	{
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * slaplace (phi_j, phi_i, trk, m.ps);
	}

	return b;
}
#endif

namespace SphereLaplace_Private
{

double
laplace_bnd1_cb ( const Polynom & phi_i,
                  const Polynom & phi_j,
                  const Triangle & trk,
                  int z,
                  const Mesh & m,
                  int point_i,
                  int point_j,
                  int, int)
{
	Polynom poly = phi_i * phi_j;
	return integrate_cos(poly, trk, z);
	//return integrate_generic (trk, (fxy_t) integrate_cos_func, &poly);
}

double
laplace_bnd2_cb ( const Polynom & phi_i,
                  const Polynom & phi_j,
                  const Triangle & trk,
                  int z,
                  const Mesh & m,
                  int point_i,
                  int point_j,
                  int, int)
{
	return -slaplace (phi_j, phi_i, trk, z);
}

double
slaplace_integrate_cb ( const Polynom & phi_i,
                        const Polynom & phi_j,
                        const Triangle & trk,
				    int z,
                        const Mesh & m,
                        int point_i,
                        int point_j,
                        int, int)
{
	double a = slaplace (phi_j, phi_i, trk, z);
	return a;
}

double id_cb (const Polynom & phi_i,
              const Polynom & phi_j,
              const Triangle & trk,
              int z,
              const Mesh & m,
              int point_i,
              int point_j,
              int, int)
{
	Polynom poly = phi_i * phi_j;
	return integrate_cos(poly, trk, z);
	//return integrate_generic (trk, z, (fxy_t) integrate_cos_func, &poly);
}

}

static double lp_rp (const Polynom & phi_i,
                     const Polynom & phi_j,
                     const Triangle & trk,
                     int z,
                     const Mesh & m,
                     int point_i,
                     int point_j,
                     int, int,
                     slaplace_right_part_cb_data * d)
{
	const double * F = d->F;
	double b = F[point_j] * slaplace (phi_j, phi_i, trk, z);;
#if 0
	if (m.ps_flags[point_j] == 1 && d->bnd)   // на границе
	{
		int j0       = m.p2io[point_j]; //номер внешней точки
		b += - d->bnd[j0] * id_cb (phi_i, phi_j,
		                           trk, m, point_i, point_j, 0);
	}
#endif
	return b;
}

