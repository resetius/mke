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

static double polynom_func (double x, double y, Polynom * p)
{
	return p->apply (x, y);
}

struct laplace_right_part_cb_data
{
	const double * F;
	const double * bnd;
};

double
laplace (const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, const Mesh::points_t & ps)
{
	Polynom poly = diff (phi_j, 0) * diff (phi_i, 0)
	               + diff (phi_j, 1) * diff (phi_i, 1);
	double d1 = -integrate(poly, trk, ps);
	return d1;
}

namespace Laplace_Private
{

double
laplace_bnd2_cb ( const Polynom & phi_i,
                  const Polynom & phi_j,
                  const Triangle & trk,
                  const Mesh & m,
                  int point_i,
                  int point_j,
                  int, int,
                  void * )
{
	return -laplace (phi_i, phi_j, trk, m.ps);
}

double id_cb (const Polynom & phi_i,
              const Polynom & phi_j,
              const Triangle & trk,
              const Mesh & m,
              int point_i, int point_j,
              int, int,
              void *)
{
	return integrate (phi_i * phi_j, trk, m.ps);
}

double
laplace_integrate_cb ( const Polynom & phi_i,
                       const Polynom & phi_j,
                       const Triangle & trk, /* номер треугольника */
                       const Mesh & m,
                       int point_i,
                       int point_j,
                       int, int,
                       void * user_data)
{
	return laplace (phi_i, phi_j, trk, m.ps);
}

} /* namespace */

