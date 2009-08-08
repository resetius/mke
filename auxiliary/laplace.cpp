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

}

static double f(double u, double x, double y, double t, double mu, double sigma)
{
	// for test
	// f(u) = \du/\dt -mu \Delta u + \sigma u
	double lapl = exp(t) * (2.0 * y * y - 2.0 * y + 2.0 * x * x - 2.0 * x);
	double uu = exp(t) * x * (x - 1) * y * (y - 1);
	return uu - mu * lapl + sigma * uu;
//	return u - mu * lapl + sigma * u;
//	return -u * u * u;
}

static double 
chafe_integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk, 
                     const Mesh & m, int point_i, int point_j,
					 int, int,
		     const Chafe * d)
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

struct chafe_right_part_cb_data
{
	const double * F;
	const double * bnd;
	const Chafe * d;
};

double 
chafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
					  int i, int j,
                      chafe_right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

//	b = F[point_j] * integrate(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += - bnd[j0] * chafe_integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, i, j, d->d);
	} 
	else {
		b += F[point_j] * integrate(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

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

Chafe::Chafe(const Mesh & m, double tau, double sigma, double mu)
	: m_(m), laplace_(m), A_((int)m.inner.size()), 
	tau_(tau), mu_(mu), sigma_(sigma)
{
	/* Матрица левой части */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, chafe_integrate_cb, this);
}

/* \frac{du}{dt} = \mu \Delta u - \sigma u + f (u) */
void Chafe::solve(double * Ans, const double * X0,
						const double * bnd, double t)
{
	int rs  = (int)m_.inner.size();
	int sz  = (int)m_.ps.size();
	vec u(rs);
	vec p(sz);
	vec delta_u(rs);
	vec rp(rs);

	// генерируем правую часть
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	u2p(&u[0], X0, m_);
	laplace_.calc2(&delta_u[0], X0);

	// u/dt + mu \Delta u / 2
	vec_sum1(&delta_u[0], &u[0], &delta_u[0], 1.0 / tau_, mu_ * 0.5, rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vec_sum1(&delta_u[0], &delta_u[0], &u[0], 1.0, -sigma_ * 0.5, rs);

	// TODO: тут надо сделать метод простой итерации
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		int point = m_.inner[i];
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();

		u[i] = delta_u[i] + f(u[i], x, y, t, mu_, sigma_);
	}

	// правую часть на границе не знаем !
	p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	chafe_right_part_cb_data data2;
	data2.F   = &p[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		chafe_right_part_cb, &data2);

	phelm::solve(Ans, bnd, &rp[0], A_, m_);
}

