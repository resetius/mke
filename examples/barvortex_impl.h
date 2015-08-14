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

#include <assert.h>
#include <vector>
#include <math.h>

#include "util.h"
#include "ver.h"

VERSION ("$Id$");

using namespace phelm;
using std::vector;

#ifndef BARVORTEX_H
#error "do not include this file!"
#endif

namespace bv_private
{

template < typename BV >
double
integrate_cb ( const Polynom & phi_i,
               const Polynom & phi_j,
               const Triangle & trk, int z,
               const Mesh & m,
               int point_i, int point_j,
               int i, int j,
               const typename BV::conf_t * d)
{
	double tau   = d->tau;
	double mu    = d->mu;
	double sigma = d->sigma;

	double pt1, pt2;

	pt1  = BV::integrate (phi_i * phi_j, trk, z);
	pt1 *= 1.0 / tau + sigma * d->theta;

	pt2  =  BV::laplace (phi_j, phi_i, trk, z);
	pt2 *= - d->theta * mu;

	return pt1 + pt2;
}

template < typename BV >
double
integrate_backward_cb ( const Polynom & phi_i,
                        const Polynom & phi_j,
                        const Triangle & trk, int z,
                        const Mesh & m,
                        int point_i, int point_j,
                        int i, int j,
                        const typename BV::conf_t * d)
{
	double tau   = d->tau;
	double mu    = d->mu;
	double sigma = d->sigma;

	double pt1, pt2;

	pt1  = BV::integrate (phi_i * phi_j, trk, z);
	pt1 *= 1.0 / tau - sigma * (1.0 - d->theta);

	pt2  =  BV::laplace (phi_j, phi_i, trk, z);
	pt2 *= (1.0 - d->theta) * mu;

	return pt1 + pt2;
}

}

template < typename L, typename J, typename N >
BarVortex < L, J, N > ::BarVortex (const Mesh & m, const BarVortexConf & conf)
		: N (m), m_ (m), conf_(conf), l_ (m), j_ (m),
		A_ (m.inner_size),
		bnd_ (m.inner_size),
		Ab_ (m.inner_size),
		bndb_ (m.inner_size)
{
	int sz = m_.size;
	int os = m_.outer_size;
	int rs = m_.inner_size;

	lh_.resize (sz);
	proj (&lh_[0], m_, conf.coriolis);

	/* Матрица левой части совпадает с Чафе-Инфантом на сфере */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/

	auto func = [&conf](
		const Polynom & phi_i, const Polynom & phi_j,
		const Triangle & tr, int zone, const Mesh & m,
		int point_i, int point_j, int i, int j
		) 
	{
		return bv_private::integrate_cb<BarVortex>(
			phi_i, phi_j, tr, zone, m,
			point_i, point_j, i, j, &conf);
	};

	auto func_backward = [&conf](
		const Polynom & phi_i, const Polynom & phi_j,
		const Triangle & tr, int zone, const Mesh & m,
		int point_i, int point_j, int i, int j
		)
	{
		return bv_private::integrate_backward_cb<BarVortex>(
			phi_i, phi_j, tr, zone, m,
			point_i, point_j, i, j, &conf);
	};

	generate_matrix (A_, m_, func);
	generate_boundary_matrix (bnd_, m_, func);

	generate_matrix (Ab_, m_, func_backward);
	generate_boundary_matrix(bndb_, m_, func_backward);
}

template < typename L, typename J, typename N >
void BarVortex < L, J, N > ::info()
{
	fprintf (stderr, "#tau:%.16lf\n", conf_.tau);
	fprintf (stderr, "#sigma:%.16lf\n", conf_.sigma);
	fprintf (stderr, "#mu:%.16lf\n", conf_.mu);
	fprintf (stderr, "#k1:%.16lf\n", conf_.k1);
	fprintf (stderr, "#k2:%.16lf\n", conf_.k2);
	fprintf (stderr, "#theta:%.16lf\n", conf_.theta);
}

/*
 * d L(phi)/dt + k1 J(psi, L(psi)) + k2 J(psi, l + h) + sigma L(psi) - mu LL(psi) = f(phi, la)
 * L = Laplace
 */
template < typename L, typename J, typename N >
void BarVortex < L, J, N > ::calc (double * u1,
                                const double * u,
                                const double * bnd_u,
                                const double * bnd_w,
                                double t)
{
	int rs = m_.inner_size; // размерность внутренней области
	int os = m_.outer_size; // размерность внешней области
	int sz = m_.size;       // размерность полная

	double nr0 = N::norm (u);

	double k1 = conf_.k1;
	double k2 = conf_.k2;
	double theta = conf_.theta;
	double mu = conf_.mu;
	double sigma = conf_.sigma;
	double tau = conf_.tau;

	vector < double > w (sz);      // w = L(u)
	vector < double > dw (sz);     // dw = L(w) = LL(u)
	vector < double > FC (sz);      // правая часть
	vector < double > FC1 (rs);     // правая часть

	// next
	vector < double > u_n (sz);
	vector < double > w_n (sz);
	vector < double > u_n1 (sz);

	// tmp
	vector < double > tmp1 (sz);
	vector < double > tmp2 (sz);
	vector < double > tmp3 (2*rs);

	// inner (!!!) jacobian
	vector < double > jac (rs);
	vector < double > jac1 (rs);
	vector < double > jac2 (rs);

	//
	vector < double > F (rs);
	vector < double > rp (2*rs);
	vector < double > ans (2*rs);

	vector < double > total_bnd (2*rs);
	if (bnd_w)
	{
		memcpy (&total_bnd[0],  bnd_w, os);
	}

	if (bnd_u)
	{
		memcpy (&total_bnd[os], bnd_u, os);
	}

	// генерируем правую часть
	// w/dt + mu (1-theta) L w - \sigma(1-theta) w -
	// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h) + f(x, y)

	// w = L (u)
	l_.calc1 (&w[0], &u[0], bnd_w);

	// dw = L (w)
	l_.calc1 (&dw[0], &w[0], 0);

	// w/dt + mu (1-theta) L w
	vec_sum1 (&FC[0], &w[0], &dw[0], 1.0 / tau,
	          mu * (1.0 - theta), sz);

	// w/dt + mu (1-theta) L w - \sigma (1-theta) w
	vec_sum1 (&FC[0], &FC[0], &w[0], 1.0,
	          -sigma * (1.0 - theta), sz);

#pragma omp parallel for
	for (int i = 0; i < sz; ++i)
	{
		int point = i;
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();

		if (conf_.rp)
		{
			FC[i] += conf_.rp (x, y, t, mu, sigma);
		}
	}

	// в FC содержится правая часть, которая не меняется при итерациях!

	for (int it = 0; it < 1000; ++it)
	{
		//   k1 J(0.5(u+u), 0.5(w+w)) + k2 J(0.5(u+u), l + h)   =
		// = J(0.5 (u+u), 0.5 k1 (w+w)) + J(0.5 (u+u), k2 (l + h)) =
		// = J(0.5 (u+u), 0.5 k1 (w+w) + k2 (l + h))
		// 0.5 k1 (w+w) + k2 (l + h) <- для вычисления Якобиана это надо знать и на границе!
		//
		vec_sum1 (&tmp1[0], &w_n[0], &w[0], k1 * theta,
		          k1 * (1.0 - theta), sz);
		vec_sum2 (&tmp1[0], &tmp1[0], &lh_[0], k2, sz);
		// 0.5(u+u)
		vec_sum1 (&tmp2[0], &u_n[0], &u[0], theta,
		          1.0 - theta, sz);
		// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h)
		j_.calc2 (&jac[0], &tmp2[0], &tmp1[0]);

		// w/dt + mu (1-theta) L w  - \sigma (1-theta) w -
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)

#pragma omp parallel for
		for (int i = 0; i < rs; ++i)
		{
			int point = m_.inner[i];
			F[i] = FC[point]  - jac[i];
		}

		memset (&rp[0], 0, rs * sizeof (double) );
		l_.idt_.mult_vector (&rp[0], &F[0]);
		if (bnd_w)
		{
			memset (&tmp1[0], 0, rs * sizeof (double) );
			bnd_.mult_vector (&tmp1[0], bnd_w);
			vec_sum (&rp[0], &rp[0], &tmp1[0], (int) rp.size() );
		}

		// тут граничное условие на омега!
		solve (&w_n[0], bnd_w, &rp[0], A_, m_);
		// а тут граничное условие на пси!
		l_.solve (&u_n1[0], &w_n[0], bnd_u);

		double nr = N::dist (&u_n1[0], &u_n[0]);

		u_n1.swap (u_n);
		if (nr / nr0 < 1e-14 || isnan (nr) )
		{
			break;
		}
	}
	memcpy (u1, &u_n[0], sz * sizeof (double) );
	//phelm::smooth1(&u1[0], &u_n[0], m_);
}

double l_rp (double x, double y, double t, double sigma, double mu)
{
	return -9*cos (y + t) *
	       ipow (cos (x), 3) *sin (x) + 15*x*cos (y + t) *
	       ipow (cos (x), 2) - 20*x*cos (y + t) *
	       ipow (cos (x), 4) - 9*sigma*sin (y + t) *
	       ipow (cos (x), 3) *sin (x) + 15*sigma*sin (y + t) *
	       ipow (cos (x), 2) *x - 20*sigma*sin (y + t) *
	       ipow (cos (x), 4) *x - 360*mu*sin (y + t) *sin (x) *
	       ipow (cos (x), 3) + 390*mu*sin (y + t) *x*
	       ipow (cos (x), 2) - 400*mu*sin (y + t) *x*
	       ipow (cos (x), 4) + 147*mu*sin (y + t) *sin (x) *cos (x) - 45*mu*sin (y + t) *x;
}

template < typename L, typename J, typename N >
void BarVortex < L, J, N >::calc_L (double * u1, const double * u, const double * z,
                                 const double * bnd_u,
                                 const double * bnd_w,
                                 double t)
{
	int rs = (int) m_.inner.size(); // размерность внутренней области
	int sz = (int) m_.ps.size();   // размерность полная
	double k1 = conf_.k1;
	double k2 = conf_.k2;
	double theta = conf_.theta;
	double mu = conf_.mu;
	double sigma = conf_.sigma;
	double tau = conf_.tau;

	vector < double > z_lapl (sz);

	vector < double > pt1 (sz); //лаплас, умноженный на коэф
	vector < double > pt2 (sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3 (sz); //якобиан, умноженный на коэф

	l_.calc1 (&z_lapl[0], z, bnd_w);
	l_.calc1 (&pt1[0], u, bnd_w); //первая часть - лаплас, умноженный на коэф,

	{
		vector < double > jac1 (sz);
		vector < double > jac2 (sz);
		vector < double > jac3 (sz);

		j_.calc1 (&jac1[0], u, &lh_[0], 0);
		j_.calc1 (&jac2[0], z, &pt1[0], 0);
		j_.calc1 (&jac3[0], u, &z_lapl[0], 0);

		vec_sum2 (&pt3[0], &pt3[0], &jac1[0], k2, sz);
		vec_sum2 (&pt3[0], &pt3[0], &jac2[0], k1, sz);
		vec_sum2 (&pt3[0], &pt3[0], &jac3[0], k1, sz);
	}

	vec_mult_scalar (&pt3[0], &pt3[0], -1.0, sz);

	l_.calc1 (&pt2[0], &pt1[0], 0);
	vec_mult_scalar (&pt2[0], &pt2[0], (1. - theta) * mu, sz);

	vec_mult_scalar (&pt1[0], &pt1[0],
	                 1.0 / tau - (1. - theta) * sigma, sz);

	memset (u1, 0, sz * sizeof (double) );
	vec_sum (u1, u1, &pt1[0], sz);
	vec_sum (u1, u1, &pt2[0], sz);
	vec_sum (u1, u1, &pt3[0], sz);

	//@{test
	/*
	for (int i = 0; i < sz; ++i) {
		double x  = m_.ps[i].x();
		double y  = m_.ps[i].y();
		u1[i] += l_rp(x, y, t, sigma_, mu_);
	}
	*/
	//@}test

	{
		vector < double > tmp (rs);
		vector < double > rp (rs);
		u2p (&tmp[0], u1, m_);
		l_.idt_.mult_vector (&rp[0], &tmp[0]);
		if (bnd_w)
		{
			bnd_.mult_vector (&tmp[0], bnd_w);
			vec_sum (&rp[0], &rp[0], &tmp[0], (int) rp.size() );
		}
		solve (&u1[0], bnd_w, &rp[0], A_, m_);
	}
	l_.solve (u1, u1, bnd_u);
}

template < typename L, typename J, typename N >
void BarVortex < L, J, N >::calc_L_1 (double * u1,
                                   const double * u, const double * z,
                                   const double * bnd_u,
                                   const double * bnd_w,
                                   double t)
{
	int rs = (int) m_.inner.size(); // размерность внутренней области
	int sz = (int) m_.ps.size();   // размерность полная
	double k1 = conf_.k1;
	double k2 = conf_.k2;
	double theta = conf_.theta;
	double mu = conf_.mu;
	double sigma = conf_.sigma;
	double tau = conf_.tau;


	vector < double > z_lapl (sz);

	vector < double > pt1 (sz); //лаплас, умноженный на коэф
	vector < double > pt2 (sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3 (sz); //якобиан, умноженный на коэф

	l_.calc1 (&z_lapl[0], z, bnd_w);
	l_.calc1 (&pt1[0], u, bnd_w); //первая часть - лаплас, умноженный на коэф,

	{
		vector < double > jac1 (sz);
		vector < double > jac2 (sz);
		vector < double > jac3 (sz);

		j_.calc1 (&jac1[0], u, &lh_[0], 0);
		j_.calc1 (&jac2[0], z, &pt1[0], 0);
		j_.calc1 (&jac3[0], u, &z_lapl[0], 0);

		vec_sum2 (&pt3[0], &pt3[0], &jac1[0], k2, sz);
		vec_sum2 (&pt3[0], &pt3[0], &jac2[0], k1, sz);
		vec_sum2 (&pt3[0], &pt3[0], &jac3[0], k1, sz);
	}

	//vec_mult_scalar(&pt3[0], &pt3[0], -1.0, sz);

	l_.calc1 (&pt2[0], &pt1[0], 0);
	vec_mult_scalar (&pt2[0], &pt2[0], - theta * mu, sz);

	vec_mult_scalar (&pt1[0], &pt1[0],
	                 1.0 / tau + theta * sigma, sz);

	vec_sum (u1, u1, &pt1[0], sz);
	vec_sum (u1, u1, &pt2[0], sz);
	vec_sum (u1, u1, &pt3[0], sz);

	{
		vector < double > tmp (rs);
		vector < double > rp (rs);
		u2p (&tmp[0], u1, m_);
		l_.idt_.mult_vector (&rp[0], &tmp[0]);
		if (bnd_w)
		{
			bndb_.mult_vector (&tmp[0], bnd_w);
			vec_sum (&rp[0], &rp[0], &tmp[0], rs);
		}
		solve (&u1[0], bnd_w, &rp[0], Ab_, m_);
	}
	l_.solve (u1, u1, bnd_u);
}

template < typename L, typename J, typename N >
void BarVortex < L, J, N >::calc_LT (double * v1, const double * v, const double * z,
                                  const double * bnd_u,
                                  const double * bnd_w, double t)
{
	int rs = (int) m_.inner.size(); // размерность внутренней области
	int sz = (int) m_.ps.size();   // размерность полная
	double k1 = conf_.k1;
	double k2 = conf_.k2;
	double theta = conf_.theta;
	double mu = conf_.mu;
	double sigma = conf_.sigma;
	double tau = conf_.tau;

	vector < double > lz (sz);

	vector < double > pt1 (sz); //лаплас, умноженный на коэф
	vector < double > pt2 (sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3 (sz); //якобиан, умноженный на коэф

	l_.calc1 (&lz[0], z, bnd_w);
	l_.solve (&v1[0], v, bnd_u);

	{
		vector < double > tmp (rs);
		vector < double > rp (rs);
		u2p (&tmp[0], v1, m_);
		l_.idt_.mult_vector (&rp[0], &tmp[0]);
		if (bnd_w)
		{
			bnd_.mult_vector (&tmp[0], bnd_w);
			vec_sum (&rp[0], &rp[0], &tmp[0], rs);
		}
		solve (&v1[0], bnd_w, &rp[0], A_, m_);
	}

	l_.calc1 (&pt1[0], v1, bnd_w); //<- wtf ?
	vec_mult_scalar (&pt1[0], &pt1[0],
	                 1.0 / tau - (1.0 - theta) * sigma, sz);

	l_.calc1 (&pt2[0], v1, bnd_w); //<- wtf ?
	l_.calc1 (&pt2[0], &pt2[0], 0);

	vec_mult_scalar (&pt2[0], &pt2[0],
	                 (1.0 - theta) * mu, sz);

	{
		// JT
		double * h1 = &pt3[0];
		vector < double > p_lapl (sz);
		vector < double > tmp (sz);

		j_.calc1t (&tmp[0], v1, &lh_[0], 0);
		j_.calc1t (&p_lapl[0], z, v1, 0);
		l_.calc1 (&p_lapl[0], &p_lapl[0], 0);
		j_.calc1t (h1, v1, &lz[0], 0);

		vec_sum1 (h1, h1, &p_lapl[0], k1, k1, sz);
		vec_sum2 (h1, h1, &tmp[0], k2, sz);
		vec_mult_scalar (&h1[0], &h1[0], -1.0, sz);
	}

	memset (v1, 0, sz * sizeof (double) );
	vec_sum (v1, v1, &pt1[0], sz);
	vec_sum (v1, v1, &pt2[0], sz);
	vec_sum (v1, v1, &pt3[0], sz);
	set_bnd (v1, bnd_u, m_);
}

template < typename L, typename J, typename N >
void BarVortex < L, J, N >::S_step (double * Ans, const double * F)
{
	calc (Ans, F, 0, 0, 0);
}

/*
 * d L(psi)/dt + J(psi, L(z)) + J(z, L(psi)) + J(psi, l + h) + sigma L(psi) - mu LL(psi) = 0
 * L = Laplace
 */
template < typename L, typename J, typename N >
void BarVortex < L, J, N >::L_step (double * Ans, const double * F, const double * z)
{
	vector < double > tmp (m_.ps.size() );
	calc_L (&tmp[0], F, z, 0, 0, 0);
	memcpy (Ans, &tmp[0], tmp.size() * sizeof (double) );
}

template < typename L, typename J, typename N >
void BarVortex < L, J, N >::L_1_step (double * Ans, const double * F, const double * z)
{
	vector < double > tmp (m_.ps.size() );
	calc_L_1 (&tmp[0], F, z, 0, 0, 0);
	memcpy (Ans, &tmp[0], tmp.size() * sizeof (double) );
}

template < typename L, typename J, typename N >
void BarVortex < L, J, N >::LT_step (double * Ans, const double * F, const double * z)
{
	vector < double > tmp (m_.ps.size() );
	calc_LT (&tmp[0], F, z, 0, 0, 0);
	memcpy (Ans, &tmp[0], tmp.size() * sizeof (double) );
}

