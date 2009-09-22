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
#include <vector>
#include <math.h>

#include "barvortex.h"
#include "util.h"
#include "ver.h"

VERSION("$Id$");

using namespace std;
using namespace phelm;

#define SCHEME_THETA 0.5
//#define SCHEME_THETA 1.0
//#define SCHEME_THETA 0.2

namespace bv_private {

template < typename BV > 
double 
integrate_cb( const Polynom & phi_i,
              const Polynom & phi_j, 
              const Triangle & trk,
              const Mesh & m,
              int point_i, int point_j,
			  int i, int j,
              BV * d)
{
	double tau   = d->tau_;
	double mu    = d->mu_;
	double sigma = d->sigma_;

	double pt1, pt2;

	pt1  = integrate_cos(phi_i * phi_j, trk, m.ps);
	pt1 *= 1.0 / tau + sigma * d->theta_;

	pt2  =  slaplace(phi_j, phi_i, trk, m.ps);
	pt2 *= - d->theta_ * mu;

	return pt1 + pt2;
}

template < typename BV > 
double 
integrate_backward_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
					 int i, int j,
                     BV * d)
{
	double tau   = d->tau_;
	double mu    = d->mu_;
	double sigma = d->sigma_;

	double pt1, pt2;

	pt1  = integrate_cos(phi_i * phi_j, trk, m.ps);
	pt1 *= 1.0 / tau - sigma * (1.0 - d->theta_);

	pt2  =  slaplace(phi_j, phi_i, trk, m.ps);
	pt2 *= (1.0 - d->theta_) * mu;

	return pt1 + pt2;
}

template < typename BV > 
struct right_part_cb_data
{
	const double * F;
	const double * bnd;
	BV  * d;
};

template < typename BV > 
double 
right_part_cb( const Polynom & phi_i,
              const Polynom & phi_j,
              const Triangle & trk,
              const Mesh & m,
              int point_i, int point_j,
			  int i, int j,
			   right_part_cb_data < BV >  * d)
{
	const double * F = d->F;
	double b = 0.0;

	//b = F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1 && d->bnd) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, i, j, d->d);
	} else if (m.ps_flags[point_j] == 0) {
		b += F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);
	}

	return b;
}

template < typename BV > 
double 
right_part_backward_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
					  int i, int j,
						right_part_cb_data < BV > * d)
{
	const double * F = d->F;
	double b = 0.0;

	//b = F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1 && d->bnd) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * integrate_backward_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, i, j, d->d);
	} else if (m.ps_flags[point_j] == 0) {
		b += F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);
	}

	return b;
}
}

template < typename L, typename J >
BarVortex < L, J > ::BarVortex(const Mesh & m, rp_t rp, coriolis_t coriolis, double tau, 
					 double sigma, double mu, double k1, double k2)
	: SphereNorm < double > (m), m_(m), l_(m), j_(m), 
	  A_((int)m.inner.size()),
	  bnd_((int)m.inner.size()),
	  Ab_((int)m.inner.size()),
	  bndb_((int)m.inner.size()),
	  tau_(tau), sigma_(sigma), mu_(mu), k1_(k1), k2_(k2), theta_(SCHEME_THETA),
	  rp_(rp), coriolis_(coriolis)
{
	int sz = (int)m_.ps.size();
	lh_.resize(sz);
	proj(&lh_[0], m_, coriolis);

	/* Матрица левой части совпадает с Чафе-Инфантом на сфере */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, bv_private::integrate_cb < BarVortex < L, J > > , this);
	generate_boundary_matrix(bnd_, m_, bv_private::integrate_cb < BarVortex < L, J > >, this);

	generate_matrix(Ab_, m_, bv_private::integrate_backward_cb < BarVortex < L, J > > , this);
	generate_boundary_matrix(bndb_, m_, bv_private::integrate_backward_cb < BarVortex < L, J > >, this);
}

template < typename L, typename J >
void BarVortex < L, J > ::info()
{
	fprintf(stderr, "#tau:%.16lf\n", tau_);
	fprintf(stderr, "#sigma:%.16lf\n", sigma_);
	fprintf(stderr, "#mu:%.16lf\n", mu_);
	fprintf(stderr, "#k1:%.16lf\n", k1_);
	fprintf(stderr, "#k2:%.16lf\n", k2_);
	fprintf(stderr, "#theta:%.16lf\n", theta_);
}

/*
 * d L(phi)/dt + k1 J(psi, L(psi)) + k2 J(psi, l + h) + sigma L(psi) - mu LL(psi) = f(phi, la)
 * L = Laplace
 */
template < typename L, typename J >
void BarVortex < L, J > ::calc(double * u1,
							   const double * u, 
							   const double * bnd_u,
							   const double * bnd_w,
							   double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	double nr0 = norm(u);

	vector < double > w(sz);       // w = L(u)
	vector < double > dw(sz);      // dw = L(w) = LL(u)
	vector < double > FC(sz);      // правая часть

	// next
	vector < double > u_n(sz);
	vector < double > w_n(sz);
	vector < double > u_n1(sz);

	// tmp
	vector < double > tmp1(sz);
	vector < double > tmp2(sz);

	// inner (!!!) jacobian
	vector < double > jac(rs);

	//
	vector < double > F(rs);
	vector < double > rp(rs);

	// генерируем правую часть
	// w/dt + mu (1-theta) L w - \sigma(1-theta) w -
	// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h) + f(x, y)

	// w = L (u)
	l_.calc1(&w[0], &u[0], bnd_w);

	// dw = L (w)
	l_.calc1(&dw[0], &w[0], 0);

	// w/dt + mu (1-theta) L w
	vec_sum1(&FC[0], &w[0], &dw[0], 1.0 / tau_, 
		mu_ * (1.0 - theta_), sz);

	// w/dt + mu (1-theta) L w - \sigma (1-theta) w
	vec_sum1(&FC[0], &FC[0], &w[0], 1.0, 
		-sigma_ * (1.0 - theta_), sz);

#pragma omp parallel for
	for (int i = 0; i < sz; ++i) {
		int point = i;
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();

		if (rp_) {
			FC[i] += rp_(x, y, t, mu_, sigma_);
		}
	}

	memcpy(&u_n[0], &u[0], sz * sizeof(double));
	memcpy(&w_n[0], &w[0], sz * sizeof(double));

	// в FC содержится правая часть, которая не меняется при итерациях!

	for (int it = 0; it < 1000; ++it) {
		//   k1 J(0.5(u+u), 0.5(w+w)) + k2 J(0.5(u+u), l + h)   =
		// = J(0.5 (u+u), 0.5 k1 (w+w)) + J(0.5 (u+u), k2 (l + h)) =
		// = J(0.5 (u+u), 0.5 k1 (w+w) + k2 (l + h))

		// 0.5 k1 (w+w) + k2 (l + h) <- для вычисления Якобиана это надо знать и на границе!
		vec_sum1(&tmp1[0], &w_n[0], &w[0], k1_ * theta_, 
				 k1_ * (1.0 - theta_), sz);
		vec_sum2(&tmp1[0], &tmp1[0], &lh_[0], k2_, sz);
		// 0.5(u+u)
		vec_sum1(&tmp2[0], &u_n[0], &u[0], theta_, 
			1.0 - theta_, sz);
		// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h)
		j_.calc2(&jac[0], &tmp2[0], &tmp1[0]);
		// w/dt + mu (1-theta) L w  - \sigma (1-theta) w -
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)
#pragma omp parallel for
		for (int i = 0; i < rs; ++i) {
			int point = m_.inner[i];
			F[i] = FC[point] - jac[i];
		}
#if 0
		right_part_cb_data data2;
		//генератор правой части учитывает то, что функция задана внутри!!!
		data2.F   = &F[0];
		data2.bnd = bnd; 
		data2.d   = this;

		generate_right_part(&rp[0], m_, 
			right_part_cb, &data2);
#endif
#if 1
		memset(&rp[0], 0, rs * sizeof(double));
		l_.idt_.mult_vector(&rp[0], &F[0]);
		if (bnd_w) {
			memset(&tmp1[0], 0, rs * sizeof(double));
			bnd_.mult_vector(&tmp1[0], bnd_w);
			vec_sum(&rp[0], &rp[0], &tmp1[0], (int)rp.size());
		}
#endif
		// тут граничное условие на омега!
		solve(&w_n[0], bnd_w, &rp[0], A_, m_);
		// а тут граничное условие на пси!
		l_.solve(&u_n1[0], &w_n[0], bnd_u);
		
		//l_.solve(&u1[0], &w_n[0], bnd);
		//phelm::smooth1(&u_n[0], &u1[0], m_);

		double nr = dist(&u_n1[0], &u_n[0]);
		u_n1.swap(u_n);
		if (nr / nr0 < 1e-14) {
			break;
		}
	}
	memcpy(u1, &u_n[0], sz * sizeof(double));
	//phelm::smooth1(&u1[0], &u_n[0], m_);
}

double l_rp(double x, double y, double t, double sigma, double mu)
{
	return -9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2)-20*x*cos(y+t)*
		ipow(cos(x),4)-9*sigma*sin(y+t)*
		ipow(cos(x),3)*sin(x)+15*sigma*sin(y+t)*
		ipow(cos(x),2)*x-20*sigma*sin(y+t)*
		ipow(cos(x),4)*x-360*mu*sin(y+t)*sin(x)*
		ipow(cos(x),3)+390*mu*sin(y+t)*x*
		ipow(cos(x),2)-400*mu*sin(y+t)*x*
		ipow(cos(x),4)+147*mu*sin(y+t)*sin(x)*cos(x)-45*mu*sin(y+t)*x;
}

template < typename L, typename J >
void BarVortex < L, J >::calc_L(double * u1, const double * u, const double * z,
								const double * bnd_u,
								const double * bnd_w,
								double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > z_lapl(sz);

	vector < double > pt1(sz); //лаплас, умноженный на коэф
	vector < double > pt2(sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3(sz); //якобиан, умноженный на коэф

	l_.calc1(&z_lapl[0], z, bnd_w);
	l_.calc1(&pt1[0], u, bnd_w); //первая часть - лаплас, умноженный на коэф, 

	{
		vector < double > jac1(sz);
		vector < double > jac2(sz);
		vector < double > jac3(sz);

		j_.calc1(&jac1[0], u, &lh_[0], 0);
		j_.calc1(&jac2[0], z, &pt1[0], 0);
		j_.calc1(&jac3[0], u, &z_lapl[0], 0);

		vec_sum2(&pt3[0], &pt3[0], &jac1[0], k2_, sz);
		vec_sum2(&pt3[0], &pt3[0], &jac2[0], k1_, sz);
		vec_sum2(&pt3[0], &pt3[0], &jac3[0], k1_, sz);
	}

	vec_mult_scalar(&pt3[0], &pt3[0], -1.0, sz);

	l_.calc1(&pt2[0], &pt1[0], 0);
	vec_mult_scalar(&pt2[0], &pt2[0], (1. - theta_) * mu_, sz);

	vec_mult_scalar(&pt1[0], &pt1[0], 
		1.0 / tau_ - (1. - theta_) * sigma_, sz);

	memset(u1, 0, sz * sizeof(double));
	vec_sum(u1, u1, &pt1[0], sz);
	vec_sum(u1, u1, &pt2[0], sz);
	vec_sum(u1, u1, &pt3[0], sz);

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
		vector < double > tmp(rs);
		vector < double > rp(rs);
		u2p(&tmp[0], u1, m_);
		l_.idt_.mult_vector(&rp[0], &tmp[0]);
		if (bnd_w) {
			bnd_.mult_vector(&tmp[0], bnd_w);
			vec_sum(&rp[0], &rp[0], &tmp[0], (int)rp.size());
		}
		solve(&u1[0], bnd_w, &rp[0], A_, m_);
	}
	l_.solve(u1, u1, bnd_u);
}

template < typename L, typename J >
void BarVortex < L, J >::calc_L_1(double * u1,
								  const double * u, const double * z,
								  const double * bnd_u,
								  const double * bnd_w,
								  double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > z_lapl(sz);

	vector < double > pt1(sz); //лаплас, умноженный на коэф
	vector < double > pt2(sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3(sz); //якобиан, умноженный на коэф

	l_.calc1(&z_lapl[0], z, bnd_w);
	l_.calc1(&pt1[0], u, bnd_w); //первая часть - лаплас, умноженный на коэф, 

	{
		vector < double > jac1(sz);
		vector < double > jac2(sz);
		vector < double > jac3(sz);

		j_.calc1(&jac1[0], u, &lh_[0], 0);
		j_.calc1(&jac2[0], z, &pt1[0], 0);
		j_.calc1(&jac3[0], u, &z_lapl[0], 0);

		vec_sum2(&pt3[0], &pt3[0], &jac1[0], k2_, sz);
		vec_sum2(&pt3[0], &pt3[0], &jac2[0], k1_, sz);
		vec_sum2(&pt3[0], &pt3[0], &jac3[0], k1_, sz);
	}

	//vec_mult_scalar(&pt3[0], &pt3[0], -1.0, sz);

	l_.calc1(&pt2[0], &pt1[0], 0);
	vec_mult_scalar(&pt2[0], &pt2[0], - theta_ * mu_, sz);

	vec_mult_scalar(&pt1[0], &pt1[0], 
		1.0 / tau_ + theta_ * sigma_, sz);

	vec_sum(u1, u1, &pt1[0], sz);
	vec_sum(u1, u1, &pt2[0], sz);
	vec_sum(u1, u1, &pt3[0], sz);

	{
		vector < double > tmp(rs);
		vector < double > rp(rs);
		u2p(&tmp[0], u1, m_);
		l_.idt_.mult_vector(&rp[0], &tmp[0]);
		if (bnd_w) {
			bndb_.mult_vector(&tmp[0], bnd_w);
			vec_sum(&rp[0], &rp[0], &tmp[0], rs);
		}
		solve(&u1[0], bnd_w, &rp[0], Ab_, m_);
	}
	l_.solve(u1, u1, bnd_u);
}

template < typename L, typename J >
void BarVortex < L, J >::calc_LT(double * v1, const double * v, const double * z,
								 const double * bnd_u,
								 const double * bnd_w, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > lz(sz);

	vector < double > pt1(sz); //лаплас, умноженный на коэф
	vector < double > pt2(sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3(sz); //якобиан, умноженный на коэф

	l_.calc1(&lz[0], z, bnd_w);
	l_.solve(&v1[0], v, bnd_u);

	{
		vector < double > tmp(rs);
		vector < double > rp(rs);
		u2p(&tmp[0], v1, m_);
		l_.idt_.mult_vector(&rp[0], &tmp[0]);
		if (bnd_w) {
			bnd_.mult_vector(&tmp[0], bnd_w);
			vec_sum(&rp[0], &rp[0], &tmp[0], rs);
		}
		solve(&v1[0], bnd_w, &rp[0], A_, m_);
	}

	l_.calc1(&pt1[0], v1, bnd_w); //<- wtf ?
	vec_mult_scalar(&pt1[0], &pt1[0], 
		1.0 / tau_ - (1.0 - theta_) * sigma_, sz);

	l_.calc1(&pt2[0], v1, bnd_w); //<- wtf ?
	l_.calc1(&pt2[0], &pt2[0], 0);

	vec_mult_scalar(&pt2[0], &pt2[0], 
		(1.0 - theta_) * mu_, sz);

	{
		// JT
		double * h1 = &pt3[0];
		vector < double > p_lapl(sz);
		vector < double > tmp(sz);

		j_.calc1t(&tmp[0], v1, &lh_[0], 0);
		j_.calc1t(&p_lapl[0], z, v1, 0);
		l_.calc1(&p_lapl[0], &p_lapl[0], 0);
		j_.calc1t(h1, v1, &lz[0], 0);

		vec_sum1(h1, h1, &p_lapl[0], k1_, k1_, sz);
		vec_sum2(h1, h1, &tmp[0], k2_, sz);
		vec_mult_scalar(&h1[0], &h1[0], -1.0, sz);
	}

	memset(v1, 0, sz * sizeof(double));
	vec_sum(v1, v1, &pt1[0], sz);
	vec_sum(v1, v1, &pt2[0], sz);
	vec_sum(v1, v1, &pt3[0], sz);
	set_bnd(v1, bnd_u, m_);
}

template < typename L, typename J >
void BarVortex < L, J >::S_step(double * Ans, const double * F)
{
	calc(Ans, F, 0, 0, 0);
}

/*
 * d L(psi)/dt + J(psi, L(z)) + J(z, L(psi)) + J(psi, l + h) + sigma L(psi) - mu LL(psi) = 0
 * L = Laplace
 */
template < typename L, typename J >
void BarVortex < L, J >::L_step(double * Ans, const double * F, const double * z)
{
	vector < double > tmp(m_.ps.size());
	calc_L(&tmp[0], F, z, 0, 0, 0);
	memcpy(Ans, &tmp[0], tmp.size() * sizeof(double));
}

template < typename L, typename J >
void BarVortex < L, J >::L_1_step(double * Ans, const double * F, const double * z)
{
	vector < double > tmp(m_.ps.size());
	calc_L_1(&tmp[0], F, z, 0, 0, 0);
	memcpy(Ans, &tmp[0], tmp.size() * sizeof(double));
}

template < typename L, typename J >
void BarVortex < L, J >::LT_step(double * Ans, const double * F, const double * z)
{
	vector < double > tmp(m_.ps.size());
	calc_LT(&tmp[0], F, z, 0, 0);
	memcpy(Ans, &tmp[0], tmp.size() * sizeof(double));
}
