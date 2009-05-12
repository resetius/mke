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

#include <assert.h>
#include <vector>
#include <math.h>

#include "barvortex.h"
#include "util.h"

using namespace std;

//#define SCHEME_THETA 0.5
#define SCHEME_THETA 1.0

static double 
integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
                     BarVortex * d)
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

static double 
integrate_backward_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
                     BarVortex * d)
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

BarVortex::BarVortex(const Mesh & m, rp_t rp, coriolis_t coriolis, double tau, 
		double sigma, double mu)
		 : m_(m), l_(m), j_(m), 
		 A_(m.inner.size()),
		 bnd_(m.inner.size()),
		 Ab_(m.inner.size()),
		 bndb_(m.inner.size()),
		 tau_(tau), sigma_(sigma), mu_(mu), theta_(SCHEME_THETA),
		 rp_(rp), coriolis_(coriolis)
{
	int sz = (int)m_.ps.size();
	lh_.resize(sz);
	mke_proj(&lh_[0], m_, coriolis);

	/* Матрица левой части совпадает с Чафе-Инфантом на сфере */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, (integrate_cb_t)integrate_cb, this);
	generate_boundary_matrix(bnd_, m_, (integrate_cb_t)integrate_cb, this);

	generate_matrix(Ab_, m_, (integrate_cb_t)integrate_backward_cb, this);
	generate_boundary_matrix(bndb_, m_, (integrate_cb_t)integrate_backward_cb, this);

	//f_.resize(sz);
	//for (int i = 0; i < sz; ++i) {
	//	double x = m_.ps[i].x();
	//	double y = m_.ps[i].y();

	//	f_[i] = rp_(x, y, mu_, sigma_);
	//}
}

struct right_part_cb_data
{
	const double * F;
	const double * bnd;
	BarVortex  * d;
};

double 
right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
                      right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

	//b = F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1 && d->bnd) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, d->d);		
	} else {
		b += F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);
	}

	return b;
}

double 
right_part_backward_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
                      right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

	//b = F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1 && d->bnd) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * integrate_backward_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, d->d);		
	} else {
		b += F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);
	}

	return b;
}

/**
 * d L(phi)/dt + J(phi, L(phi)) + J(phi, l + h) + sigma L(phi) - mu LL(phi) = f(phi, la)
 * L = Laplace
 */
void BarVortex::calc(double * psi, const double * x0, 
					 const double * bnd, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > omega_0(sz); // omega_0 = L (X_0)
	vector < double > omega_1(sz);
	vector < double > lomega(rs);  // L(omega)
	vector < double > rp(rs);

	vector < double > omega(rs);    // = omega_1 без границы
	vector < double > omega_lh(sz); // omega + l + h
	vector < double > jac(rs);      // jacobian

	vector < double > prev_psi(sz);

	vector < double > X_0(sz);
	memcpy(&X_0[0], x0, sz * sizeof(double));

	// генерируем правую часть
	// w/dt + mu \Delta w / 2 - \sigma w/2 -
	// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)

	// omega = L (u)
	l_.calc1(&omega_0[0], &X_0[0], bnd);    //TODO: а чему у нас на краях равно omega?
	memcpy(&omega_1[0], &omega_0[0], sizeof(double) * sz);
	memcpy(&psi[0], &X_0[0], sizeof(double) * sz);

	// L (omega)
	l_.calc2(&lomega[0], &omega_1[0]);
	mke_u2p(&omega[0], &omega_1[0], m_);

	// w/dt + mu \Delta w / 2
	vec_sum1(&lomega[0], &omega[0], &lomega[0], 1.0 / tau_, 
		mu_ * (1.0 - theta_), rs);

	// w/dt + mu \Delta w / 2 - \sigma w/2
	vec_sum1(&lomega[0], &lomega[0], &omega[0], 1.0, 
		-sigma_ * (1.0 - theta_), rs);

	// в lomega содержится правая часть, которая не меняется при итерациях!
	// правая часть только на границе !

	for (int it = 0; it < 20; ++it) {
		// 0.5(w+w) + l + h <- для вычисления Якобиана это надо знать и на границе!
		vec_sum1(&omega_lh[0], &omega_1[0], &omega_0[0], theta_, 
			1.0 - theta_, sz);
		vec_sum(&omega_lh[0], &omega_lh[0], &lh_[0], sz);
		vec_sum1(&prev_psi[0], &psi[0], &X_0[0], theta_, 
			1.0 - theta_, sz);
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h)
		j_.calc2(&jac[0], &prev_psi[0], &omega_lh[0]);
		// w/dt + mu \Delta w / 2 - \sigma w/2 -
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)
#pragma omp parallel for
		for (int i = 0; i < rs; ++i) {
			int point = m_.inner[i];
			double x  = m_.ps[point].x();
			double y  = m_.ps[point].y();

			omega[i] = lomega[i] - jac[i] + rp_(x, y, t, mu_, sigma_);
		}
#if 0
		right_part_cb_data data2;
		//генератор правой части учитывает то, что функция задана внутри!!!
		data2.F   = &omega[0];
		data2.bnd = bnd; //TODO: а чему у нас на краях равно omega?
		data2.d   = this;

		generate_right_part(&rp[0], m_, 
			(right_part_cb_t)right_part_cb, (void*)&data2);
#endif

		l_.idt_.mult_vector(&rp[0], &omega[0]);
		if (bnd) {
			// we use jac only as a storage !
			bnd_.mult_vector(&jac[0], bnd);
			vec_sum(&rp[0], &rp[0], &jac[0], rp.size());
		}

		//TODO: тут граничное условие на омега!
		mke_solve(&omega_1[0], bnd, &rp[0], A_, m_);
		//TODO: а тут граничное условие на пси!
		memcpy(&prev_psi[0], psi, sz * sizeof(double));
		l_.solve(psi, &omega_1[0], bnd);
		{
			double nr = mke_dist(&prev_psi[0], &psi[0], m_, sphere_scalar_cb);
			//fprintf(stdout, "%le\n", nr);
			if (nr < 1e-8) {
				break;
			}
		}
	}
}

/**
 * d L(phi)/dt + J(phi, L(z)) + J(z, L(phi)) + J(phi, l + h) + sigma L(phi) - mu LL(phi) = 0
 * L = Laplace
 */
#if 0
void BarVortex::calc_L(double * psi, const double * x0, const double * z,
					 const double * bnd, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > lz(sz);
	vector < double > lz1(sz);
	vector < double > omega_0(sz); // omega_0 = L (X_0)
	vector < double > omega_1(sz);
	vector < double > lomega(rs);  // L(omega)
	vector < double > rp(rs);

	vector < double > omega(rs);    // = omega_1 без границы
	vector < double > omega_lh(sz); // omega + l + h
	vector < double > jac1(rs);     // jacobian
	vector < double > jac2(rs);     // jacobian

	vector < double > prev_psi(sz);

	vector < double > X_0(sz);
	memcpy(&X_0[0], x0, sz * sizeof(double));

	// генерируем правую часть
	// w/dt + mu \Delta w / 2 - \sigma w/2 -
	// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)

	// omega = L (u)
	l_.calc1(&omega_0[0], &X_0[0], bnd);    //TODO: а чему у нас на краях равно omega?
	l_.calc1(&lz[0], &z[0], bnd);

	memcpy(&omega_1[0], &omega_0[0], sizeof(double) * sz);
	memcpy(&psi[0], &X_0[0], sizeof(double) * sz);

	// L (omega)
	l_.calc2(&lomega[0], &omega_1[0]);
	mke_u2p(&omega[0], &omega_1[0], m_);

	// w/dt + mu \Delta w / 2
	vector_sum1(&lomega[0], &omega[0], &lomega[0], 1.0 / tau_, 
		mu_ * (1.0 - theta_), rs);

	// w/dt + mu \Delta w / 2 - \sigma w/2
	vector_sum1(&lomega[0], &lomega[0], &omega[0], 1.0, 
		-sigma_ * (1.0 - theta_), rs);

	// в lomega содержится правая часть, которая не меняется при итерациях!
	// правая часть только на границе !

	for (int it = 0; it < 1; ++it) {
		// J(0.5(u+u), L(z)+l+h) + J(z, 0.5(w+w))
		// 0.5(w+w) <- для вычисления Якобиана это надо знать и на границе!
		vector_sum1(&omega_lh[0], &omega_1[0], &omega_0[0], 
			theta_, 1.0 - theta_, sz);
		//L(z)+l+h
		mke_vector_sum(&lz1[0], &lz[0], &lh_[0], sz);
		//0.5(u+u)
		vector_sum1(&prev_psi[0], &psi[0], &X_0[0], 
			theta_, 1.0 - theta_, sz);

		// J(0.5(u+u), L(z)+l+h)
		j_.calc2(&jac1[0], &prev_psi[0], &lz1[0]);
		// J(z, 0.5(w+w))
		j_.calc2(&jac2[0], &z[0], &omega_lh[0]);
		// w/dt + mu \Delta w / 2 - \sigma w/2 -
		// - J(0.5(u+u), L(z)+l+h) - J(z, 0.5(w+w))
#pragma omp parallel for
		for (int i = 0; i < rs; ++i) {
			int point = m_.inner[i];
			double x  = m_.ps[point].x();
			double y  = m_.ps[point].y();

			omega[i] = lomega[i] - jac1[i] - jac2[i];
		}

		l_.idt_.mult_vector(&rp[0], &omega[0]);
		if (bnd) {
			// we use jac1 only as a storage !
			bnd_.mult_vector(&jac1[0], bnd);
			mke_vector_sum(&rp[0], &rp[0], &jac1[0], rp.size());
		}

		//TODO: тут граничное условие на омега!
		mke_solve(&omega_1[0], bnd, &rp[0], A_, m_);
		//TODO: а тут граничное условие на пси!
		memcpy(&prev_psi[0], psi, sz * sizeof(double));
		l_.solve(psi, &omega_1[0], bnd);
		{
			double nr = mke_dist(&prev_psi[0], &psi[0], m_, sphere_scalar_cb);
			//fprintf(stdout, "%le\n", nr);
			if (nr < 1e-5) {
				break;
			}
		}
	}
}
#endif

void BarVortex::calc_L(double * u1, const double * u, const double * z,
					   const double * bnd, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > z_lapl(sz);

	vector < double > pt1(sz); //лаплас, умноженный на коэф
	vector < double > pt2(sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3(sz); //якобиан, умноженный на коэф

	l_.calc1(&z_lapl[0], z, bnd);
	l_.calc1(&pt1[0], u, bnd); //первая часть - лаплас, умноженный на коэф, 

	{
		vector < double > jac1(sz);
		vector < double > jac2(sz);
		vector < double > jac3(sz);

		j_.calc1(&jac1[0], u, &lh_[0], bnd);
		j_.calc1(&jac2[0], z, &pt1[0], bnd);
		j_.calc1(&jac3[0], u, &z_lapl[0], bnd);

		vec_sum(&pt3[0], &pt3[0], &jac1[0], sz);
		vec_sum(&pt3[0], &pt3[0], &jac2[0], sz);
		vec_sum(&pt3[0], &pt3[0], &jac3[0], sz);
	}

	vec_mult_scalar(&pt3[0], &pt3[0], -1.0, sz);

	l_.calc1(&pt2[0], &pt1[0], bnd);
	vec_mult_scalar(&pt2[0], &pt2[0], (1. - theta_) * mu_, sz);

	vec_mult_scalar(&pt1[0], &pt1[0], 
		1.0 / tau_ - (1. - theta_) * sigma_, sz);

	vec_sum(u1, u1, &pt1[0], sz);
	vec_sum(u1, u1, &pt2[0], sz);
	vec_sum(u1, u1, &pt3[0], sz);

	{
		vector < double > tmp(rs);
		vector < double > rp(rs);
		mke_u2p(&tmp[0], u1, m_);
		l_.idt_.mult_vector(&rp[0], &tmp[0]);
		if (bnd) {
			bnd_.mult_vector(&tmp[0], bnd);
			vec_sum(&rp[0], &rp[0], &tmp[0], rp.size());
		}
		mke_solve(&u1[0], bnd, &rp[0], A_, m_);
	}
	l_.solve(u1, u1, bnd);
}

void BarVortex::calc_L_1(double * u1, const double * u, const double * z,
					   const double * bnd, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > z_lapl(sz);

	vector < double > pt1(sz); //лаплас, умноженный на коэф
	vector < double > pt2(sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3(sz); //якобиан, умноженный на коэф

	l_.calc1(&z_lapl[0], z, bnd);
	l_.calc1(&pt1[0], u, bnd); //первая часть - лаплас, умноженный на коэф, 

	{
		vector < double > jac1(sz);
		vector < double > jac2(sz);
		vector < double > jac3(sz);

		j_.calc1(&jac1[0], u, &lh_[0], bnd);
		j_.calc1(&jac2[0], z, &pt1[0], bnd);
		j_.calc1(&jac3[0], u, &z_lapl[0], bnd);

		vec_sum(&pt3[0], &pt3[0], &jac1[0], sz);
		vec_sum(&pt3[0], &pt3[0], &jac2[0], sz);
		vec_sum(&pt3[0], &pt3[0], &jac3[0], sz);
	}

	//vec_mult_scalar(&pt3[0], &pt3[0], -1.0, sz);

	l_.calc1(&pt2[0], &pt1[0], bnd);
	vec_mult_scalar(&pt2[0], &pt2[0], - theta_ * mu_, sz);

	vec_mult_scalar(&pt1[0], &pt1[0], 
		1.0 / tau_ + theta_ * sigma_, sz);

	vec_sum(u1, u1, &pt1[0], sz);
	vec_sum(u1, u1, &pt2[0], sz);
	vec_sum(u1, u1, &pt3[0], sz);

	{
		vector < double > tmp(rs);
		vector < double > rp(rs);
		mke_u2p(&tmp[0], u1, m_);
		l_.idt_.mult_vector(&rp[0], &tmp[0]);
		if (bnd) {
			bndb_.mult_vector(&tmp[0], bnd);
			vec_sum(&rp[0], &rp[0], &tmp[0], rp.size());
		}
		mke_solve(&u1[0], bnd, &rp[0], Ab_, m_);
	}
	l_.solve(u1, u1, bnd);
}

void BarVortex::calc_LT(double * v1, const double * v, const double * z, const double * bnd, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > lz(sz);

	vector < double > pt1(sz); //лаплас, умноженный на коэф
	vector < double > pt2(sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3(sz); //якобиан, умноженный на коэф

	l_.calc1(&lz[0], z, bnd);
	l_.solve(&v1[0], v, bnd);

	{
		vector < double > tmp(rs);
		vector < double > rp(rs);
		mke_u2p(&tmp[0], v1, m_);
		l_.idt_.mult_vector(&rp[0], &tmp[0]);
		if (bnd) {
			bnd_.mult_vector(&tmp[0], bnd);
			vec_sum(&rp[0], &rp[0], &tmp[0], rp.size());
		}
		mke_solve(&v1[0], bnd, &rp[0], A_, m_);
	}

	l_.calc1(&pt1[0], v1, bnd);
	vec_mult_scalar(&pt1[0], &pt1[0], 
		1.0 / tau_ - (1.0 - theta_) * sigma_, sz);

	l_.calc1(&pt2[0], v1, bnd);
	l_.calc1(&pt2[0], &pt2[0], bnd);

	vec_mult_scalar(&pt2[0], &pt2[0], 
		(1.0 - theta_) * mu_, sz);

	{
		// JT
		double * h1 = &pt3[0];
		vector < double > p_lapl(sz);
		vector < double > tmp(sz);

		j_.calc1t(&tmp[0], v1, &lh_[0], bnd);
		j_.calc1t(&p_lapl[0], z, v1, bnd);
		l_.calc1(&p_lapl[0], &p_lapl[0], bnd);
		j_.calc1t(h1, v1, &lz[0], bnd);

		vec_sum(h1, h1, &p_lapl[0], sz);
		vec_sum(h1, h1, &tmp[0], sz);
		vec_mult_scalar(&h1[0], &h1[0], -1.0, sz);
	}

	memset(v1, 0, sz * sizeof(double));
	vec_sum(v1, v1, &pt1[0], sz);
	vec_sum(v1, v1, &pt2[0], sz);
	vec_sum(v1, v1, &pt3[0], sz);
}

void BarVortex::S_step(double * Ans, const double * F)
{
	calc(Ans, F, 0, 0);
}

/* J(phi, L(z)) + J(z, L(phi)) + J(phi, l + h) + sigma L(phi) - mu LL(phi) */
void BarVortex::L_spectr(double * u1, const double * u, const double * z, const double * bnd)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	vector < double > z_lapl(sz);

	vector < double > pt1(sz); //лаплас, умноженный на коэф
	vector < double > pt2(sz); //лаплас в квадрате, умноженный на коэф
	vector < double > pt3(sz); //якобиан, умноженный на коэф

	l_.calc1(&z_lapl[0], z, bnd);
	l_.calc1(&pt1[0], u, bnd); //первая часть - лаплас, умноженный на коэф, 

	{
		vector < double > jac1(sz);
		vector < double > jac2(sz);
		vector < double > jac3(sz);

		j_.calc1(&jac1[0], u, &lh_[0], bnd);
		j_.calc1(&jac2[0], z, &pt1[0], bnd);
		j_.calc1(&jac3[0], u, &z_lapl[0], bnd);

		vec_sum(&pt3[0], &pt3[0], &jac1[0], sz);
		vec_sum(&pt3[0], &pt3[0], &jac2[0], sz);
		vec_sum(&pt3[0], &pt3[0], &jac3[0], sz);
	}

	vec_mult_scalar(&pt3[0], &pt3[0], -1.0, sz);

	l_.calc1(&pt2[0], &pt1[0], bnd);
	vec_mult_scalar(&pt2[0], &pt2[0], mu_, sz);

	vec_mult_scalar(&pt1[0], &pt1[0], - sigma_, sz);

	vec_sum(u1, u1, &pt1[0], sz);
	vec_sum(u1, u1, &pt2[0], sz);
	vec_sum(u1, u1, &pt3[0], sz);
}

/**
 * d L(phi)/dt + J(phi, L(z)) + J(z, L(phi)) + J(phi, l + h) + sigma L(phi) - mu LL(phi) = 0
 * L = Laplace
 */
void BarVortex::L_step(double * Ans, const double * F, const double * z)
{
	vector < double > tmp(m_.ps.size());
	calc_L(&tmp[0], F, z, 0, 0);
	memcpy(Ans, &tmp[0], tmp.size() * sizeof(double));
}

void BarVortex::L_1_step(double * Ans, const double * F, const double * z)
{
	vector < double > tmp(m_.ps.size());
	calc_L_1(&tmp[0], F, z, 0, 0);
	memcpy(Ans, &tmp[0], tmp.size() * sizeof(double));
}

void BarVortex::LT_step(double * Ans, const double * F, const double * z)
{
	vector < double > tmp(m_.ps.size());
	calc_LT(&tmp[0], F, z, 0, 0);
	memcpy(Ans, &tmp[0], tmp.size() * sizeof(double));
}

