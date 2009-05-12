/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (јлексей ќзерицкий)
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

#include "baroclin.h"
#include "util.h"

using namespace std;

double 
Baroclin::integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
                     Baroclin * d)
{
	double tau   = d->tau_;
	double mu    = d->mu_;
	double sigma = d->sigma_;

	double pt1, pt2;

	pt1  = integrate_cos(phi_i * phi_j, trk, m.ps);
	pt1 *= 1.0 / tau + sigma * 0.5;

	pt2  =  laplace(phi_j, phi_i, trk, m.ps);
	pt2 *= -0.5 * mu;

	return pt1 + pt2;
}

static double coriolis(double phi, double lambda)
{
	double omg = 0.0000727000000000;
	double l = omg * 2.0 * sin(phi);
	double h = cos(2.0 * lambda) * ipow(sin(2.0 * phi), 2);;
	return l + h;
}

Baroclin::Baroclin(const Mesh & m, rp_t rp, coriolis_t coriolis,
         double tau, double sigma, double mu, double sigma1, double mu1, double alpha)
	: m_(m), l_(m), j_(m), A_(m.inner.size()), tau_(tau), sigma_(sigma), mu_(mu),
	sigma1_(sigma1), mu1_(mu1), alpha_(alpha), rp_(rp), coriolis_(coriolis)
{
	lh_.resize(m_.ps.size());
	mke_proj(&lh_[0], m_, coriolis);

	/* ћатрица левой части совпадает с „афе-»нфантом на сфере */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, (integrate_cb_t)integrate_cb, this);
}

static double f(double x, double y, double t, 
				double mu, double sigma)
{
//	double a = exp(t) * sin(y) * sin(2.0 * x);
//	double b = -6.0 * exp(t) * sin(y) * sin(2.0 * x);

//	return a - mu * b + sigma * a;
	return 0.0;
}

struct Baroclin::right_part_cb_data
{
	const double * F;
	const double * bnd;
	Baroclin  * d;
};

double 
Baroclin::right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
                      right_part_cb_data * d)
{
	const double * F = d->F;
	double b;

	b = F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * Baroclin::integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, d->d);		
	}

	return b;
}

/**
 * (u1, u2) -> (u11, u21)
 * d L(u1)/dt + J(u1, L(u1) + l + h ?) + J(u2, L(u2)) + sigma/2 L(u1 - u2) - mu LL(u1) = 0
 * d L(u2)/dt + J(u1, L(u2)) + J(u2, L(u1) + l + h?) + sigma/2 L(u1 + u2) - mu LL(u2) -
 *   - alpha^2 (d u2/dt + J(u1, u2) - mu1 L(u2) + sigma1 u2 + f(phi, lambda))= 0
 * L = Laplace
 */
void Baroclin::calc(double * u11,  double * u21, 
		const double * u1, const double * u2, 
		const double * bnd, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полна€

	vector < double > omega_0(sz); // omega_0 = L (X_0)
	vector < double > omega_1(sz);
	vector < double > lomega(rs);  // L(omega)
	vector < double > rp(rs);

	vector < double > omega(rs);    // = omega_1 без границы
	vector < double > omega_lh(sz); // omega + l + h
	vector < double > jac(rs);      // jacobian

	vector < double > p1(sz);
	vector < double > prev_psi(sz);

	vector < double > X_0(sz);
	memcpy(&X_0[0], u1, sz * sizeof(double));

	// генерируем правую часть
	// w/dt + mu \Delta w / 2 - \sigma w/2 -
	// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)

	// omega = L (u)
	l_.calc1(&omega_0[0], &X_0[0], bnd);    //TODO: а чему у нас на кра€х равно omega?
	memcpy(&omega_1[0], &omega_0[0], sizeof(double) * sz);
	memcpy(&u11[0], &X_0[0], sizeof(double) * sz);


	// L (omega)
	l_.calc2(&lomega[0], &omega_1[0]);
	mke_u2p(&omega[0], &omega_1[0], m_);

	// w/dt + mu \Delta w / 2
	vec_sum1(&lomega[0], &omega[0], &lomega[0], 1.0 / tau_, mu_ * 0.5, rs);

	// w/dt + mu \Delta w / 2 - \sigma w/2
	vec_sum1(&lomega[0], &lomega[0], &omega[0], 1.0, -sigma_ * 0.5, rs);

	// в lomega содержитс€ права€ часть, котора€ не мен€етс€ при итераци€х!
	// права€ часть только на границе !

	for (int it = 0; it < 5; ++it) {
		// 0.5(w+w) + l + h <- дл€ вычислени€ якобиана это надо знать и на границе!
		vec_sum1(&omega_lh[0], &omega_1[0], &omega_0[0], 0.5, 0.5, sz);
		vec_sum(&omega_lh[0], &omega_lh[0], &lh_[0], sz);
		vec_sum1(&prev_psi[0], &X_0[0], &u11[0], 0.5, 0.5, sz);
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h)
		j_.calc2(&jac[0], &prev_psi[0], &omega_lh[0]);
		// w/dt + mu \Delta w / 2 - \sigma w/2 -
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)
#pragma omp parallel for
		for (int i = 0; i < rs; ++i) {
			int point = m_.inner[i];
			double x  = m_.ps[point].x();
			double y  = m_.ps[point].y();

			omega[i] = lomega[i] - jac[i] + f(x, y, t, mu_, sigma_);
		}

		// значени€ правой части на границе не знаем !
		right_part_cb_data data2;
		data2.F   = &omega[0];
		data2.bnd = bnd; //TODO: а чему у нас на кра€х равно omega?
		data2.d   = this;

		generate_right_part(&rp[0], m_, 
			(right_part_cb_t)right_part_cb, (void*)&data2);

		//TODO: тут граничное условие на омега!
		mke_solve(&omega_1[0], bnd, &rp[0], A_, m_);
		//TODO: а тут граничное условие на пси!
		memcpy(&prev_psi[0], u11, sz * sizeof(double));
		l_.solve(u11, &omega_1[0], bnd);
		{
			double nr = mke_dist(&prev_psi[0], &u11[0], m_, sphere_scalar_cb);
			//fprintf(stdout, "%le\n", nr);
			if (nr < 1e-8) {
				break;
			}
		}
	}
}

