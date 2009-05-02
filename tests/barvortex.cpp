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

#include "barvortex.h"
#include "util.h"

using namespace std;

double 
BarVortex::integrate_cb( const Polynom & phi_i,
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
	pt1 *= 1.0 / tau + sigma * 0.5;

	pt2  =  laplace(phi_j, phi_i, trk, m.ps);
	pt2 *= -0.5 * mu;

	return pt1 + pt2;
}

BarVortex::BarVortex(const Mesh & m, rp_t rp, coriolis_t coriolis, double tau, 
		double sigma, double mu)
		 : m_(m), l_(m), j_(m), A_(m.inner.size()),
		 tau_(tau), sigma_(sigma), mu_(mu), rp_(rp), coriolis_(coriolis)
{
	int sz = (int)m_.ps.size();
	lh_.resize(sz);
	mke_proj(&lh_[0], m_, coriolis);

	/* ћатрица левой части совпадает с „афе-»нфантом на сфере */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, (integrate_cb_t)integrate_cb, this);

	//f_.resize(sz);
	//for (int i = 0; i < sz; ++i) {
	//	double x = m_.ps[i].x();
	//	double y = m_.ps[i].y();

	//	f_[i] = rp_(x, y, mu_, sigma_);
	//}
	//l_.calc1(&f_[0], &f_[0], 0); // < !!!!
}

struct BarVortex::right_part_cb_data
{
	const double * F;
	const double * bnd;
	BarVortex  * d;
};

double 
BarVortex::right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
                      right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

	b = F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * BarVortex::integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, d->d);		
	}

	return b;
}

#if 0
void BarVortex::calc(double * Ans, const double * x0, 
					 const double * bnd, double t)
{

	int rs  = (int)m_.inner.size();
	int sz  = (int)m_.ps.size();
	vector < double > u(rs);
	vector < double > p(sz);
	vector < double > delta_u(rs);
	vector < double > rp(rs);

	vector < double > X0(sz);

	SphereLaplace & laplace_ = l_;

	// генерируем правую часть
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	laplace_.calc1(&X0[0], x0, bnd);

	mke_u2p(&u[0], &X0[0], m_);
	laplace_.calc2(&delta_u[0], &X0[0]);

	// u/dt + mu \Delta u / 2
	vector_sum1(&delta_u[0], &u[0], &delta_u[0], 1.0 / tau_, mu_ * 0.5, rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vector_sum1(&delta_u[0], &delta_u[0], &u[0], 1.0, -sigma_ * 0.5, rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		int point = m_.inner[i];
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();

		u[i] = delta_u[i] + rp_(x, y, t, mu_, sigma_);
	}

	// правую часть на границе не знаем !!!
	mke_p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	right_part_cb_data data2;
	data2.F   = &p[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		(right_part_cb_t)right_part_cb, (void*)&data2);

	vector < double > omg(sz);
	mke_solve(&omg[0], bnd, &rp[0], A_, m_);
	laplace_.solve(Ans, &omg[0], bnd);
}
#endif

/**
 * d L(phi)/dt + J(phi, L(phi)) + J(phi, l + h) + sigma L(phi) - mu LL(phi) = f(phi, la)
 * L = Laplace
 */
#if 1
void BarVortex::calc(double * psi, const double * x0, 
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

	vector < double > prev_psi(sz);

	vector < double > X_0(sz);
	memcpy(&X_0[0], x0, sz * sizeof(double));

	// генерируем правую часть
	// w/dt + mu \Delta w / 2 - \sigma w/2 -
	// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)

	// omega = L (u)
	l_.calc1(&omega_0[0], &X_0[0], bnd);    //TODO: а чему у нас на кра€х равно omega?
	memcpy(&omega_1[0], &omega_0[0], sizeof(double) * sz);
	memcpy(&psi[0], &X_0[0], sizeof(double) * sz);

	// L (omega)
	l_.calc2(&lomega[0], &omega_1[0]);
	mke_u2p(&omega[0], &omega_1[0], m_);

	// w/dt + mu \Delta w / 2
	vector_sum1(&lomega[0], &omega[0], &lomega[0], 1.0 / tau_, mu_ * 0.5, rs);

	// w/dt + mu \Delta w / 2 - \sigma w/2
	vector_sum1(&lomega[0], &lomega[0], &omega[0], 1.0, -sigma_ * 0.5, rs);

	// в lomega содержитс€ права€ часть, котора€ не мен€етс€ при итераци€х!
	// права€ часть только на границе !

	for (int it = 0; it < 10; ++it) {
		// 0.5(w+w) + l + h <- дл€ вычислени€ якобиана это надо знать и на границе!
		vector_sum1(&omega_lh[0], &omega_1[0], &omega_0[0], 0.5, 0.5, sz);
		vector_sum(&omega_lh[0], &omega_lh[0], &lh_[0], sz);
		vector_sum1(&prev_psi[0], &X_0[0], &psi[0], 0.5, 0.5, sz);
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h)
		j_.calc2(&jac[0], &prev_psi[0], &omega_lh[0]);
		// w/dt + mu \Delta w / 2 - \sigma w/2 -
		// - J(0.5(u+u), 0.5(w+w)) - J(0.5(u+u), l + h) + f(x, y)
#pragma omp parallel for
		for (int i = 0; i < rs; ++i) {
			int point = m_.inner[i];
			double x  = m_.ps[point].x();
			double y  = m_.ps[point].y();

			//omega[i] = omega[i] / tau_ + rp_(x, y, t, mu_, sigma_);
			omega[i] = lomega[i] - jac[i] + rp_(x, y, t, mu_, sigma_);
				//f_[point];
		}

		right_part_cb_data data2;
		//генератор правой части учитывает то, что функци€ задана внутри!!!
		data2.F   = &omega[0];
		data2.bnd = bnd; //TODO: а чему у нас на кра€х равно omega?
		data2.d   = this;

		generate_right_part(&rp[0], m_, 
			(right_part_cb_t)right_part_cb, (void*)&data2);

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

void BarVortex::S_step(double * Ans, const double * F)
{
	calc(Ans, F, 0, 0);
}

void BarVortex::L_step(double * Ans, const double * F)
{
}

void BarVortex::L_1_step(double * Ans, const double * F)
{
}

void BarVortex::LT_step(double * Ans, const double * F)
{
}
