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

static double 
jacobian(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, 
	 const Mesh & m, int i, int j, void * data)
{
	Polynom pt1 = diff(phi_i, 1) * diff(phi_j, 0);
	Polynom pt2 = diff(phi_i, 0) * diff(phi_j, 1);

	Point p = m.ps[i].p[0];
	return (pt1.apply(p.x, p.y) - pt2.apply(p.x, p.y)) / cos(p.x);
}

static double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		void *)
{
	return integrate_cos(phi_i * phi_j, trk, m.ps);
}

/**
 * J(u,v)=1/cos(phi) (du/d\la dv/d\phi - du/d\phi dv/d\la)
 */
Jacobian::Jacobian(const Mesh & m): m_(m), idt_((int)m.inner.size())
{
	generate_matrix(idt_, m, id_cb, 0);
}

void Jacobian::calc1(double * Ans, const double * u, const double * v, const double * bnd)
{
	vector < double > p1(m_.inner.size());
	calc2(&p1[0], u, v);
	mke_p2u(Ans, &p1[0], bnd, m_);
}

static double diff_1_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double r = u[j] * integrate(poly, trk, m.ps);
	return r;
}

static double diff_1_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double r = u[j] * integrate_cos(poly, trk, m.ps);
	return r;
}

static double diff_2_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double r = u[j] * integrate(poly, trk, m.ps);
	return r;
}

static double diff_2_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double r = u[j] * integrate_cos(poly, trk, m.ps);
	return r;
}

void Jacobian::calc2(double * Ans, const double * u, const double * v)
{
	int rs = (int)m_.inner.size();
	int sz = (int)m_.ps.size();
#if 0
	vector < double > rp(sz);
	convolution(&rp[0], u, v, m_, (scalar_cb_t)jacobian, 0);
	mke_u2p(Ans, &rp[0], m_);
#endif

#if 1
	vector < double > rp(rs);
	vector < double > pt1(rs);
	vector < double > pt2(rs);
	double * tmp = Ans;

	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_cos_rp, (void*)u);
	idt_.solve(&pt1[0], &rp[0]);
	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_rp, (void*)v);
	idt_.solve(&tmp[0], &rp[0]);
	vector_mult(&pt1[0], &pt1[0], &tmp[0], (int)pt1.size());

	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_cos_rp, (void*)u);
	idt_.solve(&pt2[0], &rp[0]);
	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_rp, (void*)v);
	idt_.solve(&tmp[0], &rp[0]);
	vector_mult(&pt2[0], &pt2[0], &tmp[0], (int)pt1.size());

	vector_diff(Ans, &pt1[0], &pt2[0], (int)pt1.size());
#endif
}

double laplace(const Polynom & phi_i, const Polynom & phi_j, 
		const Triangle & trk, const Mesh::points_t & ps);

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

static double coriolis(double phi, double lambda)
{
	double omg = 0.0000727000000000;
	double l = omg * 2.0 * sin(phi);
	double h = cos(2.0 * lambda) * ipow(sin(2.0 * phi), 2);;
	return l + h;
}

BarVortex::BarVortex(const Mesh & m, double tau, 
					 double sigma, double mu)
					 : m_(m), l_(m), j_(m), A_(m.inner.size()),
					 tau_(tau), sigma_(sigma), mu_(mu)
{
	lh_.resize(m_.ps.size());
	mke_proj(&lh_[0], m_, coriolis);

	/* Матрица левой части совпадает с Чафе-Инфантом на сфере */
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
	double b;

	b = F[m.p2io[point_j]] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += -bnd[j0] * BarVortex::integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, d->d);		
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
	int rs = (int)m_.inner.size();
	int sz = (int)m_.ps.size();

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
	vector_sum1(&lomega[0], &omega[0], &lomega[0], 1.0 / tau_, mu_ * 0.5, rs);

	// w/dt + mu \Delta w / 2 - \sigma w/2
	vector_sum1(&lomega[0], &lomega[0], &omega[0], 1.0, -sigma_ * 0.5, rs);

	while (true) {
		// 0.5(w+w) + l + h
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

			omega[i] = lomega[i] - jac[i] + f(x, y, t, mu_, sigma_);
		}

		right_part_cb_data data2;
		data2.F   = &omega[0];
		data2.bnd = bnd; //TODO: а чему у нас на краях равно omega?
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
		//	fprintf(stdout, "%le\n", nr);
			if (nr < 1e-8) {
				break;
			}
		}
	}
}
