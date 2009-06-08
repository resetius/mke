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

#include "baroclin.h"
#include "util.h"

using namespace std;

static elements_t
integrate_cb( const Polynom & phi_i,
              const Polynom & phi_j, 
              const Triangle & trk,
              const Mesh & m,
              int point_i, int point_j,
              int i, int j,
              Baroclin * d);

Baroclin::Baroclin(const Mesh & m, rp_t f, rp_t g,
				   coriolis_t coriolis,
         double tau, double sigma, double mu, 
		 double sigma1, double mu1, double alpha)
	: SphereNorm(m), m_(m), l_(m), j_(m), 
	A_(4 * (int)m.inner.size()), 
	tau_(tau), sigma_(sigma), mu_(mu),
	sigma1_(sigma1), mu1_(mu1), 
	alpha_(alpha), f_(f), g_(g), coriolis_(coriolis)
{
	theta_ = 0.5;
	lh_.resize(m_.ps.size());
	mke_proj(&lh_[0], m_, coriolis);

	/* Матрица левой части совпадает с Чафе-Инфантом на сфере */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, integrate_cb, this);
}

struct right_part_cb_data
{
	const double * F;
	const double * G;
	const double * BW1;
	const double * BW2;
	const double * BU1;
	const double * BU2;
	Baroclin  * d;
};

elements_t
integrate_cb( const Polynom & phi_i,
              const Polynom & phi_j, 
              const Triangle & trk,
              const Mesh & m,
              int point_i, int point_j,
              int i, int j,
              Baroclin * d)
{
	elements_t r;
	int rs = (int)m.inner.size();
	double tau    = d->tau_;
	double mu     = d->mu_;
	double mu1    = d->mu1_;
	double sigma  = d->sigma_;
	double sigma1 = d->sigma1_;
	double alpha  = d->alpha_;
	double theta  = d->theta_;

	/**
	 * 4nx4n matrix:
	 * w1 - [0,   n)
	 * w2 - [n,  2n)
	 * u1 - [2n, 3n)
	 * u3 - [3n, 4n)
	 *
	 * 1: w1 / dt + 0.5 theta * sigma (w1 - w2) - theta mu L(w1) = F1
	 * 2: w2 / dt - alpha^2 u2/dt + 0.5 theta * sigma (w1 + w2) -
	 * - theta mu L(w2) - alpha^2 (- theta mu1 w2  + theta sigma1 u2) = F2
	 * 3: w1 - L(u1) = 0
	 * 4: w2 - L(u2) = 0
	 */

	double a = integrate_cos(phi_i * phi_j, trk, m.ps);
	double b = slaplace(phi_j, phi_i, trk, m.ps);
	r.reserve(8);
	// 1:
	// w1 / dt + 0.5 theta * sigma w1 - theta mu L(w1)
	r.push_back(Element(i, j, a * (1.0 / tau + 0.5 * theta * sigma) - theta * mu * b));
	// - 0.5 theta sigma w2
	r.push_back(Element(i, j + rs, -0.5 * a * sigma * theta));
	// 2:
	// 0.5 theta * sigma w1
	r.push_back(Element(i + rs, j, 0.5 * theta * sigma * a));
	// w2 / dt + 0.5 theta * sigma w2 - theta mu L(w2) + alpha^2 theta mu1 w2
	r.push_back(Element(i + rs, j + rs, a * (1.0 / tau + 0.5 * theta * sigma + alpha * alpha * theta * mu1) 
				- theta * mu * b));
	// - alpha^2 u2/dt - alpha^2 theta sigma1 u2
	r.push_back(Element(i + rs, j + 3*rs, a * (-alpha * alpha / tau - alpha * alpha * theta * sigma1)));
	// 3:
	// w1
	r.push_back(Element(i + 2*rs, j, a));
	// -L(u1)
	r.push_back(Element(i + 2*rs, j + 2*rs, -b));
	// 4:
	// w2
	r.push_back(Element(i + 3*rs, j + rs, a));
	// - L(u2)
	r.push_back(Element(i + 3*rs, j + 3*rs, -b));

	return r;
}

elements_t
right_part_cb( const Polynom & phi_i,
               const Polynom & phi_j,
               const Triangle & trk,
               const Mesh & m,
               int point_i, int point_j,
               int i, int j,
               right_part_cb_data * d)
{
	elements_t r;
	int rs = (int)m.inner.size();
	const double * F = d->F;
	const double * G = d->G;

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0        = m.p2io[point_j]; //номер внешней точки
		elements_t r1 = integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, i, j, d->d);
		double rp[] = {0, 0, 0, 0};
		for (elements_t::iterator it = r1.begin(); it != r1.end(); ++it)
		{
			Element & e = *it;
			// TODO: fixme!
			double * b1 = 0;

			// detect equation number
			if (e.i < rs) {
				b1 = &rp[0];
			} else if (e.i < 2 * rs) {
				b1 = &rp[1];
			} else if (e.i < 3 * rs) {
				b1 = &rp[2];
			} else { // if (e.i < 4 * rs)
				b1 = &rp[3];
			}

			double & b = *b1;
			if (e.j < rs) {
				// w1
				if (d->BW1) {
					b += -d->BW1[j0] * e.a;
				}
			} else if (e.j < 2 * rs) {
				// w2
				if (d->BW2) {
					b += -d->BW2[j0] * e.a;
				}
			} else if (e.j < 3 * rs) {
				// u1
				if (d->BU1) {
					b += -d->BU1[j0] * e.a;
				}
			} else { //e.j < 4 * rs
				// u2
				if (d->BU2) {
					b += -d->BU2[j0] * e.a;
				}
			}
		}

		r.push_back(Element(i,          j, rp[0]));
		r.push_back(Element(i + rs,     j, rp[1]));
		r.push_back(Element(i + 2 * rs, j, rp[2]));
		r.push_back(Element(i + 3 * rs, j, rp[3]));
	} else {
		double a = integrate_cos(phi_i * phi_j, trk, m.ps);
		// F
		r.push_back(Element(i, j, F[m.p2io[point_j]] * a));
		// G
		r.push_back(Element(i + rs, j, G[m.p2io[point_j]] * a));
	}

	return r;
}

/**
 * (u1, u2) -> (u11, u21)
 * d L(u1)/dt + J(u1, L(u1) + l + h ?) + J(u2, L(u2)) + sigma/2 L(u1 - u2) - mu LL(u1) = f(phi, lambda)
 * d L(u2)/dt + J(u1, L(u2)) + J(u2, L(u1) + l + h?) + sigma/2 L(u1 + u2) - mu LL(u2) -
 *   - alpha^2 (d u2/dt + J(u1, u2) - mu1 L(u2) + sigma1 u2 + g(phi, lambda))= 0
 * L = Laplace
 */
void Baroclin::calc(double * u11,  double * u21, 
		const double * u1, const double * u2, 
		const double * bnd, double t)
{
	int rs = (int)m_.inner.size(); // размерность внутренней области
	int sz = (int)m_.ps.size();    // размерность полная

	// правая часть 1:
	// -J(0.5(u1+u1), 0.5(w1+w1)+l+h) - J(0.5(u2+u2),w2+w2)+
	// + w1/tau - 0.5 (1-theta)sigma (w1-w2)+mu(1-theta)(\Delta w1)
	// правая часть 2:
	// -J(0.5(u1+u1), 0.5(w2+w2)) - J(0.5(u2+u2), 0.5(w1+w1)+l+h) -
	// - 0.5 (1-theta)sigma (w1 + w2) + (1-theta) mu \Delta w2
	// + w2/tau - alpha^2 u2/tau + alpha^2 J(0.5(u1+u1), 0.5(u2+u2)) -
	// - alpha^2 (1-theta) mu1 w2 +
	// + alpha^2 sigma1 (1-theta) u2 + alpha^2 f(phi, lambda)

	vector < double > w1(sz);
	vector < double > w2(sz);
	vector < double > dw1(sz);
	vector < double > dw2(sz);
	vector < double > FC(sz);
	vector < double > GC(sz);

	// next
	vector < double > u1_n(sz);
	vector < double > u2_n(sz);
	vector < double > w1_n(sz);
	vector < double > w2_n(sz);

	// tmp
	vector < double > tmp1(sz);
	vector < double > tmp2(sz);

	// jac inner!
	vector < double > jac1(rs);
	vector < double > jac2(rs);
	vector < double > jac3(rs);

	//
	vector < double > F(rs);
	vector < double > G(rs);

	vector < double > rp(4*rs);
	vector < double > ans(4*rs);

	l_.calc1(&w1[0], &u1[0], bnd);
	l_.calc1(&w2[0], &u2[0], bnd);

	l_.calc1(&dw1[0], &w1[0], bnd);
	l_.calc1(&dw2[0], &w2[0], bnd);

	// w1/tau - 0.5 (1-theta)sigma(w1-w2) + mu(1-theta)(\Delta w1)
	vec_sum1(&FC[0], &w1[0], &w2[0], 
		-0.5 * (1.0 - theta_) * sigma_, 0.5 * (1.0 - theta_) * sigma_, sz);
	vec_sum1(&FC[0], &FC[0], &dw1[0], 1.0, mu_ * (1.0 - theta_), sz);
	vec_sum1(&FC[0], &FC[0], &w1[0], 1.0, 1.0 / tau_, sz);

	// w2/tau - 0.5 (1-theta)sigma (w1 + w2) + (1-theta) mu \Delta w2 -
	// - alpha^2 u2/tau - alpha^2 (1-theta) mu1 w2 + alpha^2 sigma1 (1-theta) u2
	vec_sum1(&GC[0], &w1[0], &w2[0],
			-0.5 * (1.0 - theta_) * sigma_, -0.5 * (1.0 - theta_) * sigma_, sz);
	vec_sum1(&GC[0], &GC[0], &dw2[0], 1.0, mu_ * (1.0 - theta_), sz);
	vec_sum1(&GC[0], &GC[0], &w2[0], 1.0, 1.0 / tau_, sz);
	vec_sum1(&GC[0], &GC[0], &u2[0], 1.0, -alpha_ * alpha_ / tau_, sz);
	vec_sum1(&GC[0], &GC[0], &w2[0], 1.0, -alpha_ * alpha_ * mu1_ * (1-theta_), sz);
	vec_sum1(&GC[0], &GC[0], &u2[0], 1.0, alpha_ * alpha_ * sigma1_ * (1-theta_), sz);

	memcpy(&u1_n[0], &u1[0], sz * sizeof(double));
	memcpy(&u2_n[0], &u2[0], sz * sizeof(double));
	memcpy(&w1_n[0], &w1[0], sz * sizeof(double));
	memcpy(&w2_n[0], &w2[0], sz * sizeof(double));

	for (int i = 0; i <sz; ++i) {
		int point = i;
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();
		if (f_) {
			FC[i] += f_(x, y, t, sigma_, mu_, sigma1_, 
				mu1_, alpha_, theta_);
		}
		if (g_) {
			GC[i] += alpha_ * alpha_ 
				* g_(x, y, t, sigma_, mu_, sigma1_, mu1_, alpha_, theta_);
		}
	}

	right_part_cb_data data2;
	data2.F   = &F[0];
	data2.G   = &G[0];
	data2.BW1 = bnd;
	data2.BW2 = bnd;
	data2.BU1 = bnd;
	data2.BU2 = bnd;
	data2.d   = this;

	for (int it = 0; it < 20; ++it) {
		// - J(0.5(u1+u1), 0.5(w1+w1)+l+h) - J(0.5(u2+u2),w2+w2)
		// J(0.5(u1+u1), 0.5(w1+w1)+l+h)
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta_, theta_, sz);
		vec_sum1(&tmp2[0], &w1[0], &w1_n[0], 1.0 - theta_, theta_, sz);
		vec_sum(&tmp2[0], &tmp2[0], &lh_[0], sz);
		j_.calc2(&jac1[0], &tmp1[0], &tmp2[0]);
		// J(0.5(u2+u2),w2+w2)
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], 1.0 - theta_, theta_, sz);
		vec_sum1(&tmp2[0], &w2[0], &w2_n[0], 1.0 - theta_, theta_, sz);
		j_.calc2(&jac2[0], &tmp1[0], &tmp2[0]);

		vec_sum1(&F[0], &jac1[0], &jac2[0], -1.0, -1.0, rs);

		// -J(0.5(u1+u1), 0.5(w2+w2)) - J(0.5(u2+u2), 0.5(w1+w1)+l+h) +
		// + alpha^2 J(0.5(u1+u1), 0.5(u2+u2))
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta_, theta_, sz);
		vec_sum1(&tmp2[0], &w2[0], &w2_n[0], 1.0 - theta_, theta_, sz);
		j_.calc2(&jac1[0], &tmp1[0], &tmp2[0]);
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], 1.0 - theta_, theta_, sz);
		vec_sum1(&tmp2[0], &w1[0], &w1_n[0], 1.0 - theta_, theta_, sz);
		vec_sum(&tmp2[0], &tmp2[0], &lh_[0], sz);
		j_.calc2(&jac2[0], &tmp1[0], &tmp2[0]);
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta_, theta_, sz);
		vec_sum1(&tmp2[0], &u2[0], &u2_n[0], 1.0 - theta_, theta_, sz);
		j_.calc2(&jac3[0], &tmp1[0], &tmp2[0]);
		vec_sum1(&G[0], &jac1[0], &jac2[0], -1.0, -1.0, rs);
		vec_sum1(&G[0], &G[0], &jac3[0], 1.0, alpha_ * alpha_, rs);

		for (int i = 0; i < rs; ++i) {
			int point = m_.inner[i];
			F[i] += FC[point];
			G[i] += GC[point];
		}

		generate_right_part(&rp[0], m_, right_part_cb, &data2);
		A_.solve(&ans[0], &rp[0]);
		mke_p2u(&w1_n[0], &ans[0],    bnd, m_);
		mke_p2u(&w2_n[0], &ans[rs],   bnd, m_);
		mke_p2u(&u1_n[0], &ans[2*rs], bnd, m_);
		mke_p2u(&u2_n[0], &ans[3*rs], bnd, m_);

		double nr1 = dist(&u1[0], &u1_n[0]);
		double nr2 = dist(&u2[0], &u2_n[0]);
		double nr  = std::max(nr1, nr2);

		if (nr < 1e-8) {
			break;
		}
	}

	memcpy(u11, &u1_n[0], sz * sizeof(double));
	memcpy(u21, &u2_n[0], sz * sizeof(double));
}

