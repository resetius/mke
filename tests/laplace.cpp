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
#include <math.h>

#include "mke.h"
#include "util.h"
#include "solver.h"
#include "laplace.h"

using namespace std;

struct laplace_right_part_cb_data
{
	const double * F;
	const double * bnd;
};

static double
laplace(const Polynom & phi_i,
		const Polynom & phi_j,
		int trk_i,
		const Mesh & m
       )
{
	const Triangle & trk    = m.tr[trk_i];
	Polynom poly = diff(phi_j, 0) * diff(phi_i, 0)
		+ diff(phi_j, 1) * diff(phi_i, 1);
	return -integrate(poly, trk, m.ps);
}

static double 
laplace_right_part_cb( const Polynom & phi_i,
                       const Polynom & phi_j,
                       int point, /* номер точки */
                       int trk_i, /* номер треугольника */
                       const Mesh & m,
                       laplace_right_part_cb_data * d)
{
	const double * F = d->F;
	const Triangle & trk    = m.tr[trk_i];
	double b;

	if (m.ps_flags[point] == 1) { // на границе
		int j0       = m.p2io[point]; //номер внешней точки
		const double * bnd = d->bnd;
		b = - bnd[j0] * laplace(phi_i, phi_j, trk_i, m);
	} else {
		b = F[point] * integrate(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

static double 
laplace_integrate_cb( const Polynom & phi_i,
                      const Polynom & phi_j, 
                      int point, /* номер точки */
                      int trk_i, /* номер треугольника */
                      const Mesh & m,
                      void * user_data)
{
	return laplace(phi_i, phi_j, trk_i, m);
}

void laplace_solve(double * Ans, const Mesh & m, 
				   const double * F, const double * bnd)
{
	//пока используем первый пор€док
	int sz  = m.ps.size();
	int ntr = m.tr.size();
	int rs  = m.inner.size();     //размерность

	vector < double > b(rs);      // права€ часть
	vector < double > x(rs);      // ответ
	Matrix A(rs);

	Timer full;

	laplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = bnd;

	generate_matrix(A, m, laplace_integrate_cb, 0);
	generate_right_part(&b[0], m, (right_part_cb_t)(laplace_right_part_cb), (void*)&d);
	//A.print();
	//vector_print(&b[0], rs);
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed());
	mke_solve(Ans, bnd, &b[0], A, m);
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed()); 
}

static double f(double u, double x, double y, double t, 
				double mu, double sigma)
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
                     int point, /* номер точки */
                     int trk_i, /* номер треугольника */
                     const Mesh & m,
					 Chafe::integrate_cb_data * d)
{
	double tau   = d->tau;
	double mu    = d->mu;
	double sigma = d->sigma;

	const Triangle & trk  = m.tr[trk_i];
	double pt1, pt2;

	pt1  = integrate(phi_j * phi_i, trk, m.ps);
	pt1 *= 1.0 / tau + sigma * 0.5;

	double a = laplace(phi_i, phi_j, trk_i, m);

	pt2  = a;
	pt2 *= -0.5 * mu;

	return pt1 + pt2;
}

struct chafe_right_part_cb_data
{
	const double * F;
	const double * bnd;
	Chafe::integrate_cb_data * d2;
};

static double 
chafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      int point, /* номер точки */
                      int trk_i, /* номер треугольника */
                      const Mesh & m,
                      chafe_right_part_cb_data * d)
{
	const double * F = d->F;
	const Triangle & trk    = m.tr[trk_i];
	double b;

	if (m.ps_flags[point] == 1) { // на границе
		int j0       = m.p2io[point]; //номер внешней точки
		const double * bnd = d->bnd;
		b = - bnd[j0] * chafe_integrate_cb(phi_i, phi_j, 
			point, trk_i, m, d->d2);
	} else {
		b = F[point] * integrate(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

void laplace_calc(double * Ans, const double * F, const double * bnd, const Mesh & m)
{
	vector < double > p1(m.inner.size());
	vector < double > p2(m.inner.size());
	Matrix laplace_(m.inner.size());
	generate_matrix(laplace_, m, laplace_integrate_cb, 0);

	mke_u2p(&p1[0], &F[0], m);
	laplace_.mult_vector(&p2[0], &p1[0]);
	mke_p2u(Ans, &p2[0], bnd, m);
}

Chafe::Chafe(const Mesh & m, double tau, double sigma, double mu)
	: m_(m), laplace_(m.inner.size()), A_(m.inner.size()), 
	tau_(tau), mu_(mu), sigma_(sigma)
{
	/* Ћапласиан */
	generate_matrix(laplace_, m_, laplace_integrate_cb, 0);

	/* ћатрица левой части */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/

	data1_.tau   = tau;
	data1_.sigma = sigma;
	data1_.mu    = mu;
	generate_matrix(A_, m_, (integrate_cb_t)chafe_integrate_cb, (void*)&data1_);
}

/**
 * \f$\frac{du}{dt} = \mu \delta u - \sigma u + f (u)\f$
 */
void Chafe::solve(double * Ans, const double * X0,
						const double * bnd, double t)
{
	int rs  = m_.inner.size();
	int sz  = m_.ps.size();
	vector < double > u(rs);
	vector < double > p(sz);
	vector < double > delta_u(rs);
	vector < double > rp(rs);

	// генерируем правую часть
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	mke_u2p(&u[0], X0, m_);
	laplace_.mult_vector(&delta_u[0], &u[0]);

	// u/dt + mu \Delta u / 2
	vector_sum1(&delta_u[0], &u[0], &delta_u[0], 1.0 / tau_, mu_ * 0.5, rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vector_sum1(&delta_u[0], &delta_u[0], &u[0], 1.0, -sigma_ * 0.5, rs);

	// TODO: тут надо сделать метод простой итерации
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		int point = m_.inner[i];
		double x  = m_.ps[point].x;
		double y  = m_.ps[point].y;

		u[i] = delta_u[i] + f(u[i], x, y, t, mu_, sigma_);
	}

	mke_p2u(&p[0], &u[0], bnd, m_);
	chafe_right_part_cb_data data2;
	data2.F   = &p[0];
	data2.bnd = bnd;
	data2.d2  = &data1_;
	generate_right_part(&rp[0], m_, 
		(right_part_cb_t)chafe_right_part_cb, (void*)&data2);

	mke_solve(Ans, bnd, &rp[0], A_, m_);
}
