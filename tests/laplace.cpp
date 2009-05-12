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
laplace(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, const Mesh & m)
{
	Polynom poly = diff(phi_j, 0) * diff(phi_i, 0)
		+ diff(phi_j, 1) * diff(phi_i, 1);
	return -integrate(poly, trk, m.ps);
}

static double 
laplace_right_part_cb( const Polynom & phi_i,
                       const Polynom & phi_j,
                       const Triangle & trk, /* номер треугольника */
                       const Mesh & m,
                       int point_i,
		       int point_j,
                       laplace_right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

	//b = F[point_j] * integrate(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) {         // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += - bnd[j0] * laplace(phi_i, phi_j, trk, m);
	}
	else {
		b += F[point_j] * integrate(phi_i * phi_j, trk, m.ps);
	}

	return b;
}

static double
laplace_bnd1_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		void * d)
{
	return integrate(phi_i * phi_j, trk, m.ps);
}

static double
laplace_bnd2_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		void * )
{
	return -laplace(phi_i, phi_j, trk, m);
}

static double 
laplace_integrate_cb( const Polynom & phi_i,
                      const Polynom & phi_j, 
                      const Triangle & trk, /* номер треугольника */
                      const Mesh & m,
                      int point_i,
		      int point_j,
                      void * user_data)
{
	return laplace(phi_i, phi_j, trk, m);
}

void Laplace::solve(double * Ans, const double * F, const double * bnd)
{
	//пока используем первый порядок
	int sz  = m_.ps.size();
	int ntr = m_.tr.size();
	int rs  = m_.inner.size();     //размерность

	vector < double > b(rs);      // правая часть
	vector < double > x(rs);      // ответ

	Timer full;

#if 0
	laplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = bnd;
	generate_right_part(&b[0], m_, (right_part_cb_t)(laplace_right_part_cb), (void*)&d);
#endif

#if 1
	mke_u2p(&x[0], F, m_);
	idt_.mult_vector(&b[0], &x[0]);
	bnd2_.mult_vector(&x[0], bnd);
	vec_sum(&b[0], &b[0], &x[0], x.size());
//	vector < double > tmp(m_.outer.size());
//	mke_proj_bnd(&tmp[0], F, m_);
//	bnd1_.mult_vector(&x[0], &tmp[0]);
	vec_sum(&b[0], &b[0], &x[0], x.size());
#endif	

	mke_solve(Ans, bnd, &b[0], laplace_, m_);
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed()); 
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

double 
Chafe::chafe_integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk, 
                     const Mesh & m, int point_i, int point_j,
		     const Chafe * d)
{
	double tau   = d->tau_;
	double mu    = d->mu_;
	double sigma = d->sigma_;

	double pt1, pt2;

	pt1  = integrate(phi_j * phi_i, trk, m.ps);
	pt1 *= 1.0 / tau + sigma * 0.5;

	pt2  = laplace(phi_i, phi_j, trk, m);
	pt2 *= -0.5 * mu;

	return pt1 + pt2;
}

struct Chafe::chafe_right_part_cb_data
{
	const double * F;
	const double * bnd;
	const Chafe * d;
};

double 
Chafe::chafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
                      chafe_right_part_cb_data * d)
{
	const double * F = d->F;
	double b = 0.0;

//	b = F[point_j] * integrate(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * bnd = d->bnd;
		b += - bnd[j0] * Chafe::chafe_integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, d->d);
	} 
	else {
		b += F[point_j] * integrate(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

static double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i, int point_j,
		void *)
{
	return integrate(phi_i * phi_j, trk, m.ps);
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

	b = F[point_j] * laplace(phi_i, phi_j, trk, m);

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

Laplace::Laplace(const Mesh & m): m_(m), 
	idt_(m.inner.size()), laplace_(m.inner.size()),
//	fidt_(m.ps.size()), flaplace_(m.ps.size()),
	bnd1_(m.inner.size()), bnd2_(m.inner.size()), 
	bnd3_(m.inner.size())
{
	generate_matrix(idt_, m, id_cb, 0);
	generate_matrix(laplace_, m, laplace_integrate_cb, 0);
//	generate_full_matrix(fidt_, m, id_cb, 0);
//	generate_full_matrix(flaplace_, m, laplace_integrate_cb, 0);
	generate_boundary_matrix(bnd1_, m_, laplace_bnd1_cb, 0);
	generate_boundary_matrix(bnd2_, m_, laplace_bnd2_cb, 0);
	generate_boundary_matrix(bnd3_, m_, laplace_integrate_cb, 0);
}

void Laplace::calc2(double * Ans, const double * F)
{
	int rs = m_.inner.size();
	int os = m_.outer.size();
	int sz = m_.ps.size();

#if 1
	vector < double > in(rs);
	vector < double > out(rs);
	vector < double > tmp(os);
	mke_u2p(&in[0], F, m_);
	mke_proj_bnd(&tmp[0], F, m_);
	laplace_.mult_vector(&out[0], &in[0]);
	bnd3_.mult_vector(&in[0], &tmp[0]);
	vec_sum(&out[0], &out[0], &in[0], in.size());
	idt_.solve(Ans, &out[0]);
#endif

#if 0
	vector < double > in(sz);
	vector < double > out(sz);
	vector < double > tmp(sz);

	//flaplace_.mult_vector(&in[0], &F[0]);
	laplace_right_part_cb_data data;
	data.F = &F[0];
	generate_full_right_part(&in[0], m_, (right_part_cb_t)lp_rp, &data);
	fidt_.solve(&out[0], &in[0]);
	mke_u2p(Ans, &out[0], m_);
#endif
}

/**
 * Оператор Лапласа на границе не определен, поэтому вставляйте сюда
 * границу только если вы знаете, что делаете!
 *
 * Если этот оператор Лапласа входит в праву часть уравнения, то
 * напишите 0 вместо границы.
 */
void Laplace::calc1(double * Ans, const double * F, const double * bnd)
{
#if 0
	vector < double > p1(m_.inner.size());

	//calc2(&p1[0], F);

	vector < double > rp(m_.inner.size());

	laplace_right_part_cb_data d;
	d.F = F;
	d.bnd = 0;//bnd;
	generate_right_part(&rp[0], m_, (right_part_cb_t)lp_rp, &d);
	idt_.solve(&p1[0], &rp[0]);

	mke_p2u(Ans, &p1[0], bnd, m_);
#endif
#if 1
	vector < double > out(m_.inner.size());
	calc2(&out[0], F);
	mke_p2u(Ans, &out[0], bnd, m_);
#endif
}

Chafe::Chafe(const Mesh & m, double tau, double sigma, double mu)
	: m_(m), laplace_(m), A_(m.inner.size()), 
	tau_(tau), mu_(mu), sigma_(sigma)
{
	/* Матрица левой части */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, (integrate_cb_t)chafe_integrate_cb, this);
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
	mke_p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	chafe_right_part_cb_data data2;
	data2.F   = &p[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		(right_part_cb_t)chafe_right_part_cb, (void*)&data2);

	mke_solve(Ans, bnd, &rp[0], A_, m_);
}

