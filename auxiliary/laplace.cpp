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

static double
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

static double
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

static double 
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

void Laplace::solve(double * Ans, const double * F, const double * bnd)
{
	//пока используем первый порядок
	int sz  = (int)m_.ps.size();
	int ntr = (int)m_.tr.size();
	int rs  = (int)m_.inner.size();     //размерность

	vec b(rs);      // правая часть
	vec x(rs);      // ответ

	Timer full;

#if 0
	laplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = bnd;
	generate_right_part(&b[0], m_, (right_part_cb_t)(laplace_right_part_cb), (void*)&d);
#endif

#if 1
	u2p(&x[0], F, m_);
	idt_.mult_vector(&b[0], &x[0]);
	if (bnd) {
		bnd2_.mult_vector(&x[0], bnd);
		vec_sum(&b[0], &b[0], &x[0], (int)x.size());
	}

//	vector < double > tmp(m_.outer.size());
//	proj_bnd(&tmp[0], F, m_);
//	bnd1_.mult_vector(&x[0], &tmp[0]);
//	vec_sum(&b[0], &b[0], &x[0], x.size());
#endif	

	phelm::solve(Ans, bnd, &b[0], laplace_, m_);
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

static double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i, int point_j,
		int, int,
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

Laplace::Laplace(const Mesh & m): m_(m), 
	idt_((int)m.inner.size()), laplace_((int)m.inner.size()),
//	fidt_(m.ps.size()), flaplace_(m.ps.size()),
	bnd1_((int)m.inner.size()), bnd2_((int)m.inner.size()), 
	bnd3_((int)m.inner.size())
{
	generate_matrix(idt_, m, id_cb, (void*)0);
	generate_matrix(laplace_, m, laplace_integrate_cb, (void*)0);
//	generate_full_matrix(fidt_, m, id_cb, (void*)0);
//	generate_full_matrix(flaplace_, m, laplace_integrate_cb, (void*)0);
	generate_boundary_matrix(bnd1_, m_, laplace_bnd1_cb, (void*)0);
	generate_boundary_matrix(bnd2_, m_, laplace_bnd2_cb, (void*)0);
	generate_boundary_matrix(bnd3_, m_, laplace_integrate_cb, (void*)0);
}

void Laplace::calc2(double * Ans, const double * F)
{
	int rs = (int)m_.inner.size();
	int os = (int)m_.outer.size();
	int sz = (int)m_.ps.size();

#if 1
	vec in(rs);
	vec out(rs);
	vec tmp(os);
	u2p(&in[0], F, m_);
	proj_bnd(&tmp[0], F, m_);
	laplace_.mult_vector(&out[0], &in[0]);
	bnd3_.mult_vector(&in[0], &tmp[0]);
	vec_sum(&out[0], &out[0], &in[0], (int)in.size());
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
	u2p(Ans, &out[0], m_);
#endif
}

/*
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

	p2u(Ans, &p1[0], bnd, m_);
#endif
#if 1
	vec out(m_.inner.size());
	calc2(&out[0], F);
	p2u(Ans, &out[0], bnd, m_);
#endif
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

