/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (������� ���������)
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

/**
 * phi    = x - latitude
 * lambda = y - longitude
 */

/**
 * laplace = laplace1 + laplace2
 */

/**
 * laplace1 = (1 / cos phi d / d phi phi_i, cos phi d / phi phi_j)
 */
static double laplace1(const Polynom & phi_i, const Polynom & phi_j, 
		const Triangle & trk, const Mesh::points_t & ps)
{
	return integrate_cos(diff(phi_i, 0) * diff(phi_j, 0), trk, ps);
}

/**
 * laplace2 = (1 / cos phi d / d lambda phi_i, 1 / cos phi d / d lambda phi_j)
 */
static double laplace2(const Polynom & phi_i, const Polynom & phi_j,
		const Triangle & trk, const Mesh::points_t & ps)
{
	return integrate_1_cos(diff(phi_i, 1) * diff(phi_j, 1), trk, ps);
}

double laplace(const Polynom & phi_i, const Polynom & phi_j, 
		const Triangle & trk, const Mesh::points_t & ps)
{
	return -(laplace1(phi_i, phi_j, trk, ps) + laplace2(phi_i, phi_j, trk, ps));
}

struct slaplace_right_part_cb_data
{
	const double * F;
	const double * bnd;
};

static double 
slaplace_right_part_cb( const Polynom & phi_i,
                        const Polynom & phi_j,
                        const Triangle & trk,
                        const Mesh & m,
			int point_i,
			int point_j,
                        slaplace_right_part_cb_data * d)
{
	const double * F = d->F;
	double b;

	b = F[point_j] * integrate_cos(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1 && d->bnd) { // �� �������
		int j0       = m.p2io[point_j]; //����� ������� �����
		const double * bnd = d->bnd;
		b += -bnd[j0] * laplace(phi_j, phi_i, trk, m.ps);
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
	return integrate_cos(phi_i * phi_j, trk, m.ps);
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
	return -laplace(phi_i, phi_j, trk, m.ps);
}

static double 
slaplace_integrate_cb( const Polynom & phi_i,
                       const Polynom & phi_j, 
                       const Triangle & trk,
                       const Mesh & m,
		       int point_i,
		       int point_j,
                       void * user_data)
{
	double a = laplace(phi_j, phi_i, trk, m.ps);
	return a;
}

void SphereLaplace::solve(double * Ans,
			  const double * F, const double * bnd)
{
	//���� ���������� ������ �������
	int sz  = (int)m_.ps.size();
	int ntr = (int)m_.tr.size();
	int rs  = (int)m_.inner.size();     //�����������

	vector < double > b(rs);      // ������ �����
	vector < double > x(rs);      // �����

	Timer full;
#if 0
	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = bnd;
	generate_right_part(&b[0], m_, (right_part_cb_t)(slaplace_right_part_cb), (void*)&d);
#endif

#if 1
	mke_u2p(&x[0], F, m_);
	idt_.mult_vector(&b[0], &x[0]);
	if (bnd) {
		bnd2_.mult_vector(&x[0], bnd);
		vector_sum(&b[0], &b[0], &x[0], x.size());
	}
//	vector < double > tmp(m_.outer.size()); // not necessary !
//	mke_proj_bnd(&tmp[0], F, m_);           // not necessary !
//	bnd1_.mult_vector(&x[0], &tmp[0]);      // not necessary !
//	vector_sum(&b[0], &b[0], &x[0], x.size()); // not necessary !
#endif

	mke_solve(Ans, bnd, &b[0], laplace_, m_);
#ifdef _DEBUG
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed()); 
#endif
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

static double lp_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		slaplace_right_part_cb_data * d)
{
	const double * F = d->F;
	double b = F[point_j] * laplace(phi_j, phi_i, trk, m.ps);;
#if 0
	if (m.ps_flags[point_j] == 1 && d->bnd) { // �� �������
		int j0       = m.p2io[point_j]; //����� ������� �����
		b += - d->bnd[j0] * id_cb(phi_i, phi_j, 
				trk, m, point_i, point_j, 0);
	}
#endif
	return b;
}

SphereLaplace::SphereLaplace(const Mesh & m): m_(m), 
	idt_((int)m.inner.size()),
	laplace_((int)m.inner.size()), bnd1_(m.inner.size()), bnd2_(m.inner.size()),
	bnd3_(m.inner.size())
{
	generate_matrix(idt_, m, id_cb, 0);
	generate_matrix(laplace_, m, slaplace_integrate_cb, 0);
	generate_boundary_matrix(bnd1_, m_, laplace_bnd1_cb, 0);
	generate_boundary_matrix(bnd2_, m_, laplace_bnd2_cb, 0);
	generate_boundary_matrix(bnd3_, m_, slaplace_integrate_cb, 0);
}

void SphereLaplace::calc2(double * Ans, const double * F)
{
#if 0
	vector < double > rp(m_.inner.size());

	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = 0;
	generate_right_part(&rp[0], m_, (right_part_cb_t)lp_rp, &d);
	idt_.solve(Ans, &rp[0]);
#endif
	int rs = m_.inner.size();
        int os = m_.outer.size();
	vector < double > in(rs);
	vector < double > out(rs);
	vector < double > tmp(os);
	mke_u2p(&in[0], F, m_);
	mke_proj_bnd(&tmp[0], F, m_);
	laplace_.mult_vector(&out[0], &in[0]);
	bnd3_.mult_vector(&in[0], &tmp[0]);
	mke_vector_sum(&out[0], &out[0], &in[0], in.size());
	idt_.solve(Ans, &out[0]);
}

void SphereLaplace::calc1(double * Ans, const double * F, const double * bnd)
{
	vector < double > p1(m_.inner.size());

	calc2(&p1[0], F);
#if 0
	vector < double > rp(m_.inner.size());

	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = 0;
	generate_right_part(&rp[0], m_, (right_part_cb_t)lp_rp, &d);
	idt_.solve(&p1[0], &rp[0]);
#endif
	mke_p2u(Ans, &p1[0], bnd, m_);
}


#if 0
static double f(double u, double x, double y, double t, 
				double mu, double sigma)
{
	double a = exp(t) * sin(y) * sin(2.0 * x);
	double b = -6.0 * exp(t) * sin(y) * sin(2.0 * x);

	return a - mu * b + sigma * a;

//	return (1.0 + 6.0 * mu + sigma) * u;
//	return -u * u * u;
}
#endif

static double f(double u, double x, double y, double t, 
				double mu, double sigma)
{
	return ipow(cos(x),2)*
		(x*cos(y+t)*ipow(cos(x),2)
		+sigma*sin(y+t)*ipow(cos(x),2)*x+9*mu*sin(y+t)*sin(x)*cos(x)
		-15*mu*sin(y+t)*x+20*mu*sin(y+t)*x*ipow(cos(x),2));
}

double 
SphereChafe::schafe_integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
                     SphereChafe * d)
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

struct SphereChafe::schafe_right_part_cb_data
{
	const double * F;
	const double * bnd;
	SphereChafe  * d;
};

double 
SphereChafe::schafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
                      schafe_right_part_cb_data * d)
{
	const double * F = d->F;
	double b;

	if (m.ps_flags[point_j] == 1) { // �� �������
		int j0       = m.p2io[point_j]; //����� ������� �����
		const double * bnd = d->bnd;
		b = -bnd[j0] * SphereChafe::schafe_integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, d->d);
		//b = 0.0;
	} else {
		b = F[point_j] * integrate_cos(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

SphereChafe::SphereChafe(const Mesh & m, double tau, double sigma, double mu)
	: m_(m), laplace_(m), A_((int)m.inner.size()), 
	tau_(tau), mu_(mu), sigma_(sigma)
{
	/* ������� ����� ����� */
	/* ��������(u) = u/dt-mu \Delta u/2 + sigma u/2*/

	generate_matrix(A_, m_, (integrate_cb_t)schafe_integrate_cb, this);
}

/**
 * \f$\frac{du}{dt} = \mu \delta u - \sigma u + f (u)\f$
 */
void SphereChafe::solve(double * Ans, const double * X0,
						const double * bnd, double t)
{
	int rs  = (int)m_.inner.size();
	int sz  = (int)m_.ps.size();
	vector < double > u(rs);
	vector < double > p(sz);
	vector < double > delta_u(rs);
	vector < double > rp(rs);

	// ���������� ������ �����
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	mke_u2p(&u[0], X0, m_);
	laplace_.calc2(&delta_u[0], X0);

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

		u[i] = delta_u[i] + f(u[i], x, y, t, mu_, sigma_);
	}

	// ������ ����� �� ������� �� ����� !!!
	mke_p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	schafe_right_part_cb_data data2;
	data2.F   = &p[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		(right_part_cb_t)schafe_right_part_cb, (void*)&data2);

//	fprintf(stderr, "rp: \n");vector_print(&delta_u[0], rs);
//	fprintf(stderr, "matrix:\n");A_.print();
//	laplace_.print();
	mke_solve(Ans, bnd, &rp[0], A_, m_);
}

