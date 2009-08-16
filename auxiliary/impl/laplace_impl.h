#ifndef LAPL_PRIVATE_H
#define LAPL_PRIVATE_H
/* -*- charset: utf-8 -*- */
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

namespace Laplace_Private 
{
	double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i, int point_j,
		int, int,
		void *);

	double
	laplace_bnd1_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int, int,
		void * d);

	double
	laplace_bnd2_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int, int,
		void * );

double 
laplace_integrate_cb( const Polynom & phi_i,
                      const Polynom & phi_j, 
                      const Triangle & trk, /* ����� ������������ */
                      const Mesh & m,
                      int point_i,
                      int point_j,
                      int, int,
                      void * user_data);
}

template < typename T >
void Laplace < T >::solve(T * Ans, const T * F, const T * bnd)
{
	//���� ���������� ������ �������
	int sz  = (int)m_.ps.size();
	int ntr = (int)m_.tr.size();
	int rs  = (int)m_.inner.size();     //�����������

	Array b(rs);      // ������ �����
	Array x(rs);      // �����

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

template < typename T >
Laplace < T >::Laplace(const Mesh & m): m_(m), 
	idt_((int)m.inner.size()), laplace_((int)m.inner.size()),
//	fidt_(m.ps.size()), flaplace_(m.ps.size()),
	bnd1_((int)m.inner.size()), bnd2_((int)m.inner.size()), 
	bnd3_((int)m.inner.size())
{
	generate_matrix(idt_, m, Laplace_Private::id_cb, (void*)0);
	generate_matrix(laplace_, m, Laplace_Private::laplace_integrate_cb, (void*)0);
//	generate_full_matrix(fidt_, m, id_cb, (void*)0);
//	generate_full_matrix(flaplace_, m, laplace_integrate_cb, (void*)0);
	generate_boundary_matrix(bnd1_, m_, Laplace_Private::laplace_bnd1_cb, (void*)0);
	generate_boundary_matrix(bnd2_, m_, Laplace_Private::laplace_bnd2_cb, (void*)0);
	generate_boundary_matrix(bnd3_, m_, Laplace_Private::laplace_integrate_cb, (void*)0);
}

template < typename T >
void Laplace < T > ::calc2(T * Ans, const T * F)
{
	int rs = (int)m_.inner.size();
	int os = (int)m_.outer.size();
	int sz = (int)m_.ps.size();

#if 1
	Array in(rs);
	Array out(rs);
	Array tmp(os);
	Timer t;
	u2p(&in[0], F, m_);
	fprintf(stderr, "p2u: %lf\n", t.elapsed()); t.restart();

	proj_bnd(&tmp[0], F, m_);
	fprintf(stderr, "proj_bnd: %lf\n", t.elapsed()); t.restart();

	laplace_.mult_vector(&out[0], &in[0]);
	fprintf(stderr, "mult_vector: %lf\n", t.elapsed()); t.restart();

	bnd3_.mult_vector(&in[0], &tmp[0]);
	fprintf(stderr, "mult_vector: %lf\n", t.elapsed()); t.restart();

	vec_sum(&out[0], &out[0], &in[0], (int)in.size());
	fprintf(stderr, "vec_sum: %lf\n", t.elapsed()); t.restart();

	idt_.solve(Ans, &out[0]);
	fprintf(stderr, "solve: %lf\n", t.elapsed()); t.restart();

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
 * �������� ������� �� ������� �� ���������, ������� ���������� ����
 * ������� ������ ���� �� ������, ��� �������!
 *
 * ���� ���� �������� ������� ������ � ����� ����� ���������, ��
 * �������� 0 ������ �������.
 */
template < typename T >
void Laplace < T > ::calc1(T * Ans, const T * F, const T * bnd)
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
	Array out(m_.inner.size());
	calc2(&out[0], F);
	p2u(Ans, &out[0], bnd, m_);
#endif
}

namespace Chafe_Private {

static double f(/*double u,*/ double x, double y, double t, double mu, double sigma)
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
chafe_integrate_cb( const Polynom & phi_i,
                    const Polynom & phi_j, 
                    const Triangle & trk, 
                    const Mesh & m, int point_i, int point_j,
                    int, int,
                    const ChafeConfig * d);

template < typename T >
struct chafe_right_part_cb_data
{
	const T * F;
	const T * bnd;
	const ChafeConfig * d;
};

template < typename T >
double 
chafe_right_part_cb(  const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
					  int i, int j,
                      chafe_right_part_cb_data < T > * d)
{
	const T * F = d->F;
	double b = 0.0;

//	b = F[point_j] * integrate(phi_i * phi_j, trk, m.ps);

	if (m.ps_flags[point_j] == 1) { // �� �������
		int j0         = m.p2io[point_j]; //����� ������� �����
		const T  * bnd = d->bnd;
		b += - (double)bnd[j0] * chafe_integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, i, j, d->d);
	} 
	else {
		b += (double)F[point_j] * integrate(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

};

template < typename T >
Chafe < T > ::Chafe(const Mesh & m, double tau, double sigma, double mu)
	: ChafeConfig(tau, sigma, mu),
	m_(m), laplace_(m), A_((int)m.inner.size())
{
	/* ������� ����� ����� */
	/* ��������(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, Chafe_Private::chafe_integrate_cb, this);
}

/* \frac{du}{dt} = \mu \Delta u - \sigma u + f (u) */
template < typename T >
void Chafe < T > ::solve(T * Ans, const T * X0,
						const T * bnd, double t)
{
	int rs  = (int)m_.inner.size();
	int sz  = (int)m_.ps.size();
	ArrayDevice u(rs);
	ArrayDevice p(sz);
	ArrayHost   hp(sz);
	ArrayDevice delta_u(rs);

	ArrayHost   rp(rs);
	ArrayDevice crp(rs);

	// ���������� ������ �����
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	u2p(&u[0], X0, m_);
	laplace_.calc2(&delta_u[0], X0);

	// u/dt + mu \Delta u / 2
	vec_sum1(&delta_u[0], &u[0], &delta_u[0], (T)(1.0 / tau_), (T)(mu_ * 0.5), rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vec_sum1(&delta_u[0], &delta_u[0], &u[0], (T)1.0, (T)(-sigma_ * 0.5), rs);

	// TODO: ��� ���� ������� ����� ������� ��������
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		int point = m_.inner[i];
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();

		rp[i] = (T)Chafe_Private::f(/*u[i],*/ x, y, t, mu_, sigma_);
	}
	vec_copy_from_host(&crp[0], &rp[0], rs);
	vec_sum(&u[0], &delta_u[0], &crp[0], rs);

	// ������ ����� �� ������� �� ����� !
	p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	vec_copy_from_device(&hp[0], &p[0], sz);

	Chafe_Private::chafe_right_part_cb_data < T > data2;
	data2.F   = &hp[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		Chafe_Private::chafe_right_part_cb < T >, &data2);

	vec_copy_from_host(&crp[0], &rp[0], rs);
	phelm::solve(Ans, bnd, &crp[0], A_, m_);
}

#endif /* LAPL_PRIVATE_H */
