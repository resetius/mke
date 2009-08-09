#ifndef LAPL_PRIVATE_H
#define LAPL_PRIVATE_H
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
                      const Triangle & trk, /* номер треугольника */
                      const Mesh & m,
                      int point_i,
                      int point_j,
                      int, int,
                      void * user_data);
}

template < typename T >
void Laplace < T >::solve(T * Ans, const T * F, const T * bnd)
{
	//пока используем первый порядок
	int sz  = (int)m_.ps.size();
	int ntr = (int)m_.tr.size();
	int rs  = (int)m_.inner.size();     //размерность

	Array b(rs);      // правая часть
	Array x(rs);      // ответ

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
 * Оператор Лапласа на границе не определен, поэтому вставляйте сюда
 * границу только если вы знаете, что делаете!
 *
 * Если этот оператор Лапласа входит в праву часть уравнения, то
 * напишите 0 вместо границы.
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

#endif /* LAPL_PRIVATE_H */
