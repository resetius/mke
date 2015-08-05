#ifndef LAPL_PRIVATE_H
#define LAPL_PRIVATE_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2011 Alexey Ozeritsky (Алексей Озерицкий)
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

double id_cb (const Polynom & phi_i,
              const Polynom & phi_j,
              const Triangle & trk,
              int z,
              const Mesh & m,
              int point_i, int point_j,
              int, int,
              void *);

double
laplace_bnd2_cb ( const Polynom & phi_i,
                  const Polynom & phi_j,
                  const Triangle & trk,
                  int z,
                  const Mesh & m,
                  int point_i,
                  int point_j,
                  int, int,
                  void * );

double
laplace_integrate_cb ( const Polynom & phi_i,
                       const Polynom & phi_j,
                       const Triangle & trk, /* номер треугольника */
                       int z,
                       const Mesh & m,
                       int point_i,
                       int point_j,
                       int, int,
                       void * user_data);


template < typename T >
double
laplace_rp_cb ( const Polynom & phi_i,
                const Polynom & phi_j,
                const Triangle & trk,
                int z,
                const Mesh & m,
                int point_i,
                int point_j,
                int, int,
                const T * u)
{
	double r = -laplace (phi_i, phi_j, trk, z);
	
	if (m.is_boundary(point_j)) {
		double b = 0;
		b += integrate(diff(diff(phi_j, 0) * phi_i, 0), trk, z);
		b += integrate(diff(diff(phi_j, 1) * phi_i, 1), trk, z);
		r += b;
	}
	
	r *= u[point_j];
	return r;
}

}

template < typename T, typename Matrix >
Laplace < T, Matrix >::Laplace (const Mesh & m) : 
	m_ (m),
	idt_full_ ( m.size ), 
	idt_ ( m.inner_size ), 
	laplace_ ( m.inner_size ),
	bnd2_ ( m.inner_size ),
	bnd3_ ( m.inner_size ),
	d_(m)
{
	generate_full_matrix (idt_full_, m, Laplace_Private::id_cb, (void*) 0);

	generate_matrix (idt_, m, Laplace_Private::id_cb, (void*) 0);
	generate_matrix (laplace_, m, Laplace_Private::laplace_integrate_cb, (void*) 0);
	generate_boundary_matrix (bnd2_, m_, Laplace_Private::laplace_bnd2_cb, (void*) 0);
	generate_boundary_matrix (bnd3_, m_, Laplace_Private::laplace_integrate_cb, (void*) 0);
}

template < typename T, typename Matrix >
void Laplace < T, Matrix >::solve (T * Ans, const T * F, const T * bnd)
{
	using namespace phelm;
	//пока используем первый порядок
	int sz  = (int) m_.ps.size();
	int ntr = (int) m_.tr.size();
	int rs  = (int) m_.inner.size();    //размерность

	Array b (rs);     // правая часть
	Array x (rs);     // ответ

	u2p (&x[0], F, m_);
	idt_.mult_vector (&b[0], &x[0]);
	if (bnd)
	{
		bnd2_.mult_vector (&x[0], bnd);
		vec_sum (&b[0], &b[0], &x[0], (int) x.size() );
	}

	phelm::solve (Ans, bnd, &b[0], laplace_, m_);
}

template < typename T, typename Matrix >
void Laplace < T, Matrix > ::calc2 (T * Ans, const T * F)
{
	using namespace phelm;
	int rs = (int) m_.inner.size();
	int os = (int) m_.outer.size();
	int sz = (int) m_.ps.size();

	Array in (rs);
	Array out (rs);
	Array tmp (os);
	Timer t;
	u2p (&in[0], F, m_);

	proj_bnd (&tmp[0], F, m_);
	laplace_.mult_vector (&out[0], &in[0]);
	bnd3_.mult_vector (&in[0], &tmp[0]);
	vec_sum (&out[0], &out[0], &in[0], (int) in.size() );
	idt_.solve (Ans, &out[0]);
}

/*
 * Оператор Лапласа на границе не определен, поэтому вставляйте сюда
 * границу только если вы знаете, что делаете!
 *
 * Если этот оператор Лапласа входит в праву часть уравнения, то
 * напишите 0 вместо границы.
 */
template < typename T, typename Matrix >
void Laplace < T, Matrix > ::calc1 (T * Ans, const T * F, const T * bnd)
{
	/*
	Array dx(m_.size);
	Array dy(m_.size);

	d_.calc_x(&dx[0], F);
	d_.calc_x(&dx[0], &dx[0]);
	d_.calc_y(&dy[0], F);
	d_.calc_y(&dy[0], &dy[0]);
	vec_sum(Ans, &dx[0], &dy[0], m_.size);
	*/

	Array out (m_.inner.size() );
	calc2 (&out[0], F);
	p2u (Ans, &out[0], bnd, m_);

	/*
	Array rp(m_.size);
	generate_full_right_part(&rp[0], m_, 
		Laplace_Private::laplace_rp_cb<T>, F);

	idt_full_.print(stderr);
	idt_full_.solve(&Ans[0], &rp[0]);
	*/
}

#endif /* LAPL_PRIVATE_H */

