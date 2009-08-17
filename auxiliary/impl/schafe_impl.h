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
 
/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 */

namespace SphereChafe_Private
{

double 
schafe_integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
					 int, int, 
                     SphereChafeConfig * d);

double f(/*double u,*/ double x, double y, double t, 
				double mu, double sigma);


template < typename T >
struct schafe_right_part_cb_data
{
	const T * F;
	const T * bnd;
	SphereChafeConfig  * d;
};

template < typename T >
double 
schafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
					  int i, int j,
                      schafe_right_part_cb_data < T > * d)
{
	const T * F = d->F;
	double b;

	if (m.ps_flags[point_j] == 1) { // íà ãðàíèöå
		int j0       = m.p2io[point_j]; //íîìåð âíåøíåé òî÷êè
		const T * bnd = d->bnd;
		b = - (double)bnd[j0] * schafe_integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, i, j, d->d);
		//b = 0.0;
	} else {
		b = (double)F[point_j] * integrate_cos(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

}

template < typename T >
SphereChafe < T > ::SphereChafe(const Mesh & m, double tau, double sigma, double mu)
	: SphereChafeConfig(tau, sigma, mu), m_(m), laplace_(m), A_((int)m.inner.size())
{
	/* Ìàòðèöà ëåâîé ÷àñòè */
	/* îïåðàòîð(u) = u/dt-mu \Delta u/2 + sigma u/2*/

	generate_matrix(A_, m_, SphereChafe_Private::schafe_integrate_cb, this);
}

/*
 * \f$\frac{du}{dt} = \mu \delta u - \sigma u + f (u)\f$
 */
template < typename T >
void SphereChafe < T > ::solve(T * Ans, const T * X0,
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

	// ãåíåðèðóåì ïðàâóþ ÷àñòü
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	u2p(&u[0], X0, m_);
	laplace_.calc2(&delta_u[0], X0);

	// u/dt + mu \Delta u / 2
	vec_sum1(&delta_u[0], &u[0], &delta_u[0], (T)(1.0 / tau_), (T)(mu_ * 0.5), rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vec_sum1(&delta_u[0], &delta_u[0], &u[0], (T)(1.0), (T)(-sigma_ * 0.5), rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		int point = m_.inner[i];
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();

		rp[i] = (T) SphereChafe_Private::f(/*u[i],*/ x, y, t, mu_, sigma_);
	}
	vec_copy_from_host(&crp[0], &rp[0], rs);
	vec_sum(&u[0], &delta_u[0], &crp[0], rs);

	// ïðàâóþ ÷àñòü íà ãðàíèöå íå çíàåì !!!
	p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	vec_copy_from_device(&hp[0], &p[0], sz);

	SphereChafe_Private::schafe_right_part_cb_data < T > data2;
	data2.F   = &hp[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		SphereChafe_Private::schafe_right_part_cb < T >, &data2);

//	fprintf(stderr, "rp: \n");vector_print(&delta_u[0], rs);
//	fprintf(stderr, "matrix:\n");A_.print();
//	laplace_.print();

	vec_copy_from_host(&crp[0], &rp[0], rs);
	phelm::solve(Ans, bnd, &crp[0], A_, m_);
}
