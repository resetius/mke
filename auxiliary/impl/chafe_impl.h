/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky
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

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0         = m.p2io[point_j]; //номер внешней точки
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
	/* Матрица левой части */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/
	generate_matrix(A_, m_, Chafe_Private::chafe_integrate_cb, this);
}

/* \frac{du}{dt} = \mu \Delta u - \sigma u + f (u) */
template < typename T >
void Chafe < T > ::solve(T * Ans, const T * X0,
						const T * bnd, double t)
{
	int rs  = (int)m_.inner.size();
	int os = (int)m_.outer.size();
	int sz  = (int)m_.ps.size();
	ArrayDevice u(rs);
	ArrayDevice p(sz);
	ArrayHost   hp(sz);
	ArrayDevice delta_u(rs);

	ArrayHost   rp(rs);
	ArrayDevice crp(rs);
	ArrayHost   hbnd(os); //TODO: remove
	if (bnd)
		vec_copy_from_device(&hbnd[0], &bnd[0], os);

	// генерируем правую часть
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	u2p(&u[0], X0, m_);
	laplace_.calc2(&delta_u[0], X0);

	// u/dt + mu \Delta u / 2
	vec_sum1(&delta_u[0], &u[0], &delta_u[0], (T)(1.0 / tau_), (T)(mu_ * 0.5), rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vec_sum1(&delta_u[0], &delta_u[0], &u[0], (T)1.0, (T)(-sigma_ * 0.5), rs);

	// TODO: тут надо сделать метод простой итерации
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

	// правую часть на границе не знаем !
	p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	vec_copy_from_device(&hp[0], &p[0], sz);

	// TODO: убрать это !
	Chafe_Private::chafe_right_part_cb_data < T > data2;
	data2.F   = &hp[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		Chafe_Private::chafe_right_part_cb < T >, &data2);

	vec_copy_from_host(&crp[0], &rp[0], rs);
	phelm::solve(Ans, bnd, &crp[0], A_, m_);
}
