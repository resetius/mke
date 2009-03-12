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

#include "barvortex.h"
#include "util.h"

using namespace std;

static double 
jacobian(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, 
	 const Mesh & m, int i, int j, void * data)
{
	Polynom pt1 = diff(phi_i, 1) * diff(phi_j, 0);
	Polynom pt2 = diff(phi_i, 0) * diff(phi_j, 1);

	Point p = m.ps[i];
	return (pt1.apply(p.x, p.y) - pt2.apply(p.x, p.y)) / cos(p.x);
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

/**
 * J(u,v)=1/cos(phi) (du/d\la dv/d\phi - du/d\phi dv/d\la)
 */
Jacobian::Jacobian(const Mesh & m): m_(m), idt_(m.inner.size())
{
	generate_matrix(idt_, m, id_cb, 0);
}

void Jacobian::calc1(double * Ans, const double * u, const double * v, const double * bnd)
{
	vector < double > p1(m_.inner.size());
	calc2(&p1[0], u, v);
	mke_p2u(Ans, &p1[0], bnd, m_);
}

static double diff_1_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double r = u[j] * integrate(poly, trk, m.ps);
	return r;
}

static double diff_1_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double r = u[j] * integrate_cos(poly, trk, m.ps);
	return r;
}

static double diff_2_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double r = u[j] * integrate(poly, trk, m.ps);
	return r;
}

static double diff_2_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double r = u[j] * integrate_cos(poly, trk, m.ps);
	return r;
}

void Jacobian::calc2(double * Ans, const double * u, const double * v)
{
	int rs = m_.inner.size();
	int sz = m_.ps.size();
#if 0
	vector < double > rp(sz);
	convolution(&rp[0], u, v, m_, (scalar_cb_t)jacobian, 0);
	mke_u2p(Ans, &rp[0], m_);
#endif

#if 1
	vector < double > rp(rs);
	vector < double > pt1(rs);
	vector < double > pt2(rs);
	double * tmp = Ans;

	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_cos_rp, (void*)u);
	idt_.solve(&pt1[0], &rp[0]);
	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_rp, (void*)v);
	idt_.solve(&tmp[0], &rp[0]);
	vector_mult(&pt1[0], &pt1[0], &tmp[0], pt1.size());

	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_cos_rp, (void*)u);
	idt_.solve(&pt2[0], &rp[0]);
	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_rp, (void*)v);
	idt_.solve(&tmp[0], &rp[0]);
	vector_mult(&pt2[0], &pt2[0], &tmp[0], pt1.size());

	vector_diff(Ans, &pt1[0], &pt2[0], pt1.size());
#endif
}

BarVortex::BarVortex(const Mesh & m): m_(m), l_(m), j_(m)
{
}

/**
 * d L(phi)/dt + J(phi, L(phi)) + J(phi, l + h) + sigma L(phi) - mu LL(phi) = f(phi, la)
 * L = Laplace
 */
void BarVortex::calc(double * Ans, const double * F, const double * bnd, double t)
{
	assert(0);
}

