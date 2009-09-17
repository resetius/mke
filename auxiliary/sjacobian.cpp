/* -*- charset: utf-8 -*- */
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


#include <vector>
#include <math.h>

#include "phelm.h"
#include "sjacobian.h"
#include "util.h"
#include "ver.h"

VERSION("$Id$");

using namespace std;
using namespace phelm;

static double 
jacobian(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, 
	 const Mesh & m, int i, int j, void * data)
{
	Polynom pt1 = diff(phi_i, 1) * diff(phi_j, 0);
	Polynom pt2 = diff(phi_i, 0) * diff(phi_j, 1);

	Point p = m.ps[i].p[0];
	return (pt1.apply(p.x, p.y) - pt2.apply(p.x, p.y)) / cos(p.x);
}

static double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int i,
		int j,
		void *)
{
	return integrate_cos(phi_i * phi_j, trk, m.ps);
}

static double diff_1_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i, int j,
		int, int,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate(poly, trk, m.ps);
	return r;
}

static double diff_1_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i, int j,
		int, int,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate_cos(poly, trk, m.ps);
	return r;
}

static double diff_2_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i, int j,
		int, int,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate(poly, trk, m.ps);
	return r;
}

static double diff_2_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i, int j,
		int, int,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate_cos(poly, trk, m.ps);
	return r;
}

void SphereJacobian::calc2(double * Ans, const double * u, const double * v)
{
	int rs = (int)m_.inner.size();
	int sz = (int)m_.ps.size();
	int os = (int)m_.outer.size();
#if 0
	vector < double > rp(sz);
	convolution(&rp[0], u, v, m_, (scalar_cb_t)jacobian, 0);
	u2p(Ans, &rp[0], m_);
#endif

#if 1
	vector < double > v_in(rs);
	vector < double > u_in(rs);
	vector < double > v_in_bnd(os);
	vector < double > u_in_bnd(os);

	vector < double > rp(rs);

	vector < double > pt1(rs);
	vector < double > pt2(rs);
	vector < double > tmp(rs);

	u2p(&v_in[0], v, m_);
	u2p(&u_in[0], u, m_);
	proj_bnd(&v_in_bnd[0], v, m_);
	proj_bnd(&u_in_bnd[0], u, m_);

	//generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_cos_rp, (void*)u);
	diff2_cos_.mult_vector(&tmp[0], &u_in[0]);
	diff2_cos_rp_.mult_vector(&rp[0], &u_in_bnd[0]);
	vec_sum(&rp[0], &rp[0], &tmp[0], (int)rp.size());
	idt_.solve(&pt1[0], &rp[0]);

	//generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_rp, (void*)v);
	diff1_.mult_vector(&tmp[0], &v_in[0]);
	diff1_rp_.mult_vector(&rp[0], &v_in_bnd[0]);
	vec_sum(&rp[0], &rp[0], &tmp[0], (int)rp.size());
	idt_.solve(&tmp[0], &rp[0]);
	vec_mult(&pt1[0], &pt1[0], &tmp[0], (int)pt1.size());

	//generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_cos_rp, (void*)u);
	diff1_cos_.mult_vector(&tmp[0], &u_in[0]);
	diff1_cos_rp_.mult_vector(&rp[0], &u_in_bnd[0]);
	vec_sum(&rp[0], &rp[0], &tmp[0], (int)rp.size());
	idt_.solve(&pt2[0], &rp[0]);

	//generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_rp, (void*)v);
	diff2_.mult_vector(&tmp[0], &v_in[0]);
	diff2_rp_.mult_vector(&rp[0], &v_in_bnd[0]);
	vec_sum(&rp[0], &rp[0], &tmp[0],(int) rp.size());
	idt_.solve(&tmp[0], &rp[0]);
	vec_mult(&pt2[0], &pt2[0], &tmp[0], (int)pt1.size());

	vec_diff(Ans, &pt1[0], &pt2[0], (int)pt1.size());
#endif
}

void SphereJacobian::calc2t(double * Ans, const double * u, const double * v)
{
	calc2(Ans, u, v);
	vec_mult_scalar(Ans, Ans, -1.0, (int)m_.inner.size());
}

/*
 * J(u,v)=1/cos(phi) (du/d\la dv/d\phi - du/d\phi dv/d\la)
 */
SphereJacobian::SphereJacobian(const Mesh & m): m_(m), 
	idt_((int)m.inner.size()),
	diff1_((int)m.inner.size()),
	diff2_((int)m.inner.size()),
	diff1_cos_((int)m.inner.size()),
	diff2_cos_((int)m.inner.size()),
	diff1_rp_((int)m.inner.size()),
	diff2_rp_((int)m.inner.size()),
	diff1_cos_rp_((int)m.inner.size()),
	diff2_cos_rp_((int)m.inner.size()),

	diff1_t_((int)m.inner.size()),
	diff2_t_((int)m.inner.size()),
	diff1_cos_t_((int)m.inner.size()),
	diff2_cos_t_((int)m.inner.size()),
	diff1_rp_t_((int)m.inner.size()),
	diff2_rp_t_((int)m.inner.size()),
	diff1_cos_rp_t_((int)m.inner.size()),
	diff2_cos_rp_t_((int)m.inner.size())
{
	generate_matrix(idt_, m, id_cb, (void*)0);

	generate_matrix(diff1_, m_, diff_1_rp, (double*)0);
	generate_matrix(diff2_, m_, diff_2_rp, (double*)0);
	generate_matrix(diff1_cos_, m_, diff_1_cos_rp, (double*)0);
	generate_matrix(diff2_cos_, m_, diff_2_cos_rp, (double*)0);

	generate_boundary_matrix(diff1_rp_, m_, diff_1_rp, (double*)0);
	generate_boundary_matrix(diff2_rp_, m_, diff_2_rp, (double*)0);
	generate_boundary_matrix(diff1_cos_rp_, m_, diff_1_cos_rp, (double*)0);
	generate_boundary_matrix(diff2_cos_rp_, m_, diff_2_cos_rp, (double*)0);

	generate_matrix(diff1_t_, m_, diff_1_rp, (double*)0);
	generate_matrix(diff2_t_, m_, diff_2_rp, (double*)0);
	generate_matrix(diff1_cos_t_, m_, diff_1_cos_rp, (double*)0);
	generate_matrix(diff2_cos_t_, m_, diff_2_cos_rp, (double*)0);

	generate_boundary_matrix(diff1_rp_t_, m_, diff_1_rp, (double*)0);
	generate_boundary_matrix(diff2_rp_t_, m_, diff_2_rp, (double*)0);
	generate_boundary_matrix(diff1_cos_rp_t_, m_, diff_1_cos_rp, (double*)0);
	generate_boundary_matrix(diff2_cos_rp_t_, m_, diff_2_cos_rp, (double*)0);
}

void SphereJacobian::calc1(double * Ans, const double * u, const double * v, const double * bnd)
{
	vector < double > p1(m_.inner.size());
	calc2(&p1[0], u, v);
	p2u(Ans, &p1[0], bnd, m_);
}

void SphereJacobian::calc1t(double * Ans, const double * u, const double * v, const double * bnd)
{
	vector < double > p1(m_.inner.size());
	calc2t(&p1[0], u, v);
	p2u(Ans, &p1[0], bnd, m_);
}

