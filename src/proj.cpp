/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (јлексей ќзерицкий)
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

#include "phelm.h"
#include "ver.h"

VERSION("$Id$");

using namespace std;

namespace phelm {

/* добавл€ем краевые услови€ */
template < typename T >
void p2u_(T * u, const T * p, const T * bnd, const Mesh & m)
{
	int sz  = (int)m.ps.size();

#pragma omp parallel for
	for (int i = 0; i < sz; ++i) {
		if (m.ps_flags[i] == 1) {
			//внешн€€
			if (bnd) {
				u[i] = bnd[m.p2io[i]];
			} else {
				u[i] = 0;
			}
		} else {
			//внутренн€€
			u[i] = p[m.p2io[i]];
		}
	}
}

void p2u(double * u, const double * p, const double * bnd, const Mesh & m)
{
	p2u_(u, p, bnd, m);
}

void p2u(float * u, const float * p, const float * bnd, const Mesh & m)
{
	p2u_(u, p, bnd, m);
}

/* убираем краевые услови€ */
template < typename T >
void u2p_(T * p, const T * u, const Mesh & m)
{
	int rs = (int)m.inner.size();
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		p[i] = u[m.inner[i]];
	}
}

void u2p(double * p, const double * u, const Mesh & m)
{
	u2p_(p, u, m);
}

void u2p(float * p, const float * u, const Mesh & m)
{
	u2p_(p, u, m);
}

void set_bnd(double *u, const double * bnd, const Mesh & m)
{
	int sz  = (int)m.ps.size();
	for (int i = 0; i < sz; ++i) {
		if (m.ps_flags[i] == 1) {
			if (bnd) {
				u[i] = bnd[m.p2io[i]];
			} else {
				u[i] = 0;
			}
		}
	}
}

template < typename T >
void proj_(T * F, const Mesh & mesh, f_xy_t f)
{
	int sz = (int)mesh.ps.size();
#pragma omp parallel for
	for (int i = 0; i < sz; ++i)
	{
		const Point & p = mesh.ps[i].p[0];
		F[i] = (T)f(p.x, p.y);
	}
}

void proj(double * F, const Mesh & mesh, f_xy_t f)
{
	proj_(F, mesh, f);
}

void proj(float * F, const Mesh & mesh, f_xy_t f)
{
	proj_(F, mesh, f);
}

template < typename T >
void proj_bnd_(T * F, const Mesh & m, f_xy_t f)
{
	int sz = (int)m.outer.size();
#pragma omp parallel for
	for (int i = 0; i < sz; ++i) {
		int p0 = m.outer[i];
		const Point & p = m.ps[p0].p[0];
		F[i] = (T)f(p.x, p.y);
	}
}

void proj_bnd(double * F, const Mesh & m, f_xy_t f)
{
	proj_bnd_(F, m, f);
}

void proj_bnd(float * F, const Mesh & m, f_xy_t f)
{
	proj_bnd_(F, m, f);
}

void proj_bnd(double * F, const double * F1, const Mesh & m)
{
	int sz = (int)m.outer.size();
#pragma omp parallel for
	for (int i = 0; i < sz; ++i) {
		int p0 = m.outer[i];
		F[i] = F1[p0];
	}
}

void proj_bnd(float * F, const float * F1, const Mesh & m)
{
	int sz = (int)m.outer.size();
#pragma omp parallel for
	for (int i = 0; i < sz; ++i) {
		int p0 = m.outer[i];
		F[i] = F1[p0];
	}
}

void proj(double * F, const Mesh & mesh, f_xyt_t f, double t)
{
	int sz = (int)mesh.ps.size();
#pragma omp parallel for
	for (int i = 0; i < sz; ++i)
	{
		F[i] = f(mesh.ps[i].x(), mesh.ps[i].y(), t);
	}
}

void proj_bnd(double * F, const Mesh & m, f_xyt_t f, double t)
{
	int sz = (int)m.outer.size();
#pragma omp parallel for
	for (int i = 0; i < sz; ++i) {
		int p0 = m.outer[i];
		const Point & p = m.ps[p0].p[0];
		F[i] = f(p.x, p.y, t);
	}
}

} /* namespace */
