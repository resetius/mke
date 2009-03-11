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

#include "barvortex.h"

using namespace std;

static double 
jacobian(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, 
	 const Mesh & m, int i, int j, void * data)
{
	double pt1 = integrate_cos(diff(phi_i, 1) * diff(phi_j, 0), trk, m.ps);
	double pt2 = integrate_cos(diff(phi_i, 0) * diff(phi_j, 1), trk, m.ps);
	return pt1 - pt2;
}

/**
 * J(u,v)=1/cos(phi) (du/d\la dv/d\phi - du/d\phi dv/d\la)
 */
Jacobian::Jacobian(const Mesh & m): m_(m)
{
}

void Jacobian::calc1(double * Ans, const double * u, const double * v, const double * bnd)
{
	vector < double > p1(m_.inner.size());
	calc2(&p1[0], u, v);
	mke_p2u(Ans, &p1[0], bnd, m_);
}

void Jacobian::calc2(double * Ans, const double * u, const double * v)
{
	int sz = m_.ps.size();
	vector < double > j1(sz);
	convolution(&j1[0], u, v, m_, (scalar_cb_t)jacobian, 0);
	mke_u2p(Ans, &j1[0], m_);
}

BarVortex::BarVortex(const Mesh & m): m_(m), l_(m), j_(m)
{
}

void BarVortex::calc(double * Ans, const double * F, const double * bnd, double t)
{
	assert(0);
}

