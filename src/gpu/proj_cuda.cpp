/* -*- charset: utf-8 -*- */
/*$Id $*/

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

#include "phelm.h"
#include "ver.h"

VERSION("$Id$");

namespace phelm {

	void u2p(double * p, const double * u, const int * inner, int rs);
	void u2p(float * p, const float * u, const int * inner, int rs);

	void p2u(double * p, const double * u, const double * bnd,
		const int * p2io, const int * ps_flags, int sz);
	void p2u(float * p, const float * u, const float * bnd,
		const int * p2io, const int * ps_flags, int sz);

	void proj_bnd(double * F, const double * F1, 
		const int * outer, int os);
	void proj_bnd(float * F, const float * F1, 
		const int * outer, int os);

	void set_bnd(double * u, const double * bnd, 
		const int * p2io, const int * ps_flags, int sz);
	void set_bnd(float * u, const float * bnd, 
		const int * p2io, const int * ps_flags, int sz);

	void u2p(float * p, const float * u, const Mesh & m)
	{
		u2p(p, u, &m.d.inner[0], m.inner_size);
	}

	void u2p(double * p, const double * u, const Mesh & m)
	{
		u2p(p, u, &m.d.inner[0], m.inner_size);
	}

	void p2u(float * u, const float * p, const float * bnd, const Mesh & m)
	{
		p2u(u, p, bnd, &m.d.p2io[0], &m.d.ps_flags[0], m.size);
	}

	void p2u(double * u, const double * p, const double * bnd, const Mesh & m)
	{
		p2u(u, p, bnd, &m.d.p2io[0], &m.d.ps_flags[0], m.size);
	}

	void proj_bnd(double * F, const double * F1, const Mesh & m)
	{
		proj_bnd(F, F1, &m.d.outer[0], m.outer_size);
	}

	void proj_bnd(float * F, const float * F1, const Mesh & m)
	{
		proj_bnd(F, F1, &m.d.outer[0], m.outer_size);
	}

	void set_bnd(double *u, const double * bnd, const Mesh & m)
	{
		set_bnd(u, bnd, &m.d.p2io[0], &m.d.ps_flags[0], m.size);
	}

	void set_bnd(float *u, const float * bnd, const Mesh & m)
	{
		set_bnd(u, bnd, &m.d.p2io[0], &m.d.ps_flags[0], m.size);
	}
} /* namespace */
