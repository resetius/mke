/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2015 Alexey Ozeritsky
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

#include "mesh.h"

using namespace phelm;

const Triangle::basis_t & Triangle::elem1(int zone) const
{
	if (zone < 0) {
		zone = z;
	}

	if ((int)phik.size() < zone + 1) {
		phik.resize((size_t)(zone + 1));
	}

	if (phik[zone].empty()) {
		phik[zone] = prepare_basis(zone);
	}

	return phik[zone];
}

const Polynom & Triangle::elem1(int p1, int zone) const
{
	const basis_t & phik = elem1(zone);

	if (p1 == p[0])
	{
		return phik[0];
	} 
	else if (p1 == p[1])
	{
		return phik[1];
	} 
	else if (p1 == p[2])
	{
		return phik[2];
	}
	abort();
	return *((Polynom*)0);
}

Triangle::basis_t Triangle::prepare_basis(int z) const
{
	basis_t r;
	r.reserve(3);

	// p0
	r.push_back((P2X - x(1, z)) * (y(2, z) - y(1, z))
		- (P2Y - y(1, z)) * (x(2, z) - x(1, z)));
	// p1
	r.push_back((P2X - x(0, z)) * (y(2, z) - y(0, z))
		- (P2Y - y(0, z)) * (x(2, z) - x(0, z)));
	// p2
	r.push_back((P2X - x(0, z)) * (y(1, z) - y(0, z))
		- (P2Y - y(0, z)) * (x(1, z) - x(0, z)));

	for (uint i = 0; i < 3; ++i)
	{
		r[i] /= r[i].apply(x(i, z), y(i, z));
	}

	/*
	x[0] = X (0, ps);
	y[0] = Y (0, ps);
	x[1] = X (1, ps);
	y[1] = Y (1, ps);
	x[2] = X (2, ps);
	y[2] = Y (2, ps);
	*/
	return r;
}
