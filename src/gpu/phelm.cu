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

#include "gpu/linal_cuda.h"
#include "ver.h"

VERSION("$Id$");

namespace phelm {

template < typename T >
__global__ void u2p_(T * p, const T * u, const int * inner, int rs)
{
	int threads = gridDim.x  * blockDim.x;
	int i  = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < rs; i += threads) {
		p[i] = u[inner[i]];
	}
}

__host__ void u2p(double * p, const double * u, const int * inner, int rs)
{
	SPLAY(rs);
	u2p_ <<< blocks, threads >>> (p, u, inner, rs);
}

__host__ void u2p(float * p, const float * u, const int * inner, int rs)
{
	SPLAY(rs);
	u2p_ <<< blocks, threads >>> (p, u, inner, rs);
}

template < typename T >
__global__ void p2u_(T * u, const T * p, const T * bnd, const int * p2io, 
	const int * ps_flags, int sz)
{
	int threads = gridDim.x  * blockDim.x;
	int i  = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < sz; i += threads) {
		if (ps_flags[i] == 1) {
			//внешн€€
			if (bnd) {
				u[i] = bnd[p2io[i]];
			} else {
				u[i] = 0;
			}
		} else {
			//внутренн€€
			u[i] = p[p2io[i]];
		}
	}
}

__host__ void p2u(double * u, const double * p, const double * bnd,
	const int * p2io, const int * ps_flags, int sz)
{
	SPLAY(sz);
	p2u_ <<< blocks, threads >>> (u, p, bnd, p2io, ps_flags, sz);
}

__host__ void p2u(float * u, const float * p, const float * bnd,
	const int * p2io, const int * ps_flags, int sz)
{
	SPLAY(sz);
	p2u_ <<< blocks, threads >>> (u, p, bnd, p2io, ps_flags, sz);
}

template < typename T >
__global__ void proj_bnd_(T * F, const T * F1, 
	const int * outer, int os)
{
	int threads = gridDim.x  * blockDim.x;
	int i  = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < os; i += threads) {
		int p0 = outer[i];
		F[i] = F1[p0];
	}
}

__host__ void proj_bnd(double * F, const double * F1, 
	const int * outer, int os)
{
	SPLAY(os);
	proj_bnd_ <<< blocks, threads >>> (F, F1, outer, os);
}

__host__ void proj_bnd(float * F, const float * F1, 
	const int * outer, int os)
{
	SPLAY(os);
	proj_bnd_ <<< blocks, threads >>> (F, F1, outer, os);
}

template < typename T >
__global__ void set_bnd_(T * u, const T * bnd, 
	const int * p2io, const int * ps_flags, int sz)
{
	int threads = gridDim.x  * blockDim.x;
	int i  = blockDim.x * blockIdx.x + threadIdx.x;
	for (;i < sz; i += threads) {
		if (ps_flags[i] == 1) {
			if (bnd) {
				u[i] = bnd[p2io[i]];
			} else {
				u[i] = 0;
			}
		}
	}
}

__host__ void set_bnd(double * u, const double * bnd, 
	const int * p2io, const int * ps_flags, int sz)
{
	SPLAY(sz);
	set_bnd_ <<< blocks, threads >>> (u, bnd, p2io, ps_flags, sz);
}

__host__ void set_bnd(float * u, const float * bnd, 
	const int * p2io, const int * ps_flags, int sz)
{
	SPLAY(sz);
	set_bnd_ <<< blocks, threads >>> (u, bnd, p2io, ps_flags, sz);
}

}
