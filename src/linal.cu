/* -*- charset: utf-8 -*- */
/*$Id: linal_base.cpp,v b772905ca13d 2009/08/07 07:39:47 aozeritsky $*/

/* Copyright (c) 2009 Alexey Ozeritsky (������� ���������)
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

__global__ void _sparse_mult_vector_ld(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n)
{
}

void sparse_mult_vector_ld(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n)
{
	_sparse_mult_vector_ld <<< 4, 4 >>> (r, Ap, Ai, Ax, x, n);
}

__global__ void _sparse_mult_vector_lf(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n)
{
}

void sparse_mult_vector_lf(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n)
{
	_sparse_mult_vector_lf <<< 4, 4 >>> (r, Ap, Ai, Ax, x, n);
}
