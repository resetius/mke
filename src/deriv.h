#ifndef DERIV_H
#define DERIV_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2011 Alexey Ozeritsky
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

#include "solver.h"
#include "mesh.h"

using phelm::Mesh;

/**
 * @ingroup aux
 * @{
 */

template < typename T >
class Deriv
{
//	typedef phelm::Solver < T > Matrix;

	typedef linal::SparseSolver < T, 
		linal::StoreCSR < T , linal::Allocator > , 
		linal::StoreCSR < T , linal::Allocator > > Matrix;

	const Mesh & m_;

	Matrix diff_x_;
	Matrix diff_y_;

public:
	/**
	 * Constructor.
	 * @param m - mesh
	 */
	Deriv (const Mesh & m);

	/**
	 * Calculate derivative u'x
	 * @param Ans - the answer
	 * @param u - vector u
	 */
	void calc_x(T * Ans, const T * u);

	/**
	 * Calculate derivative u'y
	 * @param Ans - the answer
	 * @param u - vector u
	 */
	void calc_y(T * Ans, const T * u);
};

/** @} */

#include "impl/deriv_impl.h"

#endif /* DERIV_H */
