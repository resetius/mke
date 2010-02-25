#ifndef NORM_H
#define NORM_H
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

/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 */

#include <vector>

#include "phelm.h"
#include "solver.h"

using phelm::Mesh;

/**
 * A fast spherical norm calculator.
 */
template < typename T >
class SphereNorm
{
public:
	typedef phelm::Solver < T > Matrix;
	const Mesh & m_; ///< mesh
	Matrix NORM_;    ///< for fast norm calculation

	/**
	 * Constructor.
	 * @param m - mesh
	 */
	SphereNorm (const Mesh & m) : m_ (m), NORM_ ( (int) m_.ps.size() )
	{
		phelm::generate_full_matrix (NORM_, m_, phelm::sphere_scalar_cb, (void*) 0);
	}

	/**
	 * Distance between two vectors.
	 * @param u - input vector
	 * @param v - input vector
	 * @return distance between u and v
	 */
	T dist (const T * u, const T * v)
	{
		return phelm::fast_dist (u, v, m_, NORM_);
	}

	/**
	 * Norm of vector.
	 * @param u - input vector
	 * @return norm of v
	 */
	T norm (const T * u)
	{
		return phelm::fast_norm (u, m_, NORM_);
	}

	/**
	 * Inner product of two vectors.
	 * @param u - input vector
	 * @param v - input vector
	 * @return inner product of u and v
	 */
	T scalar (const T * u, const T * v)
	{
		return phelm::fast_scalar (u, v, m_, NORM_);
	}
};

/**
 * A fast flat surface norm calculator.
 */
template < typename T >
class FlatNorm
{
public:
	typedef phelm::Solver < T > Matrix;

	const Mesh & m_; ///< mesh
	Matrix NORM_;    ///< for fast norm calculation

	/**
	 * Constructor.
	 * @param m - mesh
	 */
	FlatNorm (const Mesh & m) : m_ (m), NORM_ ( (int) m_.ps.size() )
	{
		phelm::generate_full_matrix (NORM_, m_,
		                             phelm::sphere_scalar_cb, (void*) 0);
	}

	/**
	 * Distance between two vectors.
	 * @param u - input vector
	 * @param v - input vector
	 * @return distance between u and v
	 */
	T dist (const T * u, const T * v)
	{
		return phelm::fast_dist (u, v, m_, NORM_);
	}

	/**
	 * Norm of vector.
	 * @param u - input vector
	 * @return norm of v
	 */
	T norm (const T * u)
	{
		return phelm::fast_norm (u, m_, NORM_);
	}

	/**
	 * Inner product of two vectors.
	 * @param u - input vector
	 * @param v - input vector
	 * @return inner product of u and v
	 */
	T scalar (const T * u, const T * v)
	{
		return phelm::fast_scalar (u, v, m_, NORM_);
	}
};

#endif /* NORM_H */
