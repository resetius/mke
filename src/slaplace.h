#ifndef SLAPL_H
#define SLAPL_H
/* -*- charset: utf-8 -*- */
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
#include "solver.h"

using phelm::Mesh;
using phelm::Triangle;
using phelm::Polynom;

/**
 * Compute Spherical Laplacian.
 \f[
  \Delta \psi(\varphi, \lambda)= \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi
 \f]
 */
template < typename T, typename Matrix = phelm::Solver < T > >
class SphereLaplace
{
public:
	typedef phelm::ArrayDevice < T > ArrayDevice;
	typedef phelm::ArrayHost < T > ArrayHost;
	Matrix idt_;
	Matrix laplace_;
	Matrix bnd1_; // L^-1
	Matrix bnd2_; // L^-1
	Matrix bnd3_; // L
	const Mesh & m_;

public:
	/**
	 * Constructor.
	 * @param m - mesh
	 */
	SphereLaplace (const Mesh & m);

	/**
	 * Calculate the value of Laplace operator for F in the inner points of the mesh.
	 * Set the value of boundary points from bnd.
	 * @param Ans - the answer (the value of Laplace operator in all mesh points)
	 * @param F -vector F
	 * @param bnd - needed boundary condition
	 */
	void calc1 (T * Ans, const T * F, const T * bnd);

	/**
	 * Calculate the value of Laplace operator for F in the inner points of the mesh.
	 * Returns vector of inner points.
	 * @param Ans - the answer (the value of Laplace operator in the inner points)
	 * @param F -vector F
	 */
	void calc2 (T * Ans, const T * F);

	/**
	 * Solve Laplace equation.
	 \f{eqnarray*}
	   \Delta \psi &=& f(x, y) \\
	  \psi|_{\partial\Omega}&=&u_0
	  \f}
	 * @param Ans - the answer
	 * @param F - right part
	 * @param bnd - boundary condition
	 */
	void solve (T * Ans, const T * F, const T * bnd);
};

/**
 * Calculate one element of finite-element matrix
 * for spheric Laplace operator.
 * @param phi_i - basis function
 * @param phi_j - basis function
 * @param trk - integral is taken over that triangle
 * @param z - subdomain number
 */
double slaplace (const Polynom & phi_i, const Polynom & phi_j,
                 const Triangle & trk, int z);


#include "impl/slaplace_impl.h"

#endif /*SLAPL_H*/
