#ifndef CHAFE_H
#define CHAFE_H
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

#include <vector>

#include "laplace.h"

/**
 * Chafe-Infante equation on flat domain.
 * Configuration.
 */
struct ChafeConfig
{
	double tau_;   ///<time step
	double mu_;    ///<\f$\mu\f$
	double sigma_; ///<\f$\sigma\f$

	ChafeConfig (double tau, double sigma, double mu) :
			tau_ (tau), mu_ (mu), sigma_ (sigma) {}
};

/**
 * Chafe-Infante equation on flat domain.
 * \f$ \frac{du}{dt} = \mu \Delta u - \sigma u + f (u)\f$, where
 * \f$Delta\f$ is Laplace. @see Laplace
 */
template < typename T >
class Chafe: public ChafeConfig
{
private:
	typedef phelm::Solver < T > Matrix;
	typedef phelm::ArrayDevice < T > ArrayDevice;
	typedef phelm::ArrayHost < T > ArrayHost;
	const Mesh & m_;
	Laplace < T > laplace_; /* Лапласиан */
	Matrix A_;        /* Матрица левой части */
	Matrix bnd_;
	ArrayPool < ArrayDevice > idp_; //inner device array pool
	ArrayPool < ArrayHost > idh_; //inner host array pool

public:
	/**
	 * Constructor.
	 * @param m - mesh
	 * @param tau - time step
	 * @param sigma - \f$\sigma\f$
	 * @param mu - \f$\mu\f$
	 */
	Chafe (const Mesh & m, double tau, double sigma, double mu);
	~Chafe() {}

	/**
	 * Solve Chafe-Infante equation on flat domain.
	 * \f$\frac{du}{dt} = \mu \Delta u - \sigma u + f (u)\f$
	 * @param Ans - output vector
	 * @param X0 - intput vector (previous time step)
	 * @param bnd - boundary condition
	 * @param t - time
	 */
	void solve (T * Ans, const T * X0,
	            const T * bnd, double t);
};

#include "chafe_impl.h"

#endif /* CHAFE_H */
