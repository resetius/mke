#ifndef BARVORTEX_H
#define BARVORTEX_H
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
 *
 * @section DESCRIPTION
 * The Barvortex equation class.
 */

#include <vector>

#include "laplace.h"
#include "solver.h"
#include "jacobian.h"

/**
 * @ingroup aux
 * @{
 */

/**
 * Solve the barvortex equation.
 *
 \f[
 \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
 \f]
 * where \f$J(\cdot,\cdot)\f$ is SphereJacobian and \f$\Delta\f$ is SphereLaplace.
 * @see SphereJacobian, SphereLaplace
 */
class BarVortex: public SphereNorm {
public:
	/**
	 * function of right part.
	 */
	typedef double (*rp_t ) (double phi, double lambda, double t,
		double mu, double sigma);
	/**
	 * coriolis
	 */
	typedef double (*coriolis_t) (double phi, double lambda);

private:
	const Mesh & m_;
	SphereLaplace l_;
	SphereJacobian j_;
	Matrix A_;
	Matrix bnd_;

	Matrix Ab_;   // for backward
	Matrix bndb_; // for backward

	std::vector < double > lh_; // l + h
	std::vector < double > f_;  // f right part

public:
	double tau_;   ///< tau
	double sigma_; ///< sigma
	double mu_;    ///< theta
	double theta_; ///< параметр схемы от 0 до 1

private:
	rp_t rp_;
	coriolis_t coriolis_;

public:
	/**
	 * describe!
	 */
	BarVortex(const Mesh & m, rp_t rp, coriolis_t coriolis, double tau, double sigma, double mu);

	/**
 \f[
 \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi)
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
 \f]
	 * @param Ans - output value
         * @param F - значение функции на предыдущем шаге по времени
	 * @param bnd - граничное условие
	 * @param t   - время
	 */
	void calc(double * Ans, const double * F, const double * bnd, double t);


	/**
	 *
 \f[
 \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta z)  + J(z, \Delta \psi)
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = 0
 \f]
	 */
	void calc_L(double * Ans, const double * F, const double * z, const double * bnd, double t);
	/**
	 * describe!
	 */
	void calc_L_1(double * Ans, const double * F, const double * z, const double * bnd, double t);

	/**
	 * describe!
	 */
	void calc_LT(double * Ans, const double * F, const double * z, const double * bnd, double t);

	/**
	 * J(psi, L(z)) + J(z, L(psi)) + J(psi, l + h) + sigma L(psi) - mu LL(psi) 
	 */
	void L_spectr(double * u1, const double * u, const double * z, const double * bnd);
	/**
	 * describe!
	 */
	void LT_spectr(double * u1, const double * u, const double * z, const double * bnd);

	/**
	 * describe!
	 */
	void S_step(double * Ans, const double * F);
	/**
	 * describe!
	 */
	void L_step(double * Ans, const double * F, const double * z);
	/**
	 * describe!
	 */
	void L_1_step(double * Ans, const double * F, const double * z);
	/**
	 * describe!
	 */
	void LT_step(double * Ans, const double * F, const double * z);
};

/** @} */

#endif /* BARVORTEX_H */

