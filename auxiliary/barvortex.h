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
 * The Barotropic vorticity equation class.
 */

#include <vector>

#include "slaplace.h"
#include "jacobian.h"
#include "norm.h"

/**
 * @ingroup aux
 * @{
 */

/**
 * Solve the Barotropic vorticity equation.
 *
 \f[
 \frac{\partial \Delta \varphi}{\partial t} + k_1 J(\psi, \Delta \psi) 
    + k_2 J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
 \f]
 * where \f$J(\cdot,\cdot)\f$ is spherical jacobian operator 
 * and \f$\Delta\f$ is spherical Laplace operator.
 * @see SphereJacobian, SphereLaplace
 */
class BarVortex: public SphereNorm < double > {
public:
	/**
	 * Function of right part.
	 * @param phi - \f$\varphi\f$
	 * @param lambda - \f$\lambda\f$
	 * @param t - time
	 * @param mu - \f$\mu\f$
	 * @param sigma - \f$\sigma\f$
	 */
	typedef double (*rp_t ) (double phi, double lambda, double t,
		double mu, double sigma);
	/**
	 * Coriolis.
	 * l+h: coriolis function plus orographic function.
	 * @param phi - \f$\varphi\f$
	 * @param lambda - \f$\lambda\f$
	 */
	typedef double (*coriolis_t) (double phi, double lambda);

private:
	const Mesh & m_;
	SphereLaplace < double > l_;
	SphereJacobian j_;
	Matrix A_;
	Matrix bnd_;

	Matrix Ab_;   // for backward
	Matrix bndb_; // for backward

	std::vector < double > lh_; // l + h
	std::vector < double > f_;  // f right part

public:
	double tau_;   ///< time step
	double sigma_; ///< \f$sigma\f$
	double mu_;    ///< \f$theta\f$
	double k1_;    ///< \f$k_1\f$
	double k2_;    ///< \f$k_2\f$

	/**
	 * Time discretization scheme parameter \f$\theta\f$.
	 * The default value is 0.5 (Crank–Nicolson).
	 */
	double theta_;

private:
	rp_t rp_;
	coriolis_t coriolis_;

public:
	/**
	 * Constructor.
	 * @param m - mesh
	 * @param rp - function of right part (\f$f\f$)
	 * @param coriolis - h+l
	 * @param tau - time step
	 * @param sigma - \f$\sigma\f$
	 * @param mu - \f$\mu\f$
	 * @param k1 - \f$\k_1\f$
	 * @param k2 - \f$\k_2\f$
	 */
	BarVortex(const Mesh & m, rp_t rp, coriolis_t coriolis, double tau, double sigma,
			  double mu, double k1, double k2);

	/**
	 * Write Parameters.
	 */
	void info();

	/**
	 * Solve the Barotropic vorticity equation.
 \f[
 \frac{\partial \Delta \varphi}{\partial t} + k_1 J(\psi, \Delta \psi)
    + k_2 J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
 \f]
	 * @param Ans - output value
	 * @param F - input vector (previous time step)
	 * @param bnd - boundary condition
	 * @param t   - time
	 */
	void calc(double * Ans, const double * F, const double * bnd, double t);


	/**
	 * Solve the linearized Barotropic vorticity equation in a neibourhood of point (z).
 \f[
 \frac{\partial \Delta \varphi}{\partial t} + k_1 J(\psi, \Delta z)  + k_1 J(z, \Delta \psi)
    + k_2 J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = 0
 \f]
	 * @param Ans - output vector
	 * @param F - input vector (previous time step)
	 * @param z - vector z
	 * @param bnd - boundary condition
	 * @param t - time
	 */
	void calc_L(double * Ans, const double * F, const double * z, const double * bnd, double t);
	/**
	 * Solve the invert linearized Barotropic vorticity equation in a neibourhood of point (z).
	 * @param Ans - output vector
	 * @param F - input vector (previous time step)
	 * @param z - vector z
	 * @param bnd - boundary condition
	 * @param t - time
	 */
	void calc_L_1(double * Ans, const double * F, const double * z, const double * bnd, double t);

	/**
	 * Adjoint operator to calc_L.
	 * (L u, v) = (u, LT v)
	 * @param Ans - output vector
	 * @param F - input vector (previous time step)
	 * @param z - vector z
	 * @param bnd - boundary condition
	 * @param t - time
	 */
	void calc_LT(double * Ans, const double * F, const double * z, const double * bnd, double t);

	/**
	 * Calculate operator.
	 *
	 \f[
	 k_1 J(\psi, \Delta z) + k_1 J(z, \Delta \psi) + k_2 J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi
	 \f]
	 */
	void L_spectr(double * u1, const double * u, const double * z, const double * bnd);
	/**
	 * Adjoint operator to L_spectr.
	 * (L u, v) = (u, LT v)
	 */
	void LT_spectr(double * u1, const double * u, const double * z, const double * bnd);

	/**
	 * Helper-function.
	 * Call calc(Ans, F, 0, 0).
	 * @param Ans - output vector
	 * @param F - input vector (previous time step)
	 */
	void S_step(double * Ans, const double * F);
	/**
	 * Helper-function.
	 * call calc_L(Ans, F, z, 0, 0).
	 * @param Ans - output vector
	 * @param F - input vector (previous time step)
	 * @param z - vector z
	 */
	void L_step(double * Ans, const double * F, const double * z);
	/**
	 * Helper-function.
	 * Call calc_L_1(Ans, F, z, 0, 0).
	 * @param Ans - output vector
	 * @param F - input vector (previous time step)
	 * @param z - vector z
	 */
	void L_1_step(double * Ans, const double * F, const double * z);
	/**
   	 * Helper-function.
	 * Call calc_LT(Ans, F, z, 0, 0).
	 * @param Ans - output vector
	 * @param F - input vector (previous time step)
	 * @param z - vector z
	 */
	void LT_step(double * Ans, const double * F, const double * z);
};

/** @} */

#endif /* BARVORTEX_H */

