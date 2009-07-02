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
 * The two-dimensional baroclinic atmosphere equations class.
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
 * Solve the two-dimensional baroclinic atmosphere equations.
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-& \\
	- \alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
 * where \f$J(\cdot,\cdot)\f$ is SphereJacobian and \f$\Delta\f$ is SphereLaplace.
 * @see SphereJacobian, SphereLaplace
 */
class Baroclin: public SphereNorm {
public:
	/**
	 * Function of right part.
	 * @param phi - \f$\varphi\f$
	 * @param lambda - \f$\lambda\f$
	 * @param t - time
	 * @param sigma - \f$\sigma\f$
	 * @param mu - \f$\mu\f$
	 * @param sigma1 - \f$\sigma_1\f$
	 * @param mu1 - \f$\mu_1\f$
	 * @param alpha - \f$\alpha\f$
	 */
	typedef double (*rp_t ) (double phi, double lambda, double t,
	                double sigma, double mu, double sigma1,
					double mu1, double alpha);
	/**
	 * Coriolis.
	 * l+h: coriolis function plus orographic function.
	 * @param phi - \f$\varphi\f$
	 * @param lambda - \f$\lambda\f$
	 */
	typedef double (*coriolis_t) (double phi, double lambda);

private:
	const Mesh & m_;
	SphereLaplace l_;
	SphereJacobian j_;
	Matrix A_;
	Matrix Ab_;   // for backward

	std::vector < double > lh_; // l + h

public:
	double tau_;   ///<time step
	double sigma_; ///<\f$\sigma\f$
	double mu_;    ///<\f$\mu\f$
	double sigma1_;///<\f$\sigma\f$
	double mu1_;   ///<\f$\mu_1\f$
	double alpha_; ///<\f$\alpha\f$

	/**
	 * Time discretization scheme parameter \f$\theta\f$.
	 * The default value is 0.5 (Crank–Nicolson).
	 */
	double theta_; 

private:
	rp_t f_;
	rp_t g_;
	coriolis_t coriolis_;

public:
	/**
	 * Constructor.
	 * @param m - mesh
	 * @param f - function of right part (\f$f\f$)
	 * @param g - function of right part (\f$g\f$)
	 * @param coriolis - h+l
	 * @param tau - time step
	 * @param sigma - \f$\sigma\f$
	 * @param mu - \f$\mu\f$
	 * @param sigma1 - \f$\sigma_1\f$
	 * @param mu1 - \f$\mu_1\f$
	 * @param alpha - \f$\alpha\f$
	 */
	Baroclin(const Mesh & m, rp_t f, rp_t g, 
		coriolis_t coriolis, double tau, 
		double sigma, double mu, 
		double sigma1, double mu1, double alpha);

	/**
	 * Solve the two-dimensional baroclinic atmosphere equations.
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	- \alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
	 * @param u1  - output vector
	 * @param u2  - output vector
	 * @param u11 - input vector (previous time step)
	 * @param u21 - input vector (previous time step)
	 * @param bnd - boundary condition
	 * @param t   - time
	 */
	void calc(double * u11,  double * u21, 
		const double * u1, const double * u2, 
		const double * bnd, double t);

	/**
	 * Solve the linearized two-dimensional baroclinic atmosphere equations in a neibourhood of point (z1, z2)
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta z_1 + l + h) + J(z_1, \Delta u_1)
  + J(z_2, \Delta u_2) + J(u_2, \Delta z_2) +
  \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& 0\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta z_2) + J(z_1, \Delta u_2)
  + J(u_2, \Delta z_1 + l + h) + J(z_2, \Delta u_1)
  + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	- \alpha^2 (\frac{\partial u_2}{\partial t} 
    + J(z_1, u_2) + J(u_1, z_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2) &=& 0,
 \f}
 
	 * @param u1  - output vector
	 * @param u2  - output vector
	 * @param u11 - input vector (previous time step)
	 * @param u21 - input vector (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 * @param bnd - boundary condition
	 * @param t   - time
	 */
	void calc_L(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2,
		const double * bnd, double t);


	/**
	 * Solve the invert linearized two-dimensional baroclinic atmosphere equations in a neibourhood of point (z1, z2)
	 * @param u1  - output vector
	 * @param u2  - output vector
	 * @param u11 - input vector (previous time step)
	 * @param u21 - input vector (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 * @param bnd - boundary condition
	 * @param t   - time
	 */
	void calc_L_1(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2,
		const double * bnd, double t);

	/**
	 * Adjoint operator to calc_L.
	 * (L u, v) = (u, LT v)
	 * @todo implement!
	 * @param u1  - output value
	 * @param u2  - output value
	 * @param u11 - input value (previous time step)
	 * @param u21 - input value (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 * @param bnd - boundary condition
	 * @param t   - time
	 */
	void calc_LT(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2,
		const double * bnd, double t);

	/**
	 * Helper-function. 
	 * Call calc_L(u11, u21, u1, u2, z1, z2, 0, 0).
	 * @param u1  - output value
	 * @param u2  - output value
	 * @param u11 - input value (previous time step)
	 * @param u21 - input value (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 */
	void L_step(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2);

	/**
	 * Helper-function.
	 * Call calc_L_1(u11, u21, u1, u2, z1, z2, 0, 0).
	 * @param u1  - output value
	 * @param u2  - output value
	 * @param u11 - input value (previous time step)
	 * @param u21 - input value (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 */
	void L_1_step(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2);

	/**
	 * Helper-function.
	 * Call calc_L(u11, u21, u1, u2, z1, z2, 0, 0).
	 * @param u1  - output value
	 * @param u2  - output value
	 * @param u11 - input value (previous time step)
	 * @param u21 - input value (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 */
	void LT_step(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2);
};

/** @} */

#endif /* BARVORTEX_H */

