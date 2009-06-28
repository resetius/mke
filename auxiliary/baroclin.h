#ifndef BARVORTEX_H
#define BARVORTEX_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/**
 * @file baroclin.h
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * Copyright (c) 2009 Alexey Ozeritsky (Алексей Озерицкий)
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

#include <vector>

#include "laplace.h"
#include "solver.h"
#include "jacobian.h"

class Baroclin: public SphereNorm {
public:
	typedef double (*rp_t ) (double phi, double lambda, double t,
	                double sigma, double mu, double sigma1,
					double mu1, double alpha);
        typedef double (*coriolis_t) (double phi, double lambda);

private:
	const Mesh & m_;
	SphereLaplace l_;
	SphereJacobian j_;
	Matrix A_;
	Matrix Ab_;   // for backward

	std::vector < double > lh_; // l + h

public:
	double tau_;
	double sigma_;
	double mu_;
	double sigma1_;
	double mu1_;
	double alpha_;
	double theta_;

private:
	rp_t f_;
	rp_t g_;
	coriolis_t coriolis_;

public:
	Baroclin(const Mesh & m, rp_t f, rp_t g, 
		coriolis_t coriolis, double tau, 
		double sigma, double mu, 
		double sigma1, double mu1, double alpha);

	/**
	 * u1  -- ответ
	 * u2  -- ответ
	 * u11  -- значение функции на предыдущем шаге по времени
	 * u21  -- значение функции на предыдущем шаге по времени
	 * bnd -- граничное условие
	 * t   -- время
	 */
	void calc(double * u11,  double * u21, 
		const double * u1, const double * u2, 
		const double * bnd, double t);

	void calc_L(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2,
		const double * bnd, double t);

	void calc_L_1(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2,
		const double * bnd, double t);

	void calc_LT(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2,
		const double * bnd, double t);

	void L_step(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2);

	void L_1_step(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2);

	void LT_step(double * u11, double * u21, 
		const double * u1, const double * u2,
		const double * z1, const double * z2);
};

#endif /* BARVORTEX_H */

