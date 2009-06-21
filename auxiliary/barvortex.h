#ifndef BARVORTEX_H
#define BARVORTEX_H
/*$Id$*/

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

#include <vector>

#include "laplace.h"
#include "solver.h"
#include "jacobian.h"

class BarVortex: public SphereNorm {
public:
	typedef double (*rp_t ) (double phi, double lambda, double t,
		double mu, double sigma);
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
	double tau_;
	double sigma_;
	double mu_;
	double theta_; // �������� ����� �� 0 �� 1

private:
	rp_t rp_;
	coriolis_t coriolis_;

public:
	BarVortex(const Mesh & m, rp_t rp, coriolis_t coriolis, double tau, double sigma, double mu);

	/**
	 * F   -- �������� ������� �� ���������� ���� �� �������
	 * bnd -- ��������� �������
	 * t   -- �����
	 */
	void calc(double * Ans, const double * F, const double * bnd, double t);


	/**
	 * d L(psi)/dt + J(psi, L(z)) + J(z, L(psi)) + J(psi, l + h) + sigma L(psi) - mu LL(psi) = 0
	 * L = Laplace
	 */
	void calc_L(double * Ans, const double * F, const double * z, const double * bnd, double t);
	void calc_L_1(double * Ans, const double * F, const double * z, const double * bnd, double t);

	void calc_LT(double * Ans, const double * F, const double * z, const double * bnd, double t);

	/* J(psi, L(z)) + J(z, L(psi)) + J(psi, l + h) + sigma L(psi) - mu LL(psi) */
	void L_spectr(double * u1, const double * u, const double * z, const double * bnd);
	void LT_spectr(double * u1, const double * u, const double * z, const double * bnd);

	void S_step(double * Ans, const double * F);
	void L_step(double * Ans, const double * F, const double * z);
	void L_1_step(double * Ans, const double * F, const double * z);
	void LT_step(double * Ans, const double * F, const double * z);
};

#endif /* BARVORTEX_H */

