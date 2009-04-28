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

class BarVortex {
public:
	typedef double (*rp_t ) (double phi, double lambda, double mu, double sigma);
	typedef double (*coriolis_t) (double phi, double lambda);

private:
	const Mesh & m_;
	SphereLaplace l_;
	Jacobian j_;
	Matrix A_;

	std::vector < double > lh_; // l + h
	std::vector < double > f_;  // f right part

	double tau_;
	double sigma_;
	double mu_;
	rp_t rp_;
	coriolis_t coriolis_;

	static double integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
                     BarVortex * d);

	struct right_part_cb_data;
	static double 
	right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
                      right_part_cb_data * d);

public:
	BarVortex(const Mesh & m, rp_t rp, coriolis_t coriolis, double tau, double sigma, double mu);

	/**
	 * F   -- �������� ������� �� ���������� ���� �� �������
	 * bnd -- ��������� �������
	 * t   -- �����
	 */
	void calc(double * Ans, const double * F, const double * bnd, double t);
};

#endif /* BARVORTEX_H */

