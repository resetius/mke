#ifndef LAPL_H
#define LAPL_H
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

/**
 * обращает оператор лапласа на плоской области
 * @param Ans - ответ
 * @param m - сетка
 * @param F - правая часть
 * @param bnd - краевое условие
 */
void laplace_solve(double * Ans, const Mesh & m, 
				   const double * F, const double * bnd);

void sphere_laplace_solve(double * Ans, const Mesh & m, 
						  const double * F, const double * bnd);

class SphereChafe {
public:
	struct integrate_cb_data
	{
		double tau;
		double mu;
		double sigma;
	};

private:
	const Mesh & m_;
	Matrix laplace_; /* Лапласиан */
	Matrix A_;       /* Матрица левой части */
	double tau_;
	double mu_;
	double sigma_;

	integrate_cb_data data1_;

public:
	SphereChafe(const Mesh & m, double tau, double sigma, double mu);
	~SphereChafe() {}

	void solve(double * Ans, const double * X0,
						const double * bnd);
};

class Chafe {
public:
	struct integrate_cb_data
	{
		double tau;
		double mu;
		double sigma;
	};

private:
	const Mesh & m_;
	Matrix laplace_; /* Лапласиан */
	Matrix A_;       /* Матрица левой части */
	double tau_;
	double mu_;
	double sigma_;

	integrate_cb_data data1_;

public:
	Chafe(const Mesh & m, double tau, double sigma, double mu);
	~Chafe() {}

	void solve(double * Ans, const double * X0,
						const double * bnd);
};

#endif /* LAPL_H */
