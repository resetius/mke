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
 * Сферический Лаплас
 */
class SphereLaplace {
	Matrix idt_;
	Matrix laplace_;
	const Mesh & m_;

public:
	SphereLaplace(const Mesh & m);

	/**
	 * Находит оператор Лапласа функции F во внутренних
	 * точках. В точках границы просто кладет значение из bnd.
	 */
	void calc1(double * Ans, const double * F, const double * bnd);

	/**
	 * Находит оператор Лапласа функции F во внутренних точках.
	 * Возвращает вектор, содержащий ТОЛЬКО внутренние точки
	 */
	void calc2(double * Ans, const double * F);

/**
 * обращает оператор лапласа на плоской области
 * @param Ans - ответ
 * @param m - сетка
 * @param F - правая часть
 * @param bnd - краевое условие
 */
	void solve(double * Ans, const double * F, const double * bnd);
};

class SphereChafe {
private:
	const Mesh & m_;
	SphereLaplace laplace_; /* Лапласиан */
	Matrix  A_;             /* Матрица левой части */
	double tau_;
	double mu_;
	double sigma_;

	static double 
	schafe_integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point,
                     SphereChafe * d);

	struct schafe_right_part_cb_data;

	static double 
	schafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point,
                      schafe_right_part_cb_data * d);

public:
	SphereChafe(const Mesh & m, double tau, double sigma, double mu);
	~SphereChafe() {}

	void solve(double * Ans, const double * X0,
						const double * bnd);
};

/**
 * Лаплас на плоской сетке
 */
class Laplace {
	Matrix idt_;
	Matrix laplace_;
	const Mesh & m_;

public:
	Laplace(const Mesh & m);

	/**
	 * Находит оператор Лапласа функции F во внутренних
	 * точках. В точках границы просто кладет значение из bnd.
	 */
	void calc1(double * Ans, const double * F, const double * bnd);

	/**
	 * Находит оператор Лапласа функции F во внутренних точках.
	 * Возвращает вектор, содержащий ТОЛЬКО внутренние точки
	 */
	void calc2(double * Ans, const double * F);

/**
 * обращает оператор лапласа на плоской области
 * @param Ans - ответ
 * @param m - сетка
 * @param F - правая часть
 * @param bnd - краевое условие
 */
	void solve(double * Ans, const double * F, const double * bnd);
};

class Chafe {
private:
	const Mesh & m_;
	Laplace laplace_; /* Лапласиан */
	Matrix A_;        /* Матрица левой части */
	double tau_;
	double mu_;
	double sigma_;

	static double chafe_integrate_cb(const Polynom & phi_i, 
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m, int point,
		const Chafe * d);

	struct chafe_right_part_cb_data;

	static double 
	chafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point,
                      chafe_right_part_cb_data * d);

public:
	Chafe(const Mesh & m, double tau, double sigma, double mu);
	~Chafe() {}

	void solve(double * Ans, const double * X0,
						const double * bnd, double t);
};

#endif /* LAPL_H */

