#ifndef LAPL_H
#define LAPL_H
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

#include "phelm.h"
#include "solver.h"

using phelm::Matrix;
using phelm::Mesh;
using phelm::Triangle;
using phelm::Polynom;

/**
 * @defgroup aux Решатели некоторых задач математической физики.
 * @{
 */

/**
 * The fast sphere norm calculator.
 */
class SphereNorm {
public:
	const Mesh & m_; ///< mesh
	Matrix NORM_;    ///< for fast norm calculation

	/**
	 * constructor.
	 * @param m - mesh
	 */
	SphereNorm(const Mesh & m): m_(m), NORM_((int)m_.ps.size()) 
	{
		phelm::generate_full_matrix(NORM_, m_, phelm::sphere_scalar_cb, (void*)0);
	}

	/**
	 * distance between two vectors.
	 * @param u - input vector
	 * @param v - input vector
	 * @return distance between u and v
	 */
	double dist(const double * u, const double * v)
	{
		return phelm::fast_dist(u, v, m_, NORM_);
	}

	/**
	 * norm of vector
	 * @param u - input vector
	 * @return norm of v
	 */
	double norm(const double * u)
	{
		return phelm::fast_norm(u, m_, NORM_);
	}

	/**
	 * inner product of two vectors
	 * @param u - input vector
	 * @param v - input vector
	 * @return inner product of u and v
	 */
	double scalar(const double * u, const double * v)
	{
		return phelm::fast_scalar(u, v, m_, NORM_);
	}
};

/**
 * Сферический Лаплас.
 \f[
  \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi
 \f]
 */
class SphereLaplace {
public:
	Matrix idt_;
	Matrix laplace_;
	Matrix bnd1_; // L^-1
	Matrix bnd2_; // L^-1
	Matrix bnd3_; // L
	const Mesh & m_;

public:
	/**
	 * constructor.
	 * @param m - mesh
	 */
	SphereLaplace(const Mesh & m);

	/**
	 * Находит оператор Лапласа функции F во внутренних
	 * точках. В точках границы просто кладет значение из bnd.
	 * @param Ans - ответ
	 * @param F - функция лаплас от которой хотим найти
	 * @param bnd - значение на границе, которое хотим получить
	 */
	void calc1(double * Ans, const double * F, const double * bnd);

	/**
	 * Находит оператор Лапласа функции F во внутренних точках.
	 * Возвращает вектор, содержащий ТОЛЬКО внутренние точки
	 * @param Ans - ответ
	 * @param F - функция лаплас от которой хотим найти
	 */
	void calc2(double * Ans, const double * F);

/**
 * обращает оператор лапласа на плоской области
 * @param Ans - ответ
 * @param F - правая часть
 * @param bnd - краевое условие
 */
	void solve(double * Ans, const double * F, const double * bnd);
};

/**
 * Считает один элемент конечно-элементной матрицы для
 * оператора лапласа на сфере.
 * @param phi_i - базисная функция
 * @param phi_j - базисная функция
 * @param trk - треугольник по которому считаем интеграл
 * @param ps - точки сетки
 */
double slaplace(const Polynom & phi_i, const Polynom & phi_j, 
		const Triangle & trk, const Mesh::points_t & ps);

/**
 * Считает один элемент конечно-элементной матрицы для
 * оператора лапласа на плоскости.
 * @param phi_i - базисная функция
 * @param phi_j - базисная функция
 * @param trk - треугольник по которому считаем интеграл
 * @param ps - точки сетки
 */
double laplace(const Polynom & phi_i, const Polynom & phi_j,
		const Triangle & trk, const Mesh::points_t & ps);

/**
 * Chafe-Infante equation on sphere.
 * \f$ \frac{du}{dt} = \mu \Delta u - \sigma u + f (u)\f$, where
 * \f$Delta\f$ is SphereLaplace. @see SphereLaplace
 */
class SphereChafe {
private:
	const Mesh & m_;
	SphereLaplace laplace_; /* Лапласиан */
	Matrix  A_;             /* Матрица левой части */

public:
	double tau_;   ///<time step
	double mu_;    ///<\f$\mu\f$
	double sigma_; ///<\f$\sigma\f$

public:
	/**
	 * Constructor.
	 * @param m - mesh
	 * @param tau - time step
	 * @param sigma - \f$\sigma\f$
	 * @param mu - \f$\mu\f$
	 */
	SphereChafe(const Mesh & m, double tau, double sigma, double mu);
	~SphereChafe() {}

/**
 * Solve Chafe-Infante equation on sphere.
 * \f$\frac{du}{dt} = \mu \Delta u - \sigma u + f (u)\f$
 * @param Ans - output vector
 * @param X0 - intput vector (previous time step)
 * @param bnd - boundary condition
 * @param t - time
 */
	void solve(double * Ans, const double * X0,
						const double * bnd, double t);
};

/**
 * Лаплас на плоской сетке
 */
class Laplace {
	Matrix idt_;      // inner
	Matrix laplace_;  // inner

//	Matrix fidt_;     // full
//	Matrix flaplace_; // full

	Matrix bnd1_; // L^-1
	Matrix bnd2_; // L^-1
	Matrix bnd3_; // L
	const Mesh & m_;

	void init_boundary_vectors();

public:
	/**
	 * constructor.
	 * @param m - mesh
	 */
	Laplace(const Mesh & m);

	/**
	 * Находит оператор Лапласа функции F во внутренних
	 * точках. В точках границы просто кладет значение из bnd.
	 * @param Ans - ответ
	 * @param F - функция лаплас от которой хотим найти
	 * @param bnd - значение награнице, которое хотим получить
	 */
	void calc1(double * Ans, const double * F, const double * bnd);

	/**
	 * Находит оператор Лапласа функции F во внутренних точках.
	 * Возвращает вектор, содержащий ТОЛЬКО внутренние точки
	 * @param Ans - ответ
	 * @param F - функция лаплас от которой хотим найти
	 */
	void calc2(double * Ans, const double * F);

	/**
	 * обращает оператор лапласа на плоской области
	 * @param Ans - ответ
	 * @param F - правая часть
	 * @param bnd - краевое условие
	 */
	void solve(double * Ans, const double * F, const double * bnd);
};

/**
 * 2D Chafe-Infante equation.
 * \f$ \frac{du}{dt} = \mu \Delta u - \sigma u + f (u)\f$, where
 * \f$Delta\f$ is Laplace. @see Laplace
 */
class Chafe {
private:
	const Mesh & m_;
	Laplace laplace_; /* Лапласиан */
	Matrix A_;        /* Матрица левой части */

public:
	double tau_;   ///<time step
	double mu_;    ///<\f$\mu\f$
	double sigma_; ///<\f$\sigma\f$

public:
	/**
	 * constructor.
	 * @param m - mesh
	 * @param tau - time step
	 * @param sigma - \f$\sigma\f$
	 * @param mu - \f$\mu\f$
	 */
	Chafe(const Mesh & m, double tau, double sigma, double mu);
	~Chafe() {}

	/**
	 * Solve 2D Chafe-Infante equation.
	 * \f$\frac{du}{dt} = \mu \Delta u - \sigma u + f (u)\f$
	 * @param Ans - output vector
	 * @param X0 - intput vector (previous time step)
	 * @param bnd - boundary condition
	 * @param t - time
	 */
	void solve(double * Ans, const double * X0,
						const double * bnd, double t);
};

/** @} */

#endif /* LAPL_H */

