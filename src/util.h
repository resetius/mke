#ifndef UTIL_H
#define UTIL_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/**
 * @file util.h
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @section LICENSE
 *
 * <pre>
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
 * </pre>
 *
 * @section DESCRIPTION
 * Misc functions.
 */

#ifdef __cplusplus
extern "C" {
#endif


#ifdef WIN32
#define inline __inline
#endif

	/**
	 * @defgroup misc Misc Functions and Classes.
	 * @{
	 */
/**
 * Power function.
 * @return the value of x raised to the power of p
 */
inline double
ipow(double x, int p)
{
	int i;
	double r = 1;
	for (i = 0; i < p; i++) {
		r *= x;
	}
	return r;
}

/*
inline double
ipow(double x, int n)
{
	double p, y;
	y = 1.0;
	p = x;

	while (1) {
		if (n & 1) {
			y = p * y;
		}
		n = n >> 1;
		if (n == 0) {
			return y;
		}
		p = p * p;
	}
}
*/

/**
 * The Gauss method. Solves a linear system Ax=b.
 *
 * @param A - the matrix of the system.
 * @param b - the right part.
 * @param x - the answer.
 * @param n - dimension.
 * @retuen 0 on success.
 */
int gauss (double *A, double *b, double *x, int n);

/**
 * интеграл от x^k*y^n по трапеции
 * y=k1x+b1 y=k2x+b2 
 * прямые ограничивающие трапецию
 * x меняется от x1 до x3 
 */
double trapezoid_integral(int k, int n, 
			 double k1, double b1, 
			 double k2, double b2, 
			 double x1, double x3);

/**
 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n cos x) dx dy
 */
double trapezoid_integral_cos(int k, int n,
	double k1, double b1,
	double k2, double b2,
	double x1, double x3);

/**
 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n sin x) dx dy
 */
double trapezoid_integral_sin(int k, int n,
	double k1, double b1,
	double k2, double b2,
	double x1, double x3);

/**
 * int [x1,x3] int [k1x+b1,k2x+b2] (x^k y^n/cos x) dx dy
 */
double trapezoid_integral_1_cos(int k, int n,
	double k1, double b1,
	double k2, double b2,
	double x1, double x3);

/**
 * Prints NxN matrix to stdout.
 * @param A - matrix
 * @param n - dimension
 */
void mat_print(const double * A, int n);

/**
 * Prints vector to stdout.
 * @param A - vector
 * @param n - dimension
 */
void vec_print(const double * A, int n);

/**
 * Product of NxN matrix and vector
 * @param r - output vector, r = Ax
 * @param A - intput matrix
 * @param x - input vector
 * @param n - dimension of matrix and vector
 */
void mat_mult_vector(double * r, const double * A, const double * x, int n);

/**
 * Sparse Matrix structure.
 * UMFPACK format.
 * Хранение по столбцам или по строкам.
 * Если хранится по строкам, то индексы в Ap это число элементов в строках, 
 * Ai индексы элементов в строках.
 * Если хранится по стобцам, то индексы в Ap это число элементов в столбцах,
 * Ai индексы элементов в столбцах.
 */
struct Sparse {
	int * Ap;    ///< the number of elements in columns or rows
	int * Ai;    ///< column or row indices
	double * Ax; ///< matrix values
};

/**
 * Если матрица хранится по столбцам, то умножаем на вектор СЛЕВА
 * r = x A
 * Если матрица хранится по строкам, то умножаем на вектор СПРАВА
 * r = A x
 */
void sparse_mult_vector_l(double * r, const struct Sparse * A,
						  const double * x, int n);

/**
 * Если матрица хранится по столбцам, то умножаем на вектор СПРАВА
 * r = A x
 * Если матрица хранится по строкам, то умножаем на вектор СЛЕВА
 * r = x A
 */
void sparse_mult_vector_r(double * r, const struct Sparse * A,
						  const double * x, int n);

/**
 * Prints sparse matrix to file.
 * @param A - input sparse matrix
 * @param n - dimension of sparse matrix
 * @param f - output file
 */
void sparse_print(const struct Sparse * A, int n, FILE * f);

/**
 * Linear combination of two vectors.
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param k1 - coefficient
 * @param k2 - coefficient
 * @param n - dimension of vectors
 * @return r = k1 * a + k2 * b
 */
void vec_sum1(double * r, const double * a,
			  const double *b, double k1, double k2, int n);

/**
 * product of vector by number
 * @param a - output vector
 * @param b - input vector
 * @param k - input number
 * @param n - dimension of vector
 * @return a = b * k
 */
void vec_mult_scalar(double * a, const double * b, double k, int n);

/**
 * sum of two vectors
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 * @return r = a + b
 */
void vec_sum(double * r, const double * a, const double *b, int n);

/**
 * Vector norm.
  \f[
  \sqrt{\sum_{i=0}v_i^2}
  \f]
 * @param v - input vector
 * @param n - dimension of vector
 * @return vector norm
 */
double vec_norm2(const double *v, int n);

/**
 * Inner product of two vectors.
   \f[
   \sum_{i=0}{n}a_i b_i
   \f]
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 * @return inner product of a and b
 */
double vec_scalar2(const double * a, const double * b, int n);

/**
 * поэлементное умножение
 * r = a * b
 */
void vec_mult(double * r, const double * a, const double * b, int n);

/**
 * difference of two vectors.
 * @param r - the output vector
 * @param a - the input vector
 * @param b - the input vector
 * @param n - the dimension of vectors
 * @return r = a - b
 */
void vec_diff(double * r, const double * a, const double * b, int n);

/**
 * Returs the number of seconds since epoch.
 * @return the number of seconds since epoch.
 */
double get_full_time();

/**
 * Function-callback that is passed is passed to gauss_kronrod15.
 */
typedef double (*fx_t)(double x, void * data);

/**
 * Gauss Kronrod quadrature formula.
 * http://en.wikipedia.org/wiki/Gauss-Kronrod_quadrature
 *
 * @param a
 * @param b
 * @param fm - function
 * @param data - the user data that are passed to the function
 * @return the integral of fm over the segment [a, b]
 */
double gauss_kronrod15(double a, double b, fx_t fm, void * data);

/**
 * Sets FPU exceptions.
 */
void set_fpe_except();

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/**
 * Timer class.
 */
class Timer {
	double t1_;

public:
	Timer(): t1_(get_full_time()) {}

	/**
	 * @return the number of seconds from timer initialize or restart. 
	 */
	double elapsed() { return (get_full_time() - t1_) / 100.0; }

	/**
	 * Restart timer.
	 */
	void restart() { t1_ = get_full_time(); }
};

/**
 * @}
 */

#endif

#endif /* UTIL_H */

