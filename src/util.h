#ifndef UTIL_H
#define UTIL_H
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

#ifdef __cplusplus
extern "C" {
#endif


#ifdef WIN32
#define inline __inline
#endif

#if 1
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
#endif

#if 0
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
#endif

int gauss (double *A, double *b, double *x, int n);

//интеграл от x^k*y^n по трапеции
//y=k1x+b1 y=k2x+b2 
//прямые ограничивающие трапецию
//x меняется от x1 до x3 
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

void matrix_print(const double * A, int n);
void vector_print(const double * A, int n);

void matrix_mult_vector(double * r, const double * A, const double * x, int n);

struct Sparse {
	int * Ap;
	int * Ai;
	double * Ax;
};

/**
 * Если матрица хранится по столбцам, то умножаем на вектор СЛЕВА
 * r = x A
 * Если матрица хранится по строкам, то умножаем на вектор СПРАВА
 * r = A x
 */
void sparse_mult_vector_l(double * r, const struct Sparse * A, const double * x, int n);

/**
 * Если матрица хранится по столбцам, то умножаем на вектор СПРАВА
 * r = A x
 * Если матрица хранится по строкам, то умножаем на вектор СЛЕВА
 * r = x A
 */
void sparse_mult_vector_r(double * r, const struct Sparse * A, const double * x, int n);

void sparse_print(const struct Sparse * A, int n, FILE * f);

/**
 * r = k1 * a + k2 * b
 */
void vector_sum1(double * r, const double * a, const double *b, double k1, double k2, int n);

/**
 * a = b * k
 */
void vector_mult_scalar(double * a, const double * b, double k, int n);

/**
 * r = a + b
 */
void vector_sum(double * r, const double * a, const double *b, int n);

/**
 * поэлементное умножение
 * r = a * b
 */
void vector_mult(double * r, const double * a, const double * b, int n);

void vector_diff(double * r, const double * a, const double * b, int n);

double get_full_time();

typedef double (*fx_t)(double x, void * data);
double gauss_kronrod15(double a, double b, fx_t fm, void * data);

void set_fpe_except();

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
class Timer {
	double t1_;

public:
	Timer(): t1_(get_full_time()) {}
	double elapsed() { return (get_full_time() - t1_) / 100.0; }
	void restart() { t1_ = get_full_time(); }
};
#endif

#endif /* UTIL_H */

