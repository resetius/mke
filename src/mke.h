#ifndef MKE_H
#define MKE_H
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

#include <stdio.h>

#include "polynom.h"

typedef unsigned int uint;

struct Point {
	double x;
	double y;

	Point(): x(0), y(0) {}
	Point(double x1, double y1): x(x1), y(y1) {}
	Point(double *x1): x(x1[0]), y(x1[1]) {}

	Point operator / (double k)
	{
		return Point(x / k, y / k);
	}

	Point operator * (double k)
	{
		return Point(x * k, y * k);
	}
};

inline Point operator + (const Point & p1, const Point & p2)
{
	return Point(p1.x + p2.x, p1.y + p2.y);
}

/** 
 * Точка на сетке может входить одновременно в несколько областей, 
 * для поддержки этой возможности нужен этот класс
 */
struct MeshPoint {
	std::vector < Point > p;

	MeshPoint() {}
	
	MeshPoint(double x1, double y1) {
		add(Point(x1, y1));
	}

	MeshPoint(double *x1) {
		add(Point(x1));
	}

	void add(const Point & p1) {
		p.push_back(p1);
	}

	double x(int zone = 0) const {
		return p[zone].x;
	}

	double y(int zone = 0) const {
		return p[zone].y;
	}
};

struct Triangle {
	int p[3];  /* point numbers */
	int z;     /* zone number   */
	std::vector < Polynom > phik;
	double x[3];
	double y[3];

	Triangle(int p1_, int p2_, int p3_, int zone = 0)
	{
		p[0] = p1_;
		p[1] = p2_;
		p[2] = p3_;
		z    = zone;
	}

	double X(int i, const std::vector < MeshPoint > & ps) const {
		return ps[p[i]].x(z);
	}

	double Y(int i, const std::vector < MeshPoint > & ps) const {
		return ps[p[i]].y(z);
	}

	void prepare(const std::vector < MeshPoint > & ps)
	{
		std::vector < Polynom > & r = phik;
		r.reserve(3);

		// p0
		r.push_back((P2X - X(1, ps)) * (Y(2, ps) - Y(1, ps)) 
			- (P2Y - Y(1, ps)) * (X(2, ps) - X(1, ps)));
		// p1
		r.push_back((P2X - X(0, ps)) * (Y(2, ps) - Y(0, ps)) 
			- (P2Y - Y(0, ps)) * (X(2, ps) - X(0, ps)));
		// p2
		r.push_back((P2X - X(0, ps)) * (Y(1, ps) - Y(0, ps)) 
			- (P2Y - Y(0, ps)) * (X(1, ps) - X(0, ps)));

		for (uint i = 0; i < 3; ++i)
		{
			r[i] /= r[i].apply(X(i, ps), Y(i, ps));
		}

		x[0] = X(0, ps); y[0] = Y(0, ps);
		x[1] = X(1, ps); y[1] = Y(1, ps);
		x[2] = X(2, ps); y[2] = Y(2, ps);
	}

	const std::vector < Polynom > & elem1() const
	{
		return phik;
	}

	const Polynom & elem1(int p1) const
	{
		if (p1 == p[0]) {
			return phik[0];
		} else if (p1 == p[1]) {
			return phik[1];
		} else {
			return phik[2];
		}
	}
};

struct Mesh {
	typedef std::vector < Triangle > triangles_t;
	typedef std::vector < MeshPoint > points_t; 
	typedef std::vector < int > points_flags_t;
	triangles_t tr;
	points_t ps;
	//флаги
	//0 - внутренняя точка
	//1 - точка на границе области
	points_flags_t ps_flags;

	// точка -> треугольники в точке
	std::vector < std::vector < int > > adj;
	// номера внутренних точек
	std::vector < int > inner;
	// номера внешних точек
	std::vector < int > outer;
	// соответствие номера в массиве ps номеру в массиве inner или outer
	std::vector < int > p2io;

	bool load(FILE * f);
	void prepare();
	void info();
};

/**
 * to global coordinates
 */
typedef double (* x_t)(double u, double v);

void print_function(FILE * to, double * ans, const Mesh & m, 
					x_t x = 0, x_t y = 0, x_t z = 0);

void print_inner_function(FILE * to, double * ans, const Mesh & m, 
					x_t x = 0, x_t y = 0, x_t z = 0);

class Matrix;

/**
 * Callback. Вызывается для всех функций phi_i, phi_j, определенных
 * в точке point на треугольнике tr
 */
typedef double (* right_part_cb_t)
(
	const Polynom  & phi_i, 
	const Polynom  & phi_j, 
	const Triangle & tr,  /* треугольник */
	const Mesh & mesh,    /* сетка       */
	int point_i,          /* номер узла  */
	int point_j,
	void * user_data /* сюда могу входить любые данные, 
	                    например значение F на пред шаге*/
);

/**
 * Создает матрицу системы.
 * Вызывает right_part_cb_t для всех функций phi_i, phi_j, определенных
 * в общей точке point на треугольнике tr
 */
void generate_right_part(double * b, const Mesh & m, right_part_cb_t right_part_cb, void * user_data);

void generate_full_right_part(double * b, const Mesh & m, right_part_cb_t right_part_cb, void * user_data);

/**
 * Callback. Вызывается для всех функций phi_i, phi_j, определенных
 * в точке point на треугольнике tr
 */
typedef double (* integrate_cb_t)
(
	const Polynom & phi_i, 
	const Polynom & phi_j, 
	const Triangle & tr, /* номер треугольника */
	const Mesh & mesh,   /* сетка */
	int point_i,         /* номер точки */
	int point_j,
	void * user_data     /* сюда могу входить любые данные */
);

/**
 * Создает матрицу системы.
 * Вызывает integrate_cb для всех функций phi_i, phi_j, определенных
 * в общей точке point на треугольнике tr
 */
void generate_matrix(Matrix & A, const Mesh & m, integrate_cb_t integrate_cb, void * user_data);

void generate_full_matrix(Matrix & A, const Mesh & m, integrate_cb_t integrate_cb, void * user_data);

/**
 * генерирует матрицу для интеграции краевых условий в правую часть
 * inner.size() x outer.size()
 * в cb передается phi_j где j точка границы
 * phi_i, где i внутренняя точка
 */
void generate_boundary_matrix(Matrix & A, const Mesh & m, right_part_cb_t right_part_cb, void * user_data);

/**
 * Решает систему
 */
void mke_solve(double * answer, const double * bnd, double * rp, Matrix & A, const Mesh & m);

/**
 * Решает систему.
 * answer - это только значения во внутренних точках!
 */
void mke_solve2(double * answer, double * rp, Matrix & A, const Mesh & m);

/* добавляем краевые условия */
void mke_p2u(double * p, const double * u, const double * bnd, const Mesh & m);

/* убираем краевые условия */
void mke_u2p(double * u, const double * p, const Mesh & m);

/** 
 * в следующих функциях вызывается 
 * scalar_cb_t
 * для вычисления скалярного произведения двух базисных функций 
 * на треугольнике
 */
typedef double (* scalar_cb_t)
(
	const Polynom & phi_i, 
	const Polynom & phi_j, 
	const Triangle & tr, /* номер треугольника */
	const Mesh & mesh,   /* сетка */
	int point_i,
	int point_j,
	void * user_data     /* сюда могу входить любые данные */
);

/* тут вычисляется интеграл от произведения функций по треугольнику */
double generic_scalar_cb(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, const Mesh & m, int, int, void * );

/* тут вычисляется интеграл от произведения функций по треугольнику на сфере */
double sphere_scalar_cb(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, const Mesh & m, int, int, void * user_data);

void convolution(double * ans, const double * u, const double * v, const Mesh & m, scalar_cb_t cb, void * user_data);

/* сеточное скалярное произведение двух функций */
double mke_scalar(const double * u, const double * v, 
				  const Mesh & m, scalar_cb_t cb = generic_scalar_cb, void * user_data = 0);
/* сеточная норма */
double mke_norm(const double * u, const Mesh & m, scalar_cb_t cb = generic_scalar_cb, void * user_data = 0);
/* сеточное расстояние */
double mke_dist(const double * u, const double * v, const Mesh & m, scalar_cb_t cb = generic_scalar_cb, void * user_data = 0);

typedef double (* f_xy_t)(double x, double y);
typedef double (* f_xyt_t)(double x, double y, double t);

/* проектирование непрерывной функции f(x,y) на сетку */
void mke_proj(double * F, const Mesh & mesh, f_xy_t f);

/* проектирование непрерывной функции f(x,y) на границу сетки */
void mke_proj_bnd(double * F, const Mesh & m, f_xy_t f);

/* проектирование функции F1 на границу сетки */
void mke_proj_bnd(double * F, const double * F1, const Mesh & m);

/* проектирование непрерывной функции f(x,y,t) на сетку */
void mke_proj(double * F, const Mesh & mesh, f_xyt_t f, double t);

/* проектирование непрерывной функции f(x,y,t) на границу сетки */
void mke_proj_bnd(double * F, const Mesh & m, f_xyt_t f, double t);

#endif /* MKE_H */

