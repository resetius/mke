#ifndef MKE_H
#define MKE_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/**
 * @file 
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @page License
 * @section LICENSE
 *
 * @verbatim
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
 * @endverbatim
 *
 * @mainpage Phelm Documentation
 * @section into_sec Introduction
 * Эта библиотека предназначена для решения уравнений в частных производных
 * Методом Конечных Элементов.
 * Библиотека поддерживает решение уравнений на двумерных гладких 
 * многообразиях в локальных координатах. 
 *
 * Для решения задачи от пользователя требуются: 
 * - разбиение многообразия на области
 * - задание в каждой области локальной системы координат
 * - функция вычисления поверхностного интеграла по треугольнику в локальных координатах.
 *
 * Библиотека уже содержит построитель сферических поверхностей и плоских прямоугольных
 * областей.
 * 
 * Пример работы построителя сетки на сфере. 
 * В данном случае сфера разбита на четыре области.
 *
 * @image html zones.png
 * 
 * @section Examples
 * Доступны следующие примеры:
    -# @ref test_laplace.cpp "Laplace on a 2D-domain"
  \f{eqnarray*}
  \Delta u &=& f(x, y) \\
  u|_{\partial\Omega}&=&u_0
  \f}
    -# @ref test_laplace.cpp "Laplace on a sphere"
    -# @ref test_system_laplace.cpp "Double Laplace"
  \f{eqnarray*}
  \Delta u + v &=& f(x, y)\\
  u + \Delta v &=& g(x, y)\\
  u|_{\partial\Omega}&=&u_0\\
  v|_{\partial\Omega}&=&v_0\\
  \f}
    -# @ref test_chafe.cpp "Chafe-Infante equation on a 2D-Domain"
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{\partial\Omega}&=&a \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
    -# @ref test_chafe.cpp "Chafe-Infante equation on a sphere"
    -# @ref test_barvortex.cpp "BarVortex equation on sphere"
  \f{eqnarray*}
  \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi &=& f(\varphi, \lambda) \\
	\psi|_{t=0}=\psi_0
  \f}
    -# @ref test_baroclin.cpp  "Baroclin equations on sphere"
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    - \alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,\\
	u_1|{t=0}&=&u_{10}\\
	u_2|{t=0}&=&u_{20}\\
 \f}
 *
 * @page Build
 * @section build_sec Build
 * @subsection Unix-like
 * @verbatim
  mkdir build-directory
  cd build-directory
  cmake -DCMAKE_BUILD_TYPE=Debug path-to-sources   # for Debug build
  cmake -DCMAKE_BUILD_TYPE=Release path-to-sources # for Release build
  make
  @endverbatim
 * @subsection Windows
 * @verbatim
  mkdir build-directory
  cd build-directory
  cmake -G "Visual Studio 2009" #place your version of Visual Studio here
  @endverbatim
 */

/**
 * 
  @example test_laplace.cpp
  Laplace on a 2D-Domain
  \f{eqnarray*}
  \Delta u &=& f(x, y) \\
  u|_{\partial\Omega}&=&u_0
  \f}
  @example test_slaplace.cpp
  Laplace on a sphere
  @example test_system_laplace.cpp
  Double Laplace on a 2D-Domain
  \f{eqnarray*}
  \Delta u + v &=& f(x, y)\\
  u + \Delta v &=& g(x, y)\\
  u|_{\partial\Omega}&=&u_0\\
  v|_{\partial\Omega}&=&v_0\\
  \f}
  @example test_chafe.cpp
  Chafe-Infante equation on a 2D-Domain
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{\partial\Omega}&=&a \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
  @example test_schafe.cpp
  Chafe-Infante equation on a sphere
  @example test_barvortex.cpp
  BarVortex equation on sphere
  \f[
  \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
  \f]
  @example test_baroclin.cpp
  Baroclin equations on sphere
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    - \alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
  */

#include <stdio.h>

#include "polynom.h"

typedef unsigned int uint;

/**
 * @defgroup main Классы и Функции для работы с сеткой и сеточными функциями.
 * @{
 */

/**
 * @namespace phelm
 *
 * That namespace contains all Phelm functions and classes.
 */
namespace phelm {

/**
 * @defgroup main Main Functions and Classes.
 * @{
 */


/**
 * Point class represents a 2-dimensional point on a plane.
 */
struct Point {
	double x; ///< x coordinate
	double y; ///< y coordinate

	/** default constructor */
	Point(): x(0), y(0) {}
	/** initialize point by x1 and y1 */
	Point(double x1, double y1): x(x1), y(y1) {}
	/** initialize point by x[2] array */
	Point(double *x1): x(x1[0]), y(x1[1]) {}

	/** divide each coordinate by k */
	Point operator / (double k)
	{
		return Point(x / k, y / k);
	}

	/** multiply each coordinate by k */
	Point operator * (double k)
	{
		return Point(x * k, y * k);
	}
};

/**
 * @relates Point
 * sum of two points.
 * @param p1 - input point
 * @param p2 - input point
 * @return p1 + p2
 */
inline Point operator + (const Point & p1, const Point & p2)
{
	return Point(p1.x + p2.x, p1.y + p2.y);
}

/**
 * Точка на сетке может входить одновременно в несколько областей, 
 * для поддержки этой возможности нужен этот класс
 */
struct MeshPoint {
	/**
	 * p[i] координаты точки в области i
	 */
	std::vector < Point > p;

	/**
	 * default constructor.
	 */
	MeshPoint() {}
	
	/**
	 * инициализация по двум координатам
	 * @param x - координата x
	 * @param y - координата y
	 */
	MeshPoint(double x, double y) {
		add(Point(x, y));
	}

	/**
	 * инициализация по массиву x[2]
	 * @param x - массив из двух координат
	 */
	MeshPoint(double *x) {
		add(Point(x));
	}

	/**
	 * добавляет координаты для следующей зоны
	 */
	void add(const Point & p1) {
		p.push_back(p1);
	}

	/**
	 * Возвращает x координату точки в области zone
	 * @param zone - номер области
	 * @return координата x
	 */
	double x(int zone = 0) const {
		return p[zone].x;
	}

	/**
	 * Возвращает y координату точки в области zone
	 * @param zone - номер области
	 * @return координата y
	 */
	double y(int zone = 0) const {
		return p[zone].y;
	}
};

/**
 * Triangle class.
 */
struct Triangle {
	int p[3];  ///< point number
	int z;     ///< zone number
	std::vector < Polynom > phik; ///< элементы на треугольнике
	double x[3]; ///< x координата вершин
	double y[3]; ///< y координата вершин

	/**
	 * инициализация треугольника по номерам точек и номеру области.
	 * @param p1 - точка 1
	 * @param p2 - точка 2
	 * @param p3 - точка 3
	 * @param zone - номер области, которой пренадлежит треугольник
	 */
	Triangle(int p1, int p2, int p3, int zone = 0)
	{
		p[0] = p1;
		p[1] = p2;
		p[2] = p3;
		z    = zone;
	}

	/**
	 * возвращает x координату вершины по номеру вершины.
	 * @param i - номер вершины (от 0 до 2)
	 * @param ps - массив точек сетки
	 * @return x координата
	 */
	double X(int i, const std::vector < MeshPoint > & ps) const 
	{
		return ps[p[i]].x(z);
	}

	/**
	 * возвращает y координату вершины по номеру вершины.
	 * @param i - номер вершины (от 0 до 2)
	 * @param ps - массив точек сетки
	 * @return y координата
	 */
	double Y(int i, const std::vector < MeshPoint > & ps) const 
	{
		return ps[p[i]].y(z);
	}

	/**
	 * инициализирует массивы x и y по массиву точек
	 * @param ps - массив точек
	 */
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

	/**
	 * возвращает конечные элементы первого порядка.
	 * @return массив конечных элементов первого порядка
	 */
	const std::vector < Polynom > & elem1() const
	{
		return phik;
	}

	/**
	 * возвращает конечный элемент по номеру точки.
	 * @param p1 - номер точки
	 * @return конечный элемент
	 */
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

/**
 * Mesh class.
 */
struct Mesh {
	typedef std::vector < Triangle > triangles_t;///<контейнер треугольников
	typedef std::vector < MeshPoint > points_t;  ///<контейнер точек
	typedef std::vector < int > points_flags_t;  ///<контейнер свойст точек

	triangles_t tr; ///<массив треугольников
	points_t ps;    ///<массив точек

	/**
	 * флаги.
	 * 0 - внутренняя точка
	 * 1 - точка на границе области
	 */
	points_flags_t ps_flags;

	//! точка -> треугольники в точке
	std::vector < std::vector < int > > adj;
	//! номера внутренних точек
	std::vector < int > inner;
	//! номера внешних точек
	std::vector < int > outer;
	//! соответствие номера в массиве ps номеру в массиве inner или outer
	std::vector < int > p2io;

	/**
	 * загружает сетку из файла.
	 * @param f - файл
	 * @return true в случае успеха
	 */
	bool load(FILE * f);

	/**
	 * подгатавливает треугольники.
	 */
	void prepare();

	/**
	 * распечатывает в stdout информацию о сетке.
	 */
	void info();
};

/** @} */ /* main */

/**
 * @defgroup print Функции печати 
 * @ingroup print
 * @{
 */

/**
 * callback that converts local coordinates to global coordinates.
 */
typedef double (* x_t)(double u, double v);

/**
 * Prints the function @param ans to file @param to
 *
 * @param to output file
 * @param ans function to output
 * @param m mesh
 * @param x (optional) local coordinates to global 'x' converter
 * @param y (optional) local coordinates to global 'y' converter
 * @param z (optional) local coordinates to global 'z' converter
 */
void print_function(FILE * to, double * ans, const Mesh & m, 
					x_t x = 0, x_t y = 0, x_t z = 0);

/**
 * Prints the function @param ans to the file fname
 *
 * @param fname output file
 * @param ans function to output
 * @param m mesh
 * @param x (optional) local coordinates to global 'x' converter
 * @param y (optional) local coordinates to global 'y' converter
 * @param z (optional) local coordinates to global 'z' converter
 */
void print_function(const char * fname, double * ans, const Mesh & m, 
					x_t x = 0, x_t y = 0, x_t z = 0);

/**
 * Prints the inner part of the function ans to file to
 *
 * @param to output file
 * @param ans function to output
 * @param m mesh
 * @param x (optional) local coordinates to global 'x' converter
 * @param y (optional) local coordinates to global 'y' converter
 * @param z (optional) local coordinates to global 'z' converter
 */
void print_inner_function(FILE * to, double * ans, const Mesh & m, 
					x_t x = 0, x_t y = 0, x_t z = 0);

/**
 * Prints the inner part of the function ans to file to
 *
 * @param to output file name
 * @param ans function to output
 * @param m mesh
 * @param x (optional) local coordinates to global 'x' converter
 * @param y (optional) local coordinates to global 'y' converter
 * @param z (optional) local coordinates to global 'z' converter
 */
void print_inner_function(const char * to, double * ans, const Mesh & m,
		x_t x = 0, x_t y = 0, x_t z = 0);

/** @} */ /* print */

/**
 * @defgroup scalar Скалярные произведения, нормы, расстояния.
 * @ingroup scalar
 * @{
 */

/**
 * тут вычисляется интеграл от произведения функций по треугольнику.
 */
double generic_scalar_cb(const Polynom & phi_i, const Polynom & phi_j,
						 const Triangle & trk, const Mesh & m, int, int,
						 int, int, void * );

/**
 * тут вычисляется интеграл от произведения функций по треугольнику на сфере.
 */
double sphere_scalar_cb(const Polynom & phi_i, const Polynom & phi_j,
						const Triangle & trk, const Mesh & m,
						int, int, int, int, void * user_data);

class Matrix;

/**
 * fast_scalar.
 */
double fast_scalar(const double * u, const double * v,
				   const Mesh & m, Matrix & mat);

/**
 * fast norm.
 */
double fast_norm(const double * u, const Mesh & m, Matrix & mat);

/**
 * fast_dist.
 */
double fast_dist(const double * u, const double * v,
				 const Mesh & m, Matrix & mat);

/** @} */ /* scalar */

/**
 * @defgroup proj Конвертеры и проекторы функций.
 * @ingroup main
 * @{
 */

/**
 * функция от двух пространственных переменных.
 */
typedef double (* f_xy_t)(double x, double y);

/**
 * функция от двух пространственных переменных и времени.
 */
typedef double (* f_xyt_t)(double x, double y, double t);

/**
 * добавляем краевые условия
 */
void p2u(double * p, const double * u, const double * bnd, const Mesh & m);

/**
 * убираем краевые условия
 */
void u2p(double * u, const double * p, const Mesh & m);

/**
 * проектирование непрерывной функции f(x,y) на сетку.
 */
void proj(double * F, const Mesh & mesh, f_xy_t f);

/**
 * проектирование непрерывной функции f(x,y) на границу сетки.
 */
void proj_bnd(double * F, const Mesh & m, f_xy_t f);

/**
 * проектирование функции F1 на границу сетки.
 */
void proj_bnd(double * F, const double * F1, const Mesh & m);

/**
 * инициализирует граничное значение.
 */
void set_bnd(double * F, const double * bnd, const Mesh & m);

/**
 * проектирование непрерывной функции f(x,y,t) на сетку.
 */
void proj(double * F, const Mesh & mesh, f_xyt_t f, double t);

/**
 * проектирование непрерывной функции f(x,y,t) на границу сетки.
 */
void proj_bnd(double * F, const Mesh & m, f_xyt_t f, double t);

/** @} */ /* proj */

}

#include "phelm_generators.h"

#endif /* MKE_H */

