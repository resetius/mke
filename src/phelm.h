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
  Copyright (c) 2009 Alexey Ozeritsky
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  3. Redistributions in any form must be accompanied by information on
     how to obtain complete source code for the Phelm software and any
     accompanying software that uses the Phelm software.  The source code
     must either be included in the distribution or be available for no
     more than the cost of distribution plus a nominal fee, and must be
     freely redistributable under reasonable conditions.  For an
     executable file, complete source code means the source code for all
     modules it contains.  It does not include source code for modules or
     files that typically accompany the major components of the operating
     system on which the executable file runs.
 
  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  @endverbatim
 * If you want to use Phelm in closed-source software 
 * contact me by sending email to aozeritsky@gmail.com.
 *
 * @mainpage Phelm Documentation
 * @section into_sec Introduction
 * This library is meant to solve partial differential equations.
 * It implements the Finite Element Method.
 * The library is able to solve partial differential equations 
 * on two-dimensional smooth manifolds. 
 *
 * To solve a problem you need the following:
 * - triangulate the manifold;
 * - partition the manifold into subdomains;
 * - specify a local coordinate system in each subdomain;
 * - specify a function of surface integral in the local coordinates.
 *
 * The library already contains triangulation builders for spherical surfaces 
 * and flat rectangular domains.
 *
 * The sample result of mesh builder is shown below: 
 * @image html zones.png
 * In that case the sphere is built on 4 parts.
 * 
 * @section ex Usage examples
 * These examples are to give you some tips on Phelm features.
    -# @ref test_laplace.cpp "Laplace equation on a flat domain"
  \f{eqnarray*}
  \Delta u &=& f(x, y) \\
  u|_{\partial\Omega}&=&u_0
  \f}
    -# @ref test_laplace.cpp "Laplace equation on a sphere"
  \f{eqnarray*}
  \Delta \psi &=& f(\varphi, \lambda) \\
  \Delta \psi &=& \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi\\
  \psi|_{\partial\Omega}&=&\psi_0 \\
  \f}
    -# @ref test_system_laplace.cpp "Double Laplace equations on a flat domain"
  \f{eqnarray*}
  \Delta u + v &=& f(x, y)\\
  u + \Delta v &=& g(x, y)\\
  u|_{\partial\Omega}&=&u_0\\
  v|_{\partial\Omega}&=&v_0\\
  \f}
    -# @ref test_chafe.cpp "Chafe-Infante equation on a flat domain"
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{\partial\Omega}&=&a \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
    -# @ref test_chafe.cpp "Chafe-Infante equation on a sphere"
    -# @ref test_barvortex.cpp "The Barotropic vorticity equation on a sphere"
  \f{eqnarray*}
  \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi &=& f(\varphi, \lambda) \\
	\psi|_{t=0}=\psi_0
  \f}
    -# @ref test_baroclin.cpp  "The two-dimensional baroclinic atmosphere equations  on a sphere"
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	-\alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,\\
	u_1|_{t=0}&=&u_{10}\\
	u_2|_{t=0}&=&u_{20}\\
 \f}
 *
 * @page Build
 * @section build_sec Build
 * Phelm uses GMRES method to solve linear equations. If you do not like
 * that, install UMFPACK or place UMFPACK sources in 
 * /path-to-phelm/contrib/umfpack
 * 
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
  Laplace equation on a flat domain
  \f{eqnarray*}
  \Delta u &=& f(x, y) \\
  u|_{\partial\Omega}&=&u_0
  \f}
  @example test_slaplace.cpp
  Laplace equation on a sphere
  \f{eqnarray*}
  \Delta \psi &=& f(\varphi, \lambda) \\
  \Delta \psi &=& \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi\\
  \psi|_{\partial\Omega}&=&\psi_0 \\
  \f}
  @example test_system_laplace.cpp
  Double Laplace equations on a flat domain
  \f{eqnarray*}
  \Delta u + v &=& f(x, y)\\
  u + \Delta v &=& g(x, y)\\
  u|_{\partial\Omega}&=&u_0\\
  v|_{\partial\Omega}&=&v_0\\
  \f}
  @example test_chafe.cpp
  Chafe-Infante equation on a flat domain
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{\partial\Omega}&=&a \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
  @example test_schafe.cpp
  Chafe-Infante equation on a sphere
  @example test_barvortex.cpp
  the Barotropic vorticity equation
  \f[
  \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
  \f]
  @example test_baroclin.cpp
  The two-dimensional baroclinic atmosphere equations
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	-\alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
  */

#include <stdio.h>

#include "polynom.h"

typedef unsigned int uint;

#define PHELM_VERSION    1
#define PHELM_PATCHLEVEL 0
#define PHELM_SUBLEVEL   2

/**
 * @namespace phelm
 *
 * This namespace contains all Phelm functions and classes.
 */
namespace phelm {

/**
 * @defgroup main Mesh and mesh functions.
 * @{
 */


/**
 * Point class represents a 2-dimensional point on a plane.
 */
struct Point {
	double x; ///< x coordinate
	double y; ///< y coordinate

	/** Default constructor. */
	Point(): x(0), y(0) {}
	/**
	 * Initialization of x1 and y1.
	 * @param x1 - x coordinate
	 * @param y1 - y coordinate
	 */
	Point(double x1, double y1): x(x1), y(y1) {}
	/**
	 * Initialization of x[2] array.
	 * @param x1 - array x1[2]
	 */
	Point(double *x1): x(x1[0]), y(x1[1]) {}

	/**
	 * Divide each coordinate by k.
	 * @param k - a number
	 * @return new point
	 */
	Point operator / (double k)
	{
		return Point(x / k, y / k);
	}

	/**
	 * Multiply each coordinate by k.
	 * @param k - a number
	 * @return new point
	 */
	Point operator * (double k)
	{
		return Point(x * k, y * k);
	}
};

/**
 * Sum of two points.
 * @relates Point
 * @param p1 - input point
 * @param p2 - input point
 * @return p1 + p2
 */
inline Point operator + (const Point & p1, const Point & p2)
{
	return Point(p1.x + p2.x, p1.y + p2.y);
}

/**
 * A point of a manifold can be included into multiple subdomains.
 * That class supports this.
 */
struct MeshPoint {
	/**
	 * p[i] local coordinates of point in subdomain i.
	 * @todo what should we do if the point is included only in one subdomain
	 * with sequence number greater that 0 ?
	 */
	std::vector < Point > p;

	/**
	 * Default constructor.
	 */
	MeshPoint() {}
	
	/**
	 * Initialization of x and y.
	 * @param x - x coordinate
	 * @param y - y coordinate
	 */
	MeshPoint(double x, double y) {
		add(Point(x, y));
	}

	/**
	 * Initialization of array x[2].
	 * @param x - array x[2]
	 */
	MeshPoint(double *x) {
		add(Point(x));
	}

	/**
	 * Add local coordinates for next subdomain.
	 * @param p1 - point in local coordinates of next subdomain.
	 */
	void add(const Point & p1) {
		p.push_back(p1);
	}

	/**
	 * Returns x coordinate in subdomain.
	 * @param zone - sequence number of subdomain
	 * @return x coordiante
	 */
	double x(int zone = 0) const {
		return p[zone].x;
	}

	/**
	 * Returns y coordinate in subdomain.
	 * @param zone - a sequence number of subdomain
	 * @return y coordiante
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
	std::vector < Polynom > phik; ///< basis functions
	double x[3]; ///< vertices x coordinates
	double y[3]; ///< vertices y coordinates

	/**
	 * Initialization of point numbers and subdomain number.
	 * @param p1 - point 1
	 * @param p2 - point 2
	 * @param p3 - point 3
	 * @param zone - subdomain number, the triangle belongs to that subdomain
	 */
	Triangle(int p1, int p2, int p3, int zone = 0)
	{
		p[0] = p1;
		p[1] = p2;
		p[2] = p3;
		z    = zone;
	}

	/**
	 * Returns vertex x coordinate.
	 * @param i - vertex number (from 0 to 2)
	 * @param ps - mesh points
	 * @return x coordinate
	 */
	double X(int i, const std::vector < MeshPoint > & ps) const 
	{
		return ps[p[i]].x(z);
	}

	/**
	 * Returns vertex y coordinate.
	 * @param i - vertex number (from 0 to 2)
	 * @param ps - mesh points
	 * @return y coordinate
	 */
	double Y(int i, const std::vector < MeshPoint > & ps) const 
	{
		return ps[p[i]].y(z);
	}

	/**
	 * Initialize arrays x and y of mesh points array.
	 * @param ps - mesh points
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
	 * Returns first order finite elements.
	 * @return first order finite elements
	 */
	const std::vector < Polynom > & elem1() const
	{
		return phik;
	}

	/**
	 * Returns finite element in point p1.
	 * @param p1 - point number
	 * @return finite element
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
	typedef std::vector < Triangle > triangles_t;///<triangles container
	typedef std::vector < MeshPoint > points_t;  ///<points container
	typedef std::vector < int > points_flags_t;  ///<points properties container

	triangles_t tr; ///<triangles array
	points_t ps;    ///<points array

	/**
	 * Properties array/
	 *  - 0 - inner point
	 *  - 1 - boundary point
	 */
	points_flags_t ps_flags;

	/**
	 * mapping: point -> triangle in point.
	 */
	std::vector < std::vector < int > > adj;
	/**
	 * Sequence numbers of inner points.
	 */
	std::vector < int > inner;
	/**
	 * Sequence numbers of boundary points.
	 */
	std::vector < int > outer;
	/**
	 * mapping: global point number -> inner point number or outer point number.
	 */
	std::vector < int > p2io;

	/**
	 * Load mesh from file.
	 * @param f - file
	 * @return true if success
	 */
	bool load(FILE * f);

	/**
	 * Prepare all of mesh triangles.
	 */
	void prepare();

	/**
	 * Print mesh information to stdout.
	 */
	void info();
};

/** @} */ /* main */

/**
 * @defgroup print Mesh functions output.
 * @ingroup main
 * @{
 */

/**
 * callback that converts local coordinates to global coordinates.
 */
typedef double (* x_t)(double u, double v);

/**
 * Print a mesh function to file.
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
 * Print a mesh function to file.
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
 * Print the inner part of a mesh function to file.
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
 * Print the inner part of a mesh function to file.
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
 * @defgroup scalar Inner products, norms, distances.
 * @ingroup main
 * @{
 */

/**
 * Calculate inner product of two basis functions on flat domain.
 * @param phi_i - basis function
 * @param phi_j - basis function
 * @param trk - triangle
 * @param m - mesh
 * @param i1 - point number
 * @param j1 - point number
 * @param i2 - inner point number
 * @param j2 - inner point number
 * @param user_data - user data
 */
double generic_scalar_cb(const Polynom & phi_i, const Polynom & phi_j,
						 const Triangle & trk, const Mesh & m, int i1, int j1,
						 int i2, int j2, void * user_data);

/**
 * Calculate inner product of two basis functions on sphere.
 * @param phi_i - basis function
 * @param phi_j - basis function
 * @param trk - triangle
 * @param m - mesh
 * @param i1 - point number
 * @param j1 - point number
 * @param i2 - inner point number
 * @param j2 - inner point number
 * @param user_data - user data
 */
double sphere_scalar_cb(const Polynom & phi_i, const Polynom & phi_j,
						const Triangle & trk, const Mesh & m,
						int i1, int j1, int i2, int j2, void * user_data);

class Matrix;

/**
 * Fast mesh inner product calculator.
 * @see generate_scalar_matrix
 * @param u - mesh function
 * @param v - mesh function
 * @param m - mesh
 * @param mat - matrix
 * @return (u, v)
 */
double fast_scalar(const double * u, const double * v,
				   const Mesh & m, Matrix & mat);

/**
 * Fast mesh norm calculator.
 * @see generate_scalar_matrix
 * @param u - mesh function
 * @param m - mesh
 * @param mat - matrix
 * @return ||u||
 */
double fast_norm(const double * u, const Mesh & m, Matrix & mat);

/**
 * Fast mesh distance calculator.
 * @see generate_scalar_matrix
 * @param u - mesh function
 * @param v - mesh function
 * @param m - mesh
 * @param mat - matrix
 * @return distance between u and v
 */
double fast_dist(const double * u, const double * v,
				 const Mesh & m, Matrix & mat);

/** @} */ /* scalar */

/**
 * @defgroup proj Converters and projectors.
 * @ingroup main
 * @{
 */

/**
 * Function.
 * @param x - x coordinate 
 * @param y - y coordinate
 * @return value
 */
typedef double (* f_xy_t)(double x, double y);

/**
 * Function.
 * @param x - x coordinate
 * @param y - y cooddinate
 * @param t - time
 * @return value
 */
typedef double (* f_xyt_t)(double x, double y, double t);

/**
 * Add boundary conditions to p.
 * @param p - the value of u on the inner mesh points
 * @param u - (output) mesh function
 * @param bnd - the value of u on the boundary mesh points
 * @param m - mesh
 */
void p2u(double * u, const double * p, const double * bnd, const Mesh & m);

/**
 * Remove boundary conditions from u.
 * @param p - (output) the value of u on the inner mesh points
 * @param u - mesh function
 * @param m - mesh
 */
void u2p(double * p, const double * u, const Mesh & m);

/**
 * Project function f(x,y) to the mesh.
 * @param F - (output) the value of f(x,y) on mesh points
 * @param mesh - the mesh
 * @param f - function f(x, y)
 */
void proj(double * F, const Mesh & mesh, f_xy_t f);

/**
 * Project function f(x,y) to the boundary of the mesh.
 * @param F - (output) the value of f(x,y) on boundary points
 * @param m - the mesh
 * @param f - function f(x, y)
 */
void proj_bnd(double * F, const Mesh & m, f_xy_t f);

/**
 * Project vector F1 to the boundary of the mesh.
 * @param F  - (output) the value of F1 on boundary points
 * @param m  - the mesh
 * @param F1 - mesh vector 
 */
void proj_bnd(double * F, const double * F1, const Mesh & m);

/**
 * Set the boundary value of vector F.
 * @param F - (input/output) mesh vector
 * @param bnd - boundary value
 * @param m - mesh  
 */
void set_bnd(double * F, const double * bnd, const Mesh & m);

/**
 * Project function f(x,y,t) to the mesh.
 * @param F - (output) the value of f(x,y) on mesh points
 * @param mesh - the mesh
 * @param f - function f(x, y)
 * @param t - time
 */
void proj(double * F, const Mesh & mesh, f_xyt_t f, double t);

/**
 * Project function f(x,y,t) to the boundary of the mesh.
 * @param F - (output) the value of f(x,y) on boundary points
 * @param m - the mesh
 * @param f - function f(x, y)
 * @param t - time
 */
void proj_bnd(double * F, const Mesh & m, f_xyt_t f, double t);

/** @} */ /* proj */

}

#include "phelm_generators.h"

#endif /* MKE_H */

