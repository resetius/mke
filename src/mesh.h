#ifndef MESH_H
#define MESH_H

/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2015 Alexey Ozeritsky
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

#include <math.h>

#include "polynom.h"
#include "linal.h"
#include "func.h"
#include "point.h"

/**
 * @namespace phelm
 *
 * This namespace contains all Phelm functions and classes.
 */
namespace phelm
{

/**
 * @defgroup main Mesh and mesh functions.
 * @{
 */

/**
 * A point of a manifold can be included into multiple subdomains.
 * That class supports this.
 */
struct MeshPoint
{
	/**
	 * global coordinates
	 */
	Point pr;

	/**
	 * p[i] local coordinates of point in subdomain i.
	 * @todo what should we do if the point is included only in one subdomain
	 * with sequence number greater that 0 ?
	 */
	std::vector < Point > p;

	enum {
		POINT_REGULAR  = 0,
		POINT_BOUNDARY = 1
	};

	int flags;

	bool is_regular() const {
		return flags == POINT_REGULAR;
	}

	bool is_boundary() const {
		return flags & POINT_BOUNDARY;
	}

	/**
	 * Default constructor.
	 */
	MeshPoint(): flags(0) {}

	/**
	 * Initialization of x and y.
	 * @param x - x coordinate
	 * @param y - y coordinate
	 */
	MeshPoint(double x, double y) : flags(0)
	{
		add (Point (x, y) );
	}

	/**
	 * Add local coordinates for next subdomain.
	 * @param p1 - point in local coordinates of next subdomain.
	 */
	void add (const Point & p1)
	{
		p.push_back (p1);
	}

	/**
	 * Returns x coordinate in subdomain.
	 * @param zone - sequence number of subdomain
	 * @return x coordiante
	 */
	double x (int zone = 0) const
	{
		assert(zone < (int) p.size());
		return p[zone].x;
	}

	/**
	 * Returns y coordinate in subdomain.
	 * @param zone - a sequence number of subdomain
	 * @return y coordiante
	 */
	double y (int zone = 0) const
	{
		assert(zone < (int) p.size());
		return p[zone].y;
	}
};

/**
 * Triangle class.
 */
struct Triangle
{
	int p[3];  ///< point number
	int z;     ///< default zone number

	typedef Polynom elem_t;
	typedef std::vector < Polynom > basis_t;

	const std::vector < MeshPoint > & ps;

	mutable std::vector < basis_t > phik; ///< basis functions in zone 	

	struct NewElem {
		// f(x1, y1, z1), z1 = 0
		FuncPtr f;

		// h, h1 depends on triangle normale
		struct {
			// x1<-(x, y, z)
			FuncPtr hx;
			// y1<-(x, y, z)
			FuncPtr hy;
			// z1<-(x, y, z)
			FuncPtr hz;
		} h;
		
		struct {
			// x<-(x1, y1, z1)
			FuncPtr h1x;
			// y<-(x1, y1, z1)
			FuncPtr h1y;
			// z<-(x1, y1, z1)
			FuncPtr h1z;
		} h1;

		// g, g1 depends on zone
		struct {
			// x<-(phi, la)
			FuncPtr gx;
			// y<-(phi, la)
			FuncPtr gy;
			// z<-(phi, la)
			FuncPtr gz;
		} g;

		struct {
			// phi<-(x, y, z)
			FuncPtr g1phi;
			// la<-(x, y, z)
			FuncPtr g1la;
		} g1;
	};

	mutable std::vector < std::vector<NewElem> > newphi;

	/**
	 * Initialization of point numbers and subdomain number.
	 * @param p1 - point 1
	 * @param p2 - point 2
	 * @param p3 - point 3
	 * @param zone - subdomain number, the triangle belongs to that subdomain
	 */
	Triangle (
		int p1, int p2, int p3, 
		const std::vector < MeshPoint > & ps, 
		int zone = 0)
		:
		z (zone),
		ps (ps)
	{
		p[0] = p1;
		p[1] = p2;
		p[2] = p3;
		phik.reserve(10); // max 10 zones
	}
 
	Triangle (const Triangle & other):
		z(other.z),
		ps(other.ps),
		phik(other.phik)
	{
		p[0] = other.p[0];
		p[1] = other.p[1];
		p[2] = other.p[2];
	}

	Triangle & operator = (const Triangle & other)
	{
		// do not change ps link
		assert(&ps == &other.ps);
		z = other.z;
		phik = other.phik;
		p[0] = other.p[0];
		p[1] = other.p[1];
		p[2] = other.p[2];
		return *this;
	}

	/**
	 * Returns vertex x coordinate.
	 * @param i - vertex number (from 0 to 2)
	 * @param z - subdomain
	 * @param ps - mesh points
	 * @return x coordinate
	 */
	double x (int i, int z) const
	{
		return ps[p[i]].x (z);
	}

	/**
	 * Returns vertex y coordinate.
	 * @param i - vertex number (from 0 to 2)
	 * @param z - subdomain
	 * @param ps - mesh points
	 * @return y coordinate
	 */
	double y (int i, int z) const
	{
		return ps[p[i]].y (z);
	}

private:
	/**
	 * Initialize arrays x and y of mesh points array.
	 * @param ps - mesh points
	 */
	basis_t prepare_basis(int z) const;
	std::vector<NewElem> prepare_new_basis(int z) const;

public:
	/**
	 * Returns triangle's vertex number (0,1,2)
	 * @param p1 - mesh point number
	 */
	int point_number(int p1) const;

	/**
	 * Returns first order finite elements.
	 * @return first order finite elements
	 */
	const basis_t & elem1(int zone) const;

	/**
	 * Returns finite element in point p1.
	 * @param p1 - point number
	 * @return finite element
	 */
	const Polynom & elem1(int p1, int zone) const;

	const NewElem & new_elem1(int p1, int zone) const;
	std::vector<NewElem> new_elem1(int zone) const;
};

using namespace linal;

/**
 * Mesh class.
 */
struct Mesh
{
	typedef std::vector < Triangle > triangles_t;///<triangles container
	typedef std::vector < MeshPoint > points_t;  ///<points container

	triangles_t tr; ///<triangles array
	points_t ps;    ///<points array

	bool is_regular(int i) const {
		return ps[i].is_regular();
	}

	bool is_boundary(int i) const {
		return ps[i].is_boundary();
	}

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

	int inner_size;
	int outer_size;
	int size;

	struct Device
	{
		ArrayDevice < int > inner;
		ArrayDevice < int > outer;
		ArrayDevice < int > p2io;
	};

	Device d;

	/**
	 * Load mesh from file.
	 * @param f - file
	 * @return true if success
	 */
	bool load (FILE * f);

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
typedef double (* x_t) (double u, double v);

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
void print_function (FILE * to, double * ans, const Mesh & m,
                     x_t x = 0, x_t y = 0, x_t z = 0);
void print_function (FILE * to, float * ans, const Mesh & m,
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
void print_function (const char * fname, double * ans, const Mesh & m,
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
void print_inner_function (FILE * to, double * ans, const Mesh & m,
                           x_t x = 0, x_t y = 0, x_t z = 0);
void print_inner_function (FILE * to, float * ans, const Mesh & m,
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
void print_inner_function (const char * to, double * ans, const Mesh & m,
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
 * @param z - subdomain number
 * @param m - mesh
 * @param i1 - point number
 * @param j1 - point number
 * @param i2 - inner point number
 * @param j2 - inner point number
 * @param user_data - user data
 */
double generic_scalar_cb (const Polynom & phi_i, const Polynom & phi_j,
                          const Triangle & trk, int z,
                          const Mesh & m, int i1, int j1,
                          int i2, int j2, void * user_data);

/**
 * Calculate inner product of two basis functions on sphere.
 * @param phi_i - basis function
 * @param phi_j - basis function
 * @param trk - triangle
 * @param z - subdomain number
 * @param m - mesh
 * @param i1 - point number
 * @param j1 - point number
 * @param i2 - inner point number
 * @param j2 - inner point number
 * @param user_data - user data
 */
double sphere_scalar_cb (const Polynom & phi_i, const Polynom & phi_j,
                         const Triangle & trk, int z, const Mesh & m,
                         int i1, int j1, int i2, int j2, void * user_data);

/**
 * Fast mesh inner product calculator.
 * @see generate_scalar_matrix
 * @param u - mesh function
 * @param v - mesh function
 * @param m - mesh
 * @param mat - matrix
 * @return (u, v)
 */
template < typename T, typename Matrix >
T fast_scalar (const T * u, const T * v,
               const Mesh & m, Matrix & A)
{
	int sz = (int) m.ps.size();
	Array < T, Allocator < T > > tmp (sz);
	A.mult_vector (&tmp[0], v);
	return vec_scalar2 (u, &tmp[0], sz);
}

/**
 * Fast mesh norm calculator.
 * @see generate_scalar_matrix
 * @param u - mesh function
 * @param m - mesh
 * @param mat - matrix
 * @return ||u||
 */
template < typename T, typename Matrix >
T fast_norm (const T * u, const Mesh & m, Matrix & A)
{
	return sqrt (fast_scalar (u, u, m, A) );
}

/**
 * Fast mesh distance calculator.
 * @see generate_scalar_matrix
 * @param u - mesh function
 * @param v - mesh function
 * @param m - mesh
 * @param mat - matrix
 * @return distance between u and v
 */
template < typename T, typename Matrix >
T fast_dist (const T * u, const T * v,
             const Mesh & m, Matrix & A)
{
	int sz  = (int) m.ps.size(); // размерность
	Array < T , Allocator < T > > diff (sz);
	vec_diff (&diff[0], u, v, sz);
	return fast_norm (&diff[0], m, A);
}

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
typedef double (* f_xy_t) (double x, double y);

/**
 * Function.
 * @param x - x coordinate
 * @param y - y cooddinate
 * @param t - time
 * @return value
 */
typedef double (* f_xyt_t) (double x, double y, double t);

/**
 * Add boundary conditions to p.
 * Note for phelm_cu: this functions works with device memory.
 * @param p - the value of u on the inner mesh points
 * @param u - (output) mesh function
 * @param bnd - the value of u on the boundary mesh points
 * @param m - mesh
 */
void p2u (double * u, const double * p, const double * bnd, const Mesh & m);
void p2u (float * u, const float * p, const float * bnd, const Mesh & m);

/**
 * Remove boundary conditions from u.
 * Note for phelm_cu: this functions works with device memory.
 * @param p - (output) the value of u on the inner mesh points
 * @param u - mesh function
 * @param m - mesh
 */
void u2p (double * p, const double * u, const Mesh & m);
void u2p (float * p, const float * u, const Mesh & m);

/**
 * Project function f(x,y) to the mesh.
 * Note for phelm_cu: this functions works with host memory.
 * @param F - (output) the value of f(x,y) on mesh points
 * @param mesh - the mesh
 * @param f - function f(x, y)
 */
void proj (double * F, const Mesh & mesh, f_xy_t f);
void proj (float * F, const Mesh & mesh, f_xy_t f);

/**
 * Project function f(x,y) to the boundary of the mesh.
 * Note for phelm_cu: this functions works with host memory.
 * @param F - (output) the value of f(x,y) on boundary points
 * @param m - the mesh
 * @param f - function f(x, y)
 */
void proj_bnd (double * F, const Mesh & m, f_xy_t f);
void proj_bnd (float * F, const Mesh & m, f_xy_t f);

/**
 * Project function f(x,y,t) to the mesh.
 * Note for phelm_cu: this functions works with host memory.
 * @param F - (output) the value of f(x,y) on mesh points
 * @param mesh - the mesh
 * @param f - function f(x, y)
 * @param t - time
 */
void proj (double * F, const Mesh & mesh, f_xyt_t f, double t);
void proj (float * F, const Mesh & mesh, f_xyt_t f, double t);

/**
 * Project function f(x,y,t) to the boundary of the mesh.
 * Note for phelm_cu: this functions works with host memory.
 * @param F - (output) the value of f(x,y) on boundary points
 * @param m - the mesh
 * @param f - function f(x, y)
 * @param t - time
 */
void proj_bnd (double * F, const Mesh & m, f_xyt_t f, double t);
void proj_bnd (float * F, const Mesh & m, f_xyt_t f, double t);

/**
 * Project vector F1 to the boundary of the mesh.
 * Note for phelm_cu: this functions works with device memory.
 * @param F  - (output) the value of F1 on boundary points
 * @param m  - the mesh
 * @param F1 - mesh vector
 */
void proj_bnd (double * F, const double * F1, const Mesh & m);
void proj_bnd (float * F, const float * F1, const Mesh & m);

/**
 * Set the boundary value of vector F.
 * Note for phelm_cu: this functions works with device memory.
 * @param F - (input/output) mesh vector
 * @param bnd - boundary value
 * @param m - mesh
 */
void set_bnd (double * F, const double * bnd, const Mesh & m);
void set_bnd (double * F, const float * bnd, const Mesh & m);

/** @} */ /* proj */

void smooth1 (double * out, const double * in, const Mesh & m);
void smooth2 (double * out, const double * in, const Mesh & m);


/**
 * Solve the system with A matrix (Ax=rp).
 * (Helper function)
 * The function founds an answer on the inner part of the domain
 * and then sets boundary value of the answer to bnd
 *
 * @param answer - the answer
 * @param bnd - boundary
 * @param rp - right part
 * @param A - the matrix of the system
 * @param m - mesh
 */
template < typename T, typename Matrix, typename Mesh >
void solve (T * answer, const T * bnd,
            T * rp, Matrix & A, const Mesh & m)
{
        int sz  = (int) m.ps.size();
        int rs  = (int) m.inner.size();    // размерность
        Array < T, Allocator < T > > x (rs);     // ответ
        solve2 (&x[0], rp, A, m);
        p2u (answer, &x[0], bnd, m);
}

/**
 * Solve the system with A matrix (Ax=rp).
 * (Helper function)
 * Found an answer on the inner part of the domain.
 *
 * @param answer the answer
 * @param rp the right part
 * @param A the matrix of the system
 * @param m the mesh
 */
template < typename T, typename Matrix, typename Mesh >
void solve2 (T * answer, T * rp, Matrix & A, const Mesh & m)
{
        int sz  = (int) m.ps.size();
        int rs  = (int) m.inner.size();    // размерность
        A.solve (answer, &rp[0]);
}

}

#endif /* MESH_H */

