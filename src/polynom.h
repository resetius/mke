#ifndef POLYNOM_H
#define POLYNOM_H
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
 *
 * @section DESCRIPTION
 * The polynom class, methods, functions.
 */

#include <vector>
#include <string>
#include <string.h>

#include <algorithm>

typedef unsigned int uint;

namespace phelm
{

/**
 * @defgroup polynom Polynom functions and classes
 * @{
 */

/**
 * Polynom of variables (x, y).
 * Examples:
  \f[
    x^2 y + y + x + 1,
  \f]
  x degree = 2; y degree = 1,
  \f[
    x^2 + 1,
  \f]
  x degree = 2; y degree = 0.
  \f[
  P(x,y) = \sum_{i=0}^{x\_deg}\sum_{j=0}^{y\_deg}k_{ij}x^i y^j.
  \f]
*/
struct Polynom
{
	short x_deg_;                   ///<x degree
	short y_deg_;                   ///<y degree
	std::vector < double > koef_;   ///<coefficients

	/**
	 * Creates new empty polynom.
	 * @param x_deg - x degree
	 * @param y_deg - y degree
	 */
	Polynom (short x_deg, short y_deg)
			: x_deg_ (x_deg), y_deg_ (y_deg), koef_ ( (x_deg + 1) * (y_deg + 1) )
	{
	}

	/**
	 * Creates new polynom from coefficients.
	 * @param x_deg - x degree
	 * @param y_deg - y degree
	 * @param koef - the coefficients
	 * @param l - the number of elements in the coefficients array
	 */
	Polynom (short x_deg, short y_deg, double * koef, int l)
			: x_deg_ (x_deg), y_deg_ (y_deg), koef_ ( (x_deg + 1) * (y_deg + 1) )
	{
		memcpy (&koef_[0], koef, std::min (l, (int) koef_.size() )
		        * sizeof (double) );
	}

	~Polynom() {}

	/**
	 * Prints polynom into the string.
	 * @return string representation of the polynom.
	 */
	std::string print() const;

	/**
	 * Calculates Polynom(x, y).
	 * @return Polynom(x, y)
	 */
	double apply (double x, double y) const;

	/**
	 * Returns polynom coefficient.
	 * @return \f$k_{ij}\f$
	 */
	double k (int i, int j) const
	{
		if (i > x_deg_ || j > y_deg_) return 0;
		return koef_[i * (y_deg_ + 1) + j];
	}

	/**
	 * Calculates new polynom.
	 * @return \f$\frac{P(x,y)}{k}\f$
	 */
	void operator /= (double k)
	{
		uint sz   = (uint) koef_.size();
		double k1 = 1.0 / k;
		for (uint i = 0; i < sz; ++i)
		{
			koef_[i] *= k1;
		}
	}
};

/**
 * The polynom: \f$P(x,y)=x\f$.
 * @relates Polynom
 */
extern const Polynom P2X; //p(x, y) = x

/**
 * The polynom:  \f$P(x,y)=y\f$.
 * @relates Polynom
 */
extern const Polynom P2Y; //p(x, y) = y

struct Triangle;
struct Point;
struct MeshPoint;

/**
 * A derivative of the polynom with respect to x or y.
 * If i equals to 0 then calculate a derivative with respect to x.
 * If i equals to 0 then calculate a derivative with respect to y.
 * Calculates \f$ \frac{\partial P(x,y)}{\partial x} \f$ or \f$ \frac{\partial P(x,y)}{\partial y} \f$.
 *
 * @relates Polynom
 * @param p - the polynom
 * @param i - 0 or 1 (x or y)
 * @return \f$\frac{\partial P(x,y)}{\partial x}\f$
 * or \f$\frac{\partial P(x,y)}{\partial y}\f$.
 */
Polynom diff (const Polynom & p, int i);

/**
 * Takes the integral of Polynom over the Triangle t.
 \f[
  \int_t p(x,y) dx dy
 \f]
 * @relates Polynom
 * @param p  - polynom
 * @param t  - triangle
 * @param z  - subdomain number
 */
double integrate (const Polynom & p, const Triangle & t, int z);

/**
 * Takes the integral of Polynom over the Triangle t.
 \f[
  \int_t p(x,y) cos(x) dx dy
 \f]
 *
 * @relates Polynom
 * @param p  - polynom
 * @param t  - triangle
 * @param z  - subdomain number
 */
double integrate_cos (const Polynom & p, const Triangle & t, int z);

typedef double (*fxy_t) (double x, double y, void * data);
double integrate_generic (const Triangle & t, int z, fxy_t f, void * data);

/**
 * Takes the integral of Polynom over the Triangle t.
 \f[
  \int_t p(x,y) sin(x) dx dy
 \f]
 *
 * @relates Polynom
 * @param p  - polynom
 * @param t  - triangle
 * @param z  - subdomain number
 */
double integrate_sin (const Polynom & p, const Triangle & t, int z);

/**
 * Takes the integral of Polynom over the Triangle t.
 \f[
  \int_t \frac{p(x,y)}{cos(x)} dx dy
 \f]
 *
 * @relates Polynom
 * @param p  - the polynom
 * @param t  - the triangle
 * @param ps - the mesh points
 */
double integrate_1_cos (const Polynom & p, const Triangle & t, int z);

double
integrate_boundary_x (const Triangle & tr, int z, fxy_t func, void * data);

double
integrate_boundary_y (const Triangle & tr, int z, fxy_t func, void * data);

/**
 * A product of polynom p1 by polynom p2.
 *
 * @relates Polynom
 * @param p1 - the polynom
 * @param p2 - the polynom
 * @return a product of polynoms.
 */
Polynom operator * (const Polynom &p1, const Polynom &p2);

/**
 * A difference of two polynoms.
 *
 * @relates Polynom
 * @param p1 - the polynom
 * @param p2 - the polynom
 * @return p1 - p2
 */
Polynom operator - (const Polynom &p1, const Polynom &p2);

/**
 * A sum of two polynoms.
 *
 * @relates Polynom
 * @param p1 - the polynom
 * @param p2 - the polynom
 * @return p1 + p2
 */
Polynom operator + (const Polynom &p1, const Polynom &p2);

/**
 * A difference between a polynom and a number.
 *
 * @relates Polynom
 * @param p1 - the polynom
 * @param x - the number
 * @return p1 - x
 */
inline Polynom operator - (const Polynom &p1, double x)
{
	Polynom r (p1);
	r.koef_[0] -= x;
	return r;
}

/**
 * A product of a polynom p1 by a number x.
 *
 * @relates Polynom
 * @param p1 - the polynom
 * @param x - the number
 * @return a product of a polynom and a number.
 */
inline Polynom operator * (const Polynom &p1, double x)
{
	Polynom r (p1);
	uint sz = (uint) r.koef_.size();
	for (uint i = 0; i < sz; ++i)
	{
		r.koef_[i] *= x;
	}
	return r;
}

/** @} */

}

#endif /* POLYNOM_H */

