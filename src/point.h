#ifndef _HH_POINT_H
#define _HH_POINT_H

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
#include <stdio.h>

namespace phelm
{

/**
 * @ingroup main Mesh and mesh functions.
 * @{
 */

struct Matrix {
	double m[3][3];

	Matrix();
	Matrix(double m[3][3]);
	void rotate_x(double a);
	void rotate_y(double a);
	void rotate_z(double a);
	void apply(const Matrix & m);
};

/**
 * The Point class represents a point.
 */
struct Point
{
	double x; ///< x coordinate
	double y; ///< y coordinate
	double z; ///< z coordinate

	/** Default constructor. */
	Point() : x(0), y(0), z(0) {}
	/**
	 * The initialization of x1 and y1.
	 * @param x1 - x coordinate
	 * @param y1 - y coordinate
	 * @param z1 - z coordinate
	 */
	Point(double x1, double y1, double z1 = 0) : x(x1), y(y1), z(z1) {}

	Point(double * p, int size) {
		if (size >= 1) {
			x = p[0];
		}
		if (size >= 2) {
			y = p[1];
		}
		if (size >= 3) {
			z = p[2];
		}
	}

	void print(FILE * f = stdout) const;
	double len() const;
	/**
	 * counterclockwise rotation about the positive z-axis by angle a
	 */
	Point rotate_x(double a) const;
	Point rotate_y(double a) const;
	Point rotate_z(double a) const;

	Point apply(const Matrix & m) const;
};

double scalar(const Point & a, const Point & b);

/**
 * Divide each coordinate by k.
 * @param k - a number
 * @return the new point
 */
inline Point operator / (const Point & a, double k)
{
	return Point(a.x / k, a.y / k, a.z / k);
}

/**
 * Multiply each coordinate by k.
 * @param k - a number
 * @return new point
 */
inline Point operator * (const Point & a, double k)
{
	return Point(a.x * k, a.y * k, a.z * k);
}

inline Point operator * (const Point & a, const Point & b)
{
	return Point(
		a.z*b.y - a.y*b.z,
		a.x*b.z - a.z*b.x,
		a.y*b.x - a.x*b.y
		);
}

/**
 * Sum of two points.
 * @relates Point
 * @param p1 - input point
 * @param p2 - input point
 * @return p1 + p2
 */

inline Point operator + (const Point & a, const Point & b)
{
	return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Point operator - (const Point & a, const Point & b)
{
	return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}

/**
 * @}
 */
}


#endif /*_HH_POINT_H*/
