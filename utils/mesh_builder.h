#ifndef MESH_BUILDER_H
#define MESH_BUILDER_H
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

#include <vector>
#include <stdio.h>
#include <math.h>

struct Point {
	double x;
	double y;
	double z;

	Point(double x1, double y1, double z1) : x(x1), y(y1), z(z1) {}
	Point(double * x1) : x(x1[0]), y(x1[1]), z(x1[2]) {}
	Point(): x(0), y(0), z(0) {}

	void print() {
		fprintf(stdout, "[%lf, %lf, %lf]\n", x, y, z);
	}

	double abs() const {
		return sqrt(x*x + y*y + z*z);
	}

	/**
	 * counterclockwise rotation about the positive z-axis by angle a
	 */
	Point rotate_x(double a) const {
		double sina = sin(a);
		double cosa = cos(a);
		return Point(
			x,
			y * cosa - z * sina,
			y * sina + z * cosa
			);
	}

	Point rotate_y(double a) const {
		double sina = sin(a);
		double cosa = cos(a);
		return Point(
			x * cosa - z * sina,
			y,
			x * sina + z * cosa
			);
	}

	Point rotate_z(double a) const {
		double sina = sin(a);
		double cosa = cos(a);
		return Point(
			x * cosa - y * sina,
			x * sina + y * cosa,
			z
			);
	}

	double scalar(const Point &a) const {
		return x*a.x + y*a.y + z*a.z;
	}
};

inline Point operator + (const Point & a, const Point & b)
{
	return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Point operator / (const Point & a, double k)
{
	return Point(a.x / k, a.y / k, a.z / k);
}

struct Triangle {
	int v1;
	int v2;
	int v3;

	Triangle(int v1_, int v2_, int v3_) : v1(v1_), v2(v2_), v3(v3_)
	{
	}
};

typedef void (* surface_projector)(Point & p);

void sphere_orto_projector(Point & p);
void sphere_z_projector(Point & p);
void flat_projector(Point & p);

void build_icosahedron(std::vector < Triangle > & r, std::vector < Point > & p);

void normalize_mesh(std::vector < Triangle > & mesh, std::vector < Point > & points, surface_projector project);
void iterate_mesh(std::vector < Triangle > & mesh, std::vector < Point > & points, int iterations, surface_projector project);

typedef bool (* filter_condition)(Point & p);

void filter_mesh(std::vector < Triangle > & mesh,
				 std::vector < Point > & points,
				 std::vector < int > & boundary,
				 filter_condition what);

void reorder_mesh(std::vector < Triangle > & mesh,
	std::vector < Point > & points,
	std::vector < int > & boundary);

void vizualize_adj(std::vector < Triangle > & tri,
	std::vector < Point > & points,
	std::vector < int > & boundary, 
	FILE * f,
	int w,
	int h);

#endif /* MESH_BUILDER_H */
