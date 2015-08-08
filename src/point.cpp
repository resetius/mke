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

#include <string.h>
#include "point.h"

using namespace phelm;

Matrix::Matrix() {
	memset(m, 0, sizeof(m));
	m[0][0] = m[1][1] = m[2][2] = 1.0;
}

Matrix::Matrix(double m1[3][3]) {
	memcpy(m, m1, sizeof(m));
}

void Matrix::apply(const Matrix & m1) {
	double m2[3][3] = { 0 };
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double s = 0.0;
			for (int k = 0; k < 3; k++) {
				s += m1.m[i][k] * m[k][j];
			}
			m2[i][j] = s;
		}
	}
	memcpy(m, m2, sizeof(m));
}

void Matrix::rotate_z(double a) {
	double sina = sin(a);
	double cosa = cos(a);
	double m1[3][3] = {
		{ cosa, -sina, 0 },
		{ sina,  cosa, 0 },
		{ 0,     0,    1 }
	};
	apply(Matrix(m1));
}

void Matrix::rotate_y(double a) {
	double sina = sin(a);
	double cosa = cos(a);
	double m1[3][3] = {
		{ cosa, 0, -sina },
		{ 0,    1,  0    },
		{ sina, 0,  cosa }
	};
	apply(Matrix(m1));
}

void Matrix::rotate_x(double a) {
	double sina = sin(a);
	double cosa = cos(a);
	double m1[3][3] = {
		{ 1, 0,     0    },
		{ 0, cosa, -sina },
		{ 0, sina,  cosa }		
	};
	apply(Matrix(m1));
}

void Point::print(FILE * f) const {
	fprintf(f, "[%lf, %lf, %lf]\n", x, y, z);
}

double Point::len() const {
	return sqrt(scalar(*this, *this));
}

Point Point::rotate_x(double a) const {
	double sina = sin(a);
	double cosa = cos(a);
	return Point(
		x,
		y * cosa - z * sina,
		y * sina + z * cosa
		);
}

Point Point::rotate_y(double a) const {
	double sina = sin(a);
	double cosa = cos(a);
	return Point(
		x * cosa - z * sina,
		y,
		x * sina + z * cosa
		);
}

Point Point::rotate_z(double a) const {
	double sina = sin(a);
	double cosa = cos(a);
	return Point(
		x * cosa - y * sina,
		x * sina + y * cosa,
		z
		);
}

Point Point::apply(const Matrix & m) const {
	double xx[3] = { x, y, z };
	double yy[3] = { 0 };
	for (int i = 0; i < 3; ++i) {
		double s = 0;
		for (int k = 0; k < 3; k++) {
			s += m.m[i][k] * xx[k];
		}
		yy[i] = s;
	}
	return Point(yy, 3);
}

double phelm::scalar(const Point & a, const Point & b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
