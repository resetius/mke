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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <functional>

#include "polynom.h"
#include "phelm.h"
#include "util.h"
#include "ver.h"

VERSION ("$Id$");

#undef max
#undef min

using namespace std;

namespace phelm
{

#define off(i, j) ((i) * (y_deg_ + 1) + j)

double Polynom::apply (double x, double y) const
{
	short i, j;
	double r = 0.0;
	double * xi = (double*) alloca ( (x_deg_ + 1) * sizeof (double) );
	double * yj = (double*) alloca ( (y_deg_ + 1) * sizeof (double) );

	for (i = 0; i <= x_deg_; i++)
	{
		xi[i] = ipow (x, i);
	}

	for (i = 0; i <= y_deg_; i++)
	{
		yj[i] = ipow (y, i);
	}

	for (i = 0; i <= x_deg_; i++)
	{
		for (j = 0; j <= y_deg_; j++)
		{
			double k = koef_[off (i, j) ];
			r += k * xi[i] * yj[j];
		}
	}
	return r;
}

//0,0 ; 0,1; 1,0; 1,1

double p2x [] = {0.0, 1.0};
double p2y [] = {0.0, 1.0};

const Polynom P2X (1, 0, p2x, 2); //p(x, y) = x
const Polynom P2Y (0, 1, p2y, 2); //p(x, y) = y

#ifdef WIN32
#define snprintf _snprintf
#endif

std::string Polynom::print() const
{
	short i, j;
	std::string ans;
#define BUF_SIZ 32768
	char buf[BUF_SIZ];
	for (i = 0; i <= x_deg_; i++)   // x
	{
		for (j = 0; j <= y_deg_; j++)   // y
		{
			if (fabs (koef_[off (i, j) ]) > 1e-8)
			{
				if (i && j)
				{
					snprintf (buf, BUF_SIZ, "%+.1lf x^%d y^%d ",
					          koef_[off (i, j) ],
					          i, j);
				}
				else if (i)
				{
					snprintf (buf, BUF_SIZ, "%+.1lf x^%d ",
					          koef_[off (i, j) ],
					          i);
				}
				else if (j)
				{
					snprintf (buf, BUF_SIZ, "%+.1lf y^%d ",
					          koef_[off (i, j) ],
					          j);
				}
				else
				{
					snprintf (buf, BUF_SIZ, "%+.1lf ",
					          koef_[off (i, j) ]);
				}
				ans += buf;
			}
		}
	}
	fprintf (stderr, "%s\n", ans.c_str() );
	return ans;
#undef BUF_SIZ
}

Polynom diff (const Polynom & p, int d)
{
	short i, j;
	short x_deg = p.x_deg_;
	short y_deg = p.y_deg_;
	short deg   = y_deg;

	assert (d == 0 || d == 1);

	if (d == 0)   //производная по x
	{
		x_deg -= 1;
		Polynom r (x_deg, y_deg);
		for (j = 0; j <= y_deg; ++j)   //y
		{
			for (i = 0; i <= x_deg; ++i)   //x
			{
				r.koef_[i * (y_deg + 1) + j] =
				    (double) (i + 1) * p.koef_[ (i + 1) * (deg + 1) + j];
			}
		}
		return r;
	}
	else   //d == 1 производная по y
	{
		y_deg -= 1;
		Polynom r (x_deg, y_deg);
		for (i = 0; i <= x_deg; ++i)   //x
		{
			for (j = 0; j <= y_deg; ++j)   //y
			{
				r.koef_[i * (y_deg + 1) + j] =
				    (double) (j + 1) * p.koef_[i * (deg + 1) + j + 1];
			}
		}
		return r;
	}
}

typedef double (*t_int) (int, int, double, double, double, double, double, double);

template < typename T >
double
integrate1 (const Polynom & p, const T & tr, int z, t_int trapezoid)
{
	double k1, b1, k2, b2, k3, b3, t;

	double x1, x2, x3;
	double y1, y2, y3;
	short k, n, x_deg, y_deg;

	double int1 = 0.0, int2 = 0.0;

	x1 = tr.x(0, z);
	x2 = tr.x(1, z);
	x3 = tr.x(2, z);
	y1 = tr.y(0, z);
	y2 = tr.y(1, z);
	y3 = tr.y(2, z);

	if (x1 <= x2 && x2 <= x3)
	{
		t = x2;
		x2 = x3;
		x3 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}
	else if (x2 <= x1 && x1 <= x3)
	{
		t = x1;
		x1 = x2;
		x2 = t;
		t = x2;
		x2 = x3;
		x3 = t;

		t = y1;
		y1 = y2;
		y2 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}
	else if (x2 <= x3 && x3 <= x1)
	{
		t = x1;
		x1 = x2;
		x2 = t;
		t = y1;
		y1 = y2;
		y2 = t;
	}
	else if (x3 <= x1 && x1 <= x2)
	{
		t = x1;
		x1 = x3;
		x3 = t;
		t = y1;
		y1 = y3;
		y3 = t;
	}
	else if (x3 <= x2 && x2 <= x1)
	{
		t = x1;
		x1 = x3;
		x3 = t;
		t = x2;
		x2 = x3;
		x3 = t;

		t = y1;
		y1 = y3;
		y3 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}

	k1 = (y2 - y1) / (x2 - x1);
	b1 = (y1 * x2 - x1 * y2) / (x2 - x1);

	x_deg = p.x_deg_;
	y_deg = p.y_deg_;

	if (fabs (x1 - x3) > 1e-15)
	{
		k3 = (y3 - y1) / (x3 - x1);
		b3 = (y1 * x3 - x1 * y3) / (x3 - x1);

		for (k = 0; k <= x_deg; ++k)   //x
		{
			for (n = 0; n <= y_deg; ++n)   //y
			{
				int1 += p.koef_[k * (y_deg + 1) + n] *
				        trapezoid (k, n, k1, b1, k3, b3, x1, x3);
			}
		}
	}
	else
	{
		int1 = 0;
	}

	if (fabs (x3 - x2) > 1e-15)
	{
		k2 = (y2 - y3) / (x2 - x3);
		b2 = (y3 * x2 - x3 * y2) / (x2 - x3);

		for (k = 0; k <= x_deg; ++k)   //x
		{
			for (n = 0; n <= y_deg; ++n)   //y
			{
				int2 += p.koef_[k * (y_deg + 1) + n] *
				        trapezoid (k, n, k1, b1, k2, b2, x3, x2);
			}
		}
	}
	else
	{
		int2 = 0;
	}

	if (y3 >= y1 && y3 >= y2) return int1 + int2;
	else if (y3 <= y1 && y3 <= y2) return -int1 - int2;
	else if (y3 >= k1 * x3 + b1) return int1 + int2;
	else return -int1 - int2;
}

static double integrate_ (const Polynom & p, const Triangle & tr, int z, t_int trapezoid)
{
	return integrate1 (p, tr, z, trapezoid);
}

double integrate (const Polynom & p, const Triangle & tr, int z)
{
	return integrate_ (p, tr, z, trapezoid_integral);
}

double integrate_cos (const Polynom & p, const Triangle & tr, int z)
{
	return integrate_ (p, tr, z, trapezoid_integral_cos);
}

double integrate_sin (const Polynom & p, const Triangle & tr, int z)
{
	return integrate_ (p, tr, z, trapezoid_integral_sin);
}

double integrate_1_cos (const Polynom & p, const Triangle & tr, int z)
{
	return integrate_ (p, tr, z, trapezoid_integral_1_cos);
}

struct mybind2d
{
	fxy_t func;

	double k1;
	double b1;

	double k2;
	double b2;

	void * data;

	mybind2d(fxy_t func, double k1, double b1, double k2, double b2, void * data):
		func (func), k1 (k1), b1(b1), k2(k2), b2(b2), data (data)
	{
	}

	double apply_x(double x)
	{
		return func(x, k2 * x + b2, data) - func(x, k1 * x + b1, data);
	}

	double apply_y(double y)
	{
		return func(k2 * y + b2, y, data) - func(k1 * y + b1, y, data);
	}
};

static double boundary_func_x (double x, mybind2d * b)
{
	return b->apply_x(x);
}

static double boundary_func_y (double y, mybind2d * b)
{
	return b->apply_y(y);
}

double
integrate_boundary_x (const Triangle & tr, int z, fxy_t func, void * data)
{
	double k1, b1, k2, b2, k3, b3, t;

	double x1, x2, x3;
	double y1, y2, y3;

	double int1 = 0.0, int2 = 0.0;

	x1 = tr.x(0, z);
	x2 = tr.x(1, z);
	x3 = tr.x(2, z);
	y1 = tr.y(0, z);
	y2 = tr.y(1, z);
	y3 = tr.y(2, z);

	if (x1 <= x2 && x2 <= x3)
	{
		t = x2;
		x2 = x3;
		x3 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}
	else if (x2 <= x1 && x1 <= x3)
	{
		t = x1;
		x1 = x2;
		x2 = t;
		t = x2;
		x2 = x3;
		x3 = t;

		t = y1;
		y1 = y2;
		y2 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}
	else if (x2 <= x3 && x3 <= x1)
	{
		t = x1;
		x1 = x2;
		x2 = t;
		t = y1;
		y1 = y2;
		y2 = t;
	}
	else if (x3 <= x1 && x1 <= x2)
	{
		t = x1;
		x1 = x3;
		x3 = t;
		t = y1;
		y1 = y3;
		y3 = t;
	}
	else if (x3 <= x2 && x2 <= x1)
	{
		t = x1;
		x1 = x3;
		x3 = t;
		t = x2;
		x2 = x3;
		x3 = t;

		t = y1;
		y1 = y3;
		y3 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}

	k1 = (y2 - y1) / (x2 - x1);
	b1 = (y1 * x2 - x1 * y2) / (x2 - x1);

	if (fabs (x1 - x3) > 1e-15)
	{
		k3 = (y3 - y1) / (x3 - x1);
		b3 = (y1 * x3 - x1 * y3) / (x3 - x1);

		// int [k1 x + b1, k3 x + b3][x1, x3] dy dx
		mybind2d bb(func, k1, b1, k3, b3, data);
		int1 = gauss_kronrod15(x1, x3, (fx_t)boundary_func_x, &bb);
	}
	else
	{
		int1 = 0;
	}

	if (fabs (x3 - x2) > 1e-15)
	{
		k2 = (y2 - y3) / (x2 - x3);
		b2 = (y3 * x2 - x3 * y2) / (x2 - x3);

		// int [x3, x2][k1 x + b1, k2 x + b2]
		mybind2d bb(func, k1, b1, k2, b2, data);
		int2 = gauss_kronrod15(x3, x2, (fx_t)boundary_func_x, &bb);
	}
	else
	{
		int2 = 0;
	}

	if (y3 >= y1 && y3 >= y2) return int1 + int2;
	else if (y3 <= y1 && y3 <= y2) return -int1 - int2;
	else if (y3 >= k1 * x3 + b1) return int1 + int2;
	else return -int1 - int2;
}

double
integrate_boundary_y (const Triangle & tr, int z, fxy_t func, void * data)
{
	double k1, b1, k2, b2, k3, b3, t;

	double x1, x2, x3;
	double y1, y2, y3;

	double int1 = 0.0, int2 = 0.0;

	x1 = tr.y(0, z);
	x2 = tr.y(1, z);
	x3 = tr.y(2, z);
	y1 = tr.x(0, z);
	y2 = tr.x(1, z);
	y3 = tr.x(2, z);

	if (x1 <= x2 && x2 <= x3)
	{
		t = x2;
		x2 = x3;
		x3 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}
	else if (x2 <= x1 && x1 <= x3)
	{
		t = x1;
		x1 = x2;
		x2 = t;
		t = x2;
		x2 = x3;
		x3 = t;

		t = y1;
		y1 = y2;
		y2 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}
	else if (x2 <= x3 && x3 <= x1)
	{
		t = x1;
		x1 = x2;
		x2 = t;
		t = y1;
		y1 = y2;
		y2 = t;
	}
	else if (x3 <= x1 && x1 <= x2)
	{
		t = x1;
		x1 = x3;
		x3 = t;
		t = y1;
		y1 = y3;
		y3 = t;
	}
	else if (x3 <= x2 && x2 <= x1)
	{
		t = x1;
		x1 = x3;
		x3 = t;
		t = x2;
		x2 = x3;
		x3 = t;

		t = y1;
		y1 = y3;
		y3 = t;
		t = y2;
		y2 = y3;
		y3 = t;
	}

	k1 = (y2 - y1) / (x2 - x1);
	b1 = (y1 * x2 - x1 * y2) / (x2 - x1);

	if (fabs (x1 - x3) > 1e-15)
	{
		k3 = (y3 - y1) / (x3 - x1);
		b3 = (y1 * x3 - x1 * y3) / (x3 - x1);

		// int [k1 x + b1, k3 x + b3][x1, x3] dy dx
		mybind2d bb(func, k1, b1, k3, b3, data);
		int1 = gauss_kronrod15(x1, x3, (fx_t)boundary_func_y, &bb);
	}
	else
	{
		int1 = 0;
	}

	if (fabs (x3 - x2) > 1e-15)
	{
		k2 = (y2 - y3) / (x2 - x3);
		b2 = (y3 * x2 - x3 * y2) / (x2 - x3);

		// int [x3, x2][k1 x + b1, k2 x + b2]
		mybind2d bb(func, k1, b1, k2, b2, data);
		int2 = gauss_kronrod15(x3, x2, (fx_t)boundary_func_y, &bb);
	}
	else
	{
		int2 = 0;
	}

	if (y3 >= y1 && y3 >= y2) return int1 + int2;
	else if (y3 <= y1 && y3 <= y2) return -int1 - int2;
	else if (y3 >= k1 * x3 + b1) return int1 + int2;
	else return -int1 - int2;
}

extern "C" double cubature7 (double x1, double y1, double x2, double y2, double x3, double y3, fxy_t f, void * d);

double integrate_generic (const Triangle & tr, int z, fxy_t f, void * data)
{
	double x1 = tr.x(0, z), x2 = tr.x(1, z), x3 = tr.x(2, z);
	double y1 = tr.y(0, z), y2 = tr.y(1, z), y3 = tr.y(2, z);
	return cubature7 (x1, y1, x2, y2, x3, y3, f, data);
}

struct ff_data {
	std::function<double(double, double)> f;
};

static double ff(double x, double y, void * d) {
	ff_data * data = (ff_data*)d;
	return data->f(x, y);
}

static double integrate_generic_new_0(
	const Point & p1, const Point & p2, const Point & p3,
	std::function<double(double, double)> f)
{
	double x1 = p1.x, x2 = p2.x, x3 = p3.x;
	double y1 = p1.y, y2 = p2.y, y3 = p3.y;

	ff_data data;
	data.f = f;
	double r = cubature7(x1, y1, x2, y2, x3, y3, ff, &data);
	return r;
}

double integrate_generic_new(
	const Point & p1, const Point & p2, const Point & p3,
	std::function<double(double, double)> f, int iter)
{
	if (iter == 0) {
		return integrate_generic_new_0(p1, p2, p3, f);
	}
	else {
		double s = 0;
		iter--;
		Point p12 = (p1 + p2) / 2;
		Point p23 = (p2 + p3) / 2;
		Point p31 = (p3 + p1) / 2;
		s += integrate_generic_new(p1, p12, p31, f, iter);
		s += integrate_generic_new(p2, p12, p23, f, iter);
		s += integrate_generic_new(p3, p23, p31, f, iter);
		s += integrate_generic_new(p12, p23, p31, f, iter);
		return s;
	}
}

double integrate_generic_new(
	const Triangle & tr, 
	std::function<double(double, double)> f, int iter)
{
	return integrate_generic_new(tr.pp[0], tr.pp[1], tr.pp[2], f, iter);
}

Polynom operator * (const Polynom &p1, const Polynom &p2)
{
	short i, j, k, m;
	short x_deg1 = p1.x_deg_;
	short y_deg1 = p1.y_deg_;
	short x_deg2 = p2.x_deg_;
	short y_deg2 = p2.y_deg_;
	short deg;

	Polynom p (p1.x_deg_ + p2.x_deg_, p1.y_deg_ + p2.y_deg_);
	deg = p.y_deg_;

	// ?
	for (i = 0; i <= x_deg1; ++i)
	{
		for (j = 0; j <= y_deg1; ++j)
		{
			for (k = 0; k <= x_deg2; ++k)
			{
				for (m = 0; m <= y_deg2; ++m)
				{
					p.koef_[ (i + k) * (deg + 1) + (j + m) ] +=
					    p1.koef_[i * (y_deg1 + 1) + j] * p2.koef_[k * (y_deg2 + 1) + m];
				}
			}
		}
	}

	return p;
}

Polynom operator - (const Polynom &p1, const Polynom &p2)
{
	short i, j;
	Polynom r (std::max (p1.x_deg_, p2.x_deg_), std::max (p1.y_deg_, p2.y_deg_) );

	for (i = 0; i <= r.x_deg_; ++i)   //x
	{
		for (j = 0; j <= r.y_deg_; ++j)   //y
		{
			r.koef_[i * (r.y_deg_ + 1) + j]
			= p1.k (i, j) - p2.k (i, j);
		}
	}

	return r;
}

Polynom operator + (const Polynom &p1, const Polynom &p2)
{
	short i, j;
	Polynom r (std::max (p1.x_deg_, p2.x_deg_), std::max (p1.y_deg_, p2.y_deg_) );

	for (i = 0; i <= r.x_deg_; ++i)   //x
	{
		for (j = 0; j <= r.y_deg_; ++j)   //y
		{
			r.koef_[i * (r.y_deg_ + 1) + j]
			= p1.k (i, j) + p2.k (i, j);
		}
	}
	return r;
}
}

