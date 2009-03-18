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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "polynom.h"
#include "mke.h"
#include "util.h"

#undef max
#undef min

using namespace std;

static double
Polynom_apply_1(const Polynom * p, const double * x)
{
	int i;
	double r = 0.0;
	int deg = p->deg_;
	for (i = 0; i <= deg; i++) {
		r += p->koef_[i] * ipow(x[0], i);
	}
	return r;
}

static double
Polynom_apply_2(const Polynom * p, const double * x)
{
	int i, j;
	double r = 0.0;
	int deg = p->deg_;
	for (i = 0; i <= deg; i++) {
		for (j = 0; j <= deg; j++) {
			double k = p->koef_[i * (deg + 1) + j];
			r += k * ipow(x[0], i) * ipow(x[1], j);
		}
	}
	return r;
}

double Polynom::apply(const double * x) const
{
	switch (n_) {
	case 1:
		return Polynom_apply_1(this, x);
		break;
	case 2:
		return Polynom_apply_2(this, x);
		break;
	default:
		fprintf(stderr, "dimension %d is not supported \n", n_);
		exit(1);
	}
	return 0;
}

double Polynom::apply(double x, double y) const
{
	double x1 [] = {x, y};
	return apply(x1);
}

//0,0 ; 0,1; 1,0; 1,1

double p2x [] = {0.0, 0.0, 1.0};
double p2y [] = {0.0, 1.0, 0.0};

const Polynom P2X(1, 2, p2x, 3); //p(x, y) = x
const Polynom P2Y(1, 2, p2y, 3); //p(x, y) = y

static void
print_1(const Polynom & p)
{
	int i;
	for (i = 0; i <= p.deg_; i++) {
		if (fabs(p.koef_[i]) > 1e-8) {
			fprintf(stderr, "%+.1lf x^%d + ", p.koef_[i], i);
		}
	}
	fprintf(stderr, "\n");
}

static void
print_2(const Polynom & p)
{
	int i, j;
	for (i = 0; i <= p.deg_; i++) { // x
		for (j = 0; j <= p.deg_; j++) { // y
			if (fabs(p.koef_[i * (p.deg_ + 1) + j]) > 1e-8) {
				if (i && j) {
					fprintf(stderr, "%+.1lf x^%d y^%d ", 
						p.koef_[i * (p.deg_ + 1) + j], 
						i, j);
				} else if (i) {
					fprintf(stderr, "%+.1lf x^%d ",
							p.koef_[i * (p.deg_ + 1) + j],
							i);
				} else if (j) {
					fprintf(stderr, "%+.1lf y^%d ",
							p.koef_[i * (p.deg_ + 1) + j],
							j);
				} else {
					fprintf(stderr, "%+.1lf ",
							p.koef_[i * (p.deg_ + 1) + j]);
				}
			}
		}
	}
	fprintf(stderr, "\n");
}

void Polynom::print() const
{
	switch (n_) {
	case 1:
		print_1(*this);
		break;
	case 2:
		print_2(*this);
		break;
	default:
		fprintf(stderr, "dimension %d is not supported \n", n_);
		exit(1);
	}
}

Polynom diff(const Polynom & p, int d)
{
	assert(p.n_ == 2);
	int i, j;
	int deg = p.deg_;

	// ?
	Polynom r(p.deg_, p.n_);
	if (d == 0) { //производная по x
		for (j = 0; j <= deg; ++j) { //y
			for (i = 0; i < deg; ++i) { //x
				r.koef_[i * (deg + 1) + j] = 
					(double)(i + 1) * p.koef_[(i + 1) * (deg + 1) + j];
			}
			//степень уменьшилась
			r.koef_[deg * (deg + 1) + j] = 0;
		}
	} else { //d == 1 производная по y
		for (i = 0; i <= deg; ++i) { //x
			for (j = 0; j < deg; ++j) { //y
				r.koef_[i * (deg + 1) + j] = 
					(double)(j + 1) * p.koef_[i * (deg + 1) + j + 1];
			}
			//степень уменьшилась
			r.koef_[i * (deg + 1) + deg] = 0;
		}
	}
	return r;
}

typedef double (*t_int)(int, int, double, double, double, double, double, double);

struct TriangleX {
	Point p[3];

	TriangleX(Point p1_, Point p2_, Point p3_)
	{
		p[0] = p1_;
		p[1] = p2_;
		p[2] = p3_;
	}
};

static double 
integrate1(const Polynom & p, const TriangleX & tr, t_int trapezoid)
{
	double k1, b1, k2, b2, k3, b3, t;

	double x1, x2, x3;
	double y1, y2, y3;
	int k, n, deg;

	double int1 = 0.0, int2 = 0.0;

	x1 = tr.p[0].x; x2 = tr.p[1].x; x3 = tr.p[2].x;
	y1 = tr.p[0].y; y2 = tr.p[1].y; y3 = tr.p[2].y;

	if (x1 <= x2 && x2 <= x3)
	{
		t = x2; x2 = x3; x3 = t;
		t = y2; y2 = y3; y3 = t;
	}
	else if (x2 <= x1 && x1 <= x3)
	{
		t = x1; x1 = x2; x2 = t;
		t = x2; x2 = x3; x3 = t;
		
		t = y1; y1 = y2; y2 = t;
		t = y2; y2 = y3; y3 = t;
	}
	else if (x2 <= x3 && x3 <= x1)
	{
		t = x1; x1 = x2; x2 = t;
		t = y1; y1 = y2; y2 = t;
	}
	else if (x3 <= x1 && x1 <= x2)
	{
		t = x1; x1 = x3; x3 = t;
		t = y1; y1 = y3; y3 = t;
	}
	else if (x3 <= x2 && x2 <= x1)
	{
		t = x1; x1 = x3; x3 = t;
		t = x2; x2 = x3; x3 = t;

		t = y1; y1 = y3; y3 = t;
		t = y2; y2 = y3; y3 = t;
	}

	k1 = (y2 - y1) / (x2 - x1);
	b1 = (y1 * x2 - x1 * y2) / (x2 - x1);

	k2 = (y2 - y3) / (x2 - x3);
	b2 = (y3 * x2 - x3 * y2) / (x2 - x3);

 	k3 = (y3 - y1) / (x3 - x1);
 	b3 = (y1 * x3 - x1 * y3) / (x3 - x1);

	//fprintf(stderr, "    (%.2lf,%.2lf)-(%.2lf,%.2lf)-(%.2lf,%.2lf)\n", x1, y1, x2, y2, x3, y3);
	
	deg = p.deg_;
	for (k = 0; k <= deg; ++k) { //x
		for (n = 0; n <= deg; ++n) { //y
			int1 += p.koef_[k * (deg + 1) + n] * 
				trapezoid(k, n, k1, b1, k3, b3, x1, x3);
		}
	}

	for (k = 0; k <= deg; ++k) { //x
		for (n = 0; n <= deg; ++n) { //y
			int2 += p.koef_[k * (deg + 1) + n] * 
				trapezoid(k, n, k1, b1, k2, b2, x3, x2);
		}
	}

	if (fabs(x1 - x3) < 1e-15) int1 = 0;
	if (fabs(x3 - x2) < 1e-15) int2 = 0;
	
	if (y3 >= y1 && y3 >= y2) return int1 + int2;
	else if (y3 <= y1 && y3 <= y2) return -int1 - int2;
	else if (y3 >= k1 * x3 + b1) return int1 + int2;
	else return -int1 - int2;
}

static double integrate_(const Polynom & p, const Triangle & tr, const vector < MeshPoint > & ps, t_int trapezoid)
{
#if 0
	double x1, x2, x3;
	double y1, y2, y3;

	double k;
	Point p1, p2, p3, p4;

	x1 = tr.x(0, ps); x2 = tr.x(1, ps); x3 = tr.x(2, ps);
	y1 = tr.y(0, ps); y2 = tr.y(1, ps); y3 = tr.y(2, ps);

	if (y1 <= y2 && y2 <= y3) {
		p1 = ps[tr.p[1]]; p2 = ps[tr.p[0]];	p3 = ps[tr.p[2]];
	} else if (y2 <= y1 && y1 <= y3) {
		p1 = ps[tr.p[0]]; p2 = ps[tr.p[1]];	p3 = ps[tr.p[2]];
	} else {
		p1 = ps[tr.p[2]]; p2 = ps[tr.p[1]];	p3 = ps[tr.p[0]];
	}

	p4.y = p1.y;
	k = fabs(p2.y - p4.y) / (p3.y - p4.y);
	p4.x = (p2.x - k * p3.x) / (1.0 - k);

	if (p2.x < p3.x) {
		if (p2.x <= p4.x && p4.x <= p3.x) {
			;
		} else {
			p4.x = -p4.x;
		}
		assert(p2.x <= p4.x && p4.x <= p3.x);
	} else {
		if (p3.x <= p4.x && p4.x <= p2.x) {
			;
		} else {
			p4.x = -p4.x;
		}
		assert(p3.x <= p4.x && p4.x <= p2.x);
	}

	TriangleX t1(p1, p3, p4), t2(p1, p2, p4);

	return integrate1(p, t1) + integrate1(p, t2);
#endif
	int zone = tr.z;
	TriangleX t1(ps[tr.p[0]].p_[zone], ps[tr.p[1]].p_[zone], ps[tr.p[2]].p_[zone]);
	return integrate1(p, t1, trapezoid);
}

double integrate(const Polynom & p, const Triangle & tr, const vector < MeshPoint > & ps)
{
	return integrate_(p, tr, ps, trapezoid_integral);
}

double integrate_cos(const Polynom & p, const Triangle & tr, const vector < MeshPoint > & ps)
{
	return integrate_(p, tr, ps, trapezoid_integral_cos);
}

double integrate_sin(const Polynom & p, const Triangle & tr, const vector < MeshPoint > & ps)
{
	return integrate_(p, tr, ps, trapezoid_integral_sin);
}

double integrate_1_cos(const Polynom & p, const Triangle & tr, const vector < MeshPoint > & ps)
{
	return integrate_(p, tr, ps, trapezoid_integral_1_cos);
}

Polynom polinom_mult_1(const Polynom &p1, const Polynom &p2)
{
	assert(p1.n_ == 1);

	int i, j;
	int deg1 = p1.deg_;
	int deg2 = p2.deg_;
	int deg;

	Polynom p(p1.deg_ + p2.deg_, p1.n_);
	deg = p.deg_;

	// ?
	for (i = 0; i <= deg1; ++i) {
		for (j = 0; j <= deg2; ++j) {
			p.koef_[i + j] += p1.koef_[i] * p2.koef_[j];
		}
	}

	return p;
}

Polynom polinom_mult_2(const Polynom &p1, const Polynom &p2)
{
	assert(p1.n_ == 2);

	int i, j, k, m;
	int deg1 = p1.deg_;
	int deg2 = p2.deg_;
	int deg;

	Polynom p(p1.deg_ + p2.deg_, p1.n_);
	deg = p.deg_;

	// ?
	for (i = 0; i <= deg1; ++i) {
		for (j = 0; j <= deg1; ++j) {
			for (k = 0; k <= deg2; ++k) {
				for (m = 0; m <= deg2; ++m) {
					p.koef_[(i + k) * (deg + 1) + (j + m)] += 
						p1.koef_[i * (deg1 + 1) + j] * p2.koef_[k * (deg2 + 1) + m];
				}
			}
		}
	}

	return p;
}

Polynom operator * (const Polynom &p1, const Polynom &p2)
{
	assert(p1.n_ == p2.n_);
	assert(p1.n_ == 1 || p1.n_ == 2);

	if (p1.n_ == 1) {
		return polinom_mult_1(p1, p2);
	} else {
		return polinom_mult_2(p1, p2);
	}
}

Polynom operator - (const Polynom &p1, const Polynom &p2)
{
	assert(p1.n_ == p2.n_);

	Polynom r(std::max(p1.deg_, p2.deg_), p1.n_);
	for (uint i = 0; i < r.koef_.size(); ++i) {
		r.koef_[i] = p1.koef_[i] - p2.koef_[i];
	}
	return r;
}

Polynom operator + (const Polynom &p1, const Polynom &p2)
{
	assert(p1.n_ == p2.n_);

	Polynom r(std::max(p1.deg_, p2.deg_), p1.n_);
	for (uint i = 0; i < r.koef_.size(); ++i) {
		r.koef_[i] = p1.koef_[i] + p2.koef_[i];
	}
	return r;
}
