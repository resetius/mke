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
 
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vector>

#include "mke.h"
#include "util.h"
#include "laplace.h"

using namespace std;
using namespace phelm;

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [mesh.txt|-]\n", name);
	exit(1);
}

double u(double x, double y)
{
	return sin(M_PI * x) * sin(M_PI * y) + 1.0;
}

double v(double x, double y)
{
	return x * (x - 1) * y * (y - 1) + 1.0;
}

double f(double x, double y)
{
	return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) 
		+ v(x, y);
}

double g(double x, double y)
{
	return 2.0 * y * y - 2.0 * y + 2.0 * x * x - 2.0 * x
		+ u(x, y);
}

static elements_t
laplace_integrate_cb( const Polynom & phi_i,
                      const Polynom & phi_j, 
                      const Triangle & trk, /* номер треугольника */
                      const Mesh & m,
                      int point_i, int point_j,
                      int i, int j,
                      void * user_data)
{
	int rs = (int)m.inner.size();
	elements_t r;
	double a;

	a = laplace(phi_i, phi_j, trk, m.ps);
	// Delta u
	r.push_back(Element(i,      j,      a));
	// Delta v
	r.push_back(Element(i + rs, j + rs, a));

	a = integrate(phi_i * phi_j, trk, m.ps);
	// v
	r.push_back(Element(i,      j + rs, a));
	// u
	r.push_back(Element(i + rs, j,      a));
	return r;
}

struct laplace_right_part_cb_data
{
	double * F;
	double * G;

	double * BU;
	double * BV;
};

static elements_t 
laplace_right_part_cb( const Polynom & phi_i,
                       const Polynom & phi_j,
                       const Triangle & trk, /* номер треугольника */
                       const Mesh & m,
                       int point_i, int point_j,
                       int i, int j,
                       laplace_right_part_cb_data * d)
{
	int rs = (int)m.inner.size();
	elements_t r;
	const double * F = d->F;
	const double * G = d->G;

	if (m.ps_flags[point_j] == 1) {     // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const double * BU = d->BU;
		const double * BV = d->BV;
		double a = laplace(phi_i, phi_j, trk, m.ps);
		double b = integrate(phi_i * phi_j, trk, m.ps);

		r.push_back(Element(i,      j, -(BU[j0] * a + BV[j0] * b)));
		r.push_back(Element(i + rs, j, -(BV[j0] * a + BU[j0] * b)));
	} else {
		double a = integrate(phi_i * phi_j, trk, m.ps);
		// F
		r.push_back(Element(i,      j, F[point_j] * a));
		// G
		r.push_back(Element(i + rs, j, G[point_j] * a));
	}

	return r;
}

/**
 * \Delta u + v = f(x, y)
 * u + \Delta v = g(x, y)
 */
void test_invert(Mesh & m)
{
	int sz = (int)m.ps.size();
	int rs = (int)m.inner.size();
	int os = (int)m.outer.size();

	// right part
	vector < double > F(sz);
	vector < double > G(sz);

	// real answer
	vector < double > RU(sz);
	vector < double > RV(sz);

	// answer
	vector < double > U(sz);
	vector < double > V(sz);

	// boundary
	vector < double > BU(os);
	vector < double > BV(os);

	// right part
	proj(&F[0], m, f);
	proj(&G[0], m, g);

	// real answer
	proj(&RU[0], m, u);
	proj(&RV[0], m, v);

	// boundary
	proj_bnd(&BU[0], m, u);
	proj_bnd(&BV[0], m, v);

	Matrix A(2 * rs);
	vector < double > RP(2 * rs);
	vector < double > Ans(2 * rs);

	generate_matrix(A, m, laplace_integrate_cb, (void*)0);

	laplace_right_part_cb_data data;
	data.F = &F[0];
	data.G = &G[0];
	data.BU = &BU[0];
	data.BV = &BV[0];
	generate_right_part(&RP[0], m, laplace_right_part_cb, &data);

	A.solve(&Ans[0], &RP[0]);
	p2u(&U[0], &Ans[0],  &BU[0], m);
	p2u(&V[0], &Ans[rs], &BV[0], m);

	fprintf(stdout, "answer nev: U = %le\n", dist(&U[0], &RU[0], m));
	fprintf(stdout, "answer nev: V = %le\n", dist(&V[0], &RV[0], m));
}

int main(int argc, char *argv[])
{
	Mesh mesh;
	if (argc > 1) {
		FILE * f = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "rb");
		if (!f) {
			usage(argv[0]);
		}
		if (!mesh.load(f)) {
			usage(argv[0]);
		}
		fclose(f);
	} else {
		usage(argv[0]);
	}

	mesh.info();
	test_invert(mesh);
	return 0;
}

