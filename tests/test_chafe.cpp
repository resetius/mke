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

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [mesh.txt|-]\n", name);
	exit(1);
}

double ans(double x, double y, double t)
{
	return exp(t) * x * (x - 1) * y * (y - 1);
}

double bnd(double x, double y, double t)
{
	return ans(x, y, t);
}

template < typename T >
void init_bnd(Mesh & m, vector < double > & F, T f, double t)
{
	F.resize(m.outer.size());
	for (size_t i = 0; i < m.outer.size(); ++i) {
		int p0 = m.outer[i];
		Point & p = m.ps[p0];
		F[i] = f(p.x, p.y, t);
	}
}

template < typename T >
void init_func(Mesh & mesh, vector < double > & F, T f, double t)
{
	F.resize(mesh.ps.size());
	for (size_t i = 0; i < mesh.ps.size(); ++i)
	{
		F[i] = f(mesh.ps[i].x, mesh.ps[i].y, t);
	}
}

double nr2(double * a, double * b, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt(sum);
}

int main(int argc, char *argv[])
{
	Mesh mesh;
	int i, steps = 1;
	double tau   = 0.01;
	double mu    = 1.0;
	double sigma = 0.0;//-70;

	vector < double > U;
	vector < double > B;
	vector < double > Ans;
	vector < double > P;

	if (argc > 1) {
		FILE * f = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "rb");
		if (!f) {
			usage(argv[0]);
		}
		mesh.load(f);
		fclose(f);
	} else {
		usage(argv[0]);
	}

	init_func(mesh, U, ans, 0.0);
	Ans.resize(U.size());
	P.resize(mesh.inner.size());

//	print_function(stdout, &U[0], mesh, x, y, z);
//	fflush(stdout);

	Chafe chafe(mesh, tau, sigma, mu);

	for (i = 0; i < steps; ++i) {
		init_bnd(mesh, B, bnd, tau * (i + 1));
		chafe.solve(&U[0], &U[0], &B[0]);

		// check
		{
			init_func(mesh, Ans, ans, tau * (i + 1));
			fprintf(stderr, "time %lf/ norm %le\n", tau * (i + 1), nr2(&U[0], &Ans[0], U.size()));
//			vector_print(&U[0], U.size());
//			vector_print(&Ans[0], U.size());

			//mke_u2p(&P[0], &U[0], mesh);
			//vector_print(&P[0], P.size());
			//mke_u2p(&P[0], &Ans[0], mesh);
			//vector_print(&P[0], P.size());
		}

//		print_function(stdout, &F[0], mesh, x, y, z);
//		print_function(stdout, &Ans[0], mesh, x, y, z);
//		fflush(stdout);
	}

	return 0;
}

