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
	int i, steps = 10000;
	double tau   = 0.001; //r5 => h = 0.03
	double mu    = 1.0;
	double sigma = -70;

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

	int sz = (int)mesh.ps.size();
	int os = (int)mesh.outer.size();
	int rs = (int)mesh.outer.size();

	vector < double > U(sz);
	vector < double > B(os);
	vector < double > Ans(sz);
	vector < double > P(rs);

	mke_proj(&U[0], mesh, ans, 0.0);

//	print_function(stdout, &U[0], mesh, x, y, z);
//	fflush(stdout);

	Chafe chafe(mesh, tau, sigma, mu);

	for (i = 0; i < steps; ++i) {
		mke_proj_bnd(&B[0], mesh, bnd, tau * (i + 1));
		chafe.solve(&U[0], &U[0], &B[0],  tau * (i));

		// check
		{
			mke_proj(&Ans[0], mesh, ans, tau * (i + 1));
			fprintf(stderr, "time %lf/ norm %le\n", tau * (i + 1), 
				mke_dist(&U[0], &Ans[0], mesh));
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

