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

double rp(double x, double y)
{
	//return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
	return 2.0 * y * y - 2.0 * y + 2.0 * x * x - 2.0 * x;
}

double ans(double x, double y)
{
	//return sin(M_PI * x) * sin(M_PI * y);
	return x * (x - 1) * y * (y - 1) + 1.0;
}

double bnd(double x, double y)
{
	//return 0.0;
	return 1.0;
}

template < typename T >
void init_bnd(Mesh & m, vector < double > & F, T f)
{
	F.resize(m.outer.size());
	for (size_t i = 0; i < m.outer.size(); ++i) {
		int p0 = m.outer[i];
		Point & p = m.ps[p0];
		F[i] = f(p.x, p.y);
	}
}

template < typename T >
void init_func(Mesh & mesh, vector < double > & F, T f)
{
	F.resize(mesh.ps.size());
	for (size_t i = 0; i < mesh.ps.size(); ++i)
	{
		F[i] = f(mesh.ps[i].x, mesh.ps[i].y);
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

void test_invert(Mesh & mesh)
{
	vector < double > F;
	vector < double > B;
	vector < double > Ans;
	vector < double > rans;

	init_func(mesh, F, rp);
	init_func(mesh, rans, ans);
	init_bnd(mesh, B, bnd);
	Ans.resize(F.size());
	laplace_solve(&Ans[0], mesh, &F[0], &B[0]);
//	fprintf(stderr, "1\n");
//	vector_print(&Ans[0], Ans.size());
//	fprintf(stderr, "2\n");
//	vector_print(&rans[0], rans.size());
	fprintf(stderr, "invert  err=%.2le\n", nr2(&Ans[0], &rans[0], rans.size()));
}

void test_laplace(Mesh & mesh)
{
	vector < double > U;
	vector < double > LU;
	vector < double > LU1;

	init_func(mesh, U, ans);
	init_func(mesh, LU, rp);
	init_func(mesh, LU1, rp);

	laplace_calc(&LU1[0], &LU[0], mesh);

	fprintf(stderr, "laplace err=%.2le\n", nr2(&LU[0], &LU1[0], LU.size()));
}

int main(int argc, char *argv[])
{
	Mesh mesh;
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

	test_invert(mesh);
	test_laplace(mesh);
	return 0;
}

