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
	//return -16.0 * sin(y) * sin(4.0 * x) 
	//	- sin(y) * sin(4.0 * x) / cos(x) / cos(x);

	return -6.0 * sin(y) * sin(2.0 * x);

	//return -1.0 / cos(x) / cos(x) * sin(y) * sin(2.0 * x) -
	//	2.0 * sin(x) / cos(x) * sin(y) * cos(2.0 * x) -
	//	4.0 * sin(y) * sin(2.0 * x);

	//return -1.0 / cos(x) / cos(x) * sin(y) * sin(4.0 * x) -
	//	4.0 * sin(x) / cos(x) * sin(y) * cos(4.0 * x) -
	//	16.0 * sin(y) * sin(4.0 * x);
}

double ans(double x, double y)
{
	return sin(y) * sin(2.0 * x);
	//return sin(y) * sin(4.0 * x);
	//return 0.0;
}

double bnd(double x, double y)
{
	return ans(x, y);
	//return 0.0;
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

static double x(double u, double v)
{
	return cos(u) * cos(v);
}

static double y(double u, double v)
{
	return cos(u) * sin(v);
}

static double z(double u, double v)
{
	return sin(u);
}

int main(int argc, char *argv[])
{
	Mesh mesh;
	vector < double > F;
	vector < double > B;
	vector < double > Ans;
	vector < double > rans;

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

	init_func(mesh, F, rp);
	init_func(mesh, rans, ans);
	init_bnd(mesh, B, bnd);
	Ans.resize(F.size());
	sphere_laplace_solve(&Ans[0], mesh, &F[0], &B[0]);
//	fprintf(stderr, "answer: 1\n");
//	vector_print(&Ans[0], Ans.size());
//	fprintf(stderr, "real:   2\n");
//	vector_print(&rans[0], rans.size());
	fprintf(stderr, "err=%.2le\n", nr2(&Ans[0], &rans[0], rans.size()));

	{
		FILE * f = fopen("answer.txt", "wb");
		print_function(f, &Ans[0], mesh, x, y, z);
		fclose(f);

		fprintf(stderr, "answer saved to 'answer.txt'\n");
	}
	return 0;
}
