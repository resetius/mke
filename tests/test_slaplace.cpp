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
	//return -6.0 * sin(y) * sin(2.0 * x);
	
	double sx = sin(x);
	double cx = cos(x);
	double sy = sin(y);
	double cy = cos(y);

	double ans = 0.0;

	// 0.5 * ipow(cx * cy - 1.0, 2)
	ans += (sx * sx * (cx * cy - 1) * cy +
		cx * sx * sx * cy * cy -
		cx * cx * (cx * cy - 1) * cy) / cx;
	ans += sy * sy - (cx * cy - 1) * cy / cx;

	// ipow(cx * cy - sx, 2)
	ans += -2 * sx * (cx * cy - sx) * (-sx * cy - cx) / cx +
		2 * ipow(-sx * cy - cx, 2) +
		2 * (cx * cy - sx) * (-cx * cy + sx);
	       
	ans += 2 * sy * sy - 2 * (cx * cy - sx) * cy / cx;

	// 2.0 * ipow(sx - 3, 2)
	ans += -8 * (sx - 3.0) * sx + 4.0 * cx * cx;
	return ans;
}

double ans(double x, double y)
{
	//return sin(y) * sin(2.0 * x);

	double sx = sin(x);
	double cx = cos(x);
	//double sy = sin(y);
	double cy = cos(y);

	return 0.5 * ipow(cx * cy - 1.0, 2) +
		ipow(cx * cy - sx, 2) +
		2.0 * ipow(sx - 3, 2);
}

double bnd(double x, double y)
{
	return ans(x, y);
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

void test_invert(const Mesh & mesh)
{
	int sz = mesh.ps.size();
	int rs = mesh.inner.size();
	int os = mesh.outer.size();

	vector < double > F(sz);
	vector < double > B(os);
	vector < double > Ans(sz);
	vector < double > rans(sz);

	mke_proj(&F[0], mesh, rp);
	mke_proj(&rans[0], mesh, ans);
	mke_proj_bnd(&B[0], mesh, bnd);

	SphereLaplace l(mesh);
	l.solve(&Ans[0], &F[0], &B[0]);

	fprintf(stderr, "err=%.2le\n", mke_dist(&Ans[0], &rans[0], mesh, sphere_scalar_cb));

	{
		FILE * f = fopen("invert_answer.txt", "wb");
		print_function(f, &Ans[0], mesh, x, y, z);
		fclose(f);

		fprintf(stderr, "invert answer saved to 'invert_answer.txt'\n");

		f = fopen("answer.txt", "wb");
		print_function(f, &rans[0], mesh, x, y, z);
		fclose(f);

		f = fopen("rp.txt", "wb");
		print_function(f, &F[0], mesh, x, y, z);
		fclose(f);
	}
}

void test_laplace(const Mesh & mesh)
{
	int sz = mesh.ps.size();
	int rs = mesh.inner.size();
	int os = mesh.outer.size();

	vector < double > U(sz);
	vector < double > LU(sz);
	vector < double > LU1(sz);
	vector < double > B(os);
	vector < double > P(rs);
	vector < double > P1(rs);

	mke_proj(&U[0], mesh, ans);
	mke_proj(&LU[0], mesh, rp);
	mke_proj_bnd(&B[0], mesh, rp);

	SphereLaplace l(mesh);
	l.calc1(&LU1[0], &U[0], &B[0]);

	fprintf(stderr, "laplace err=%.2le\n", mke_dist(&LU[0], &LU1[0], mesh));

	mke_proj_bnd(&B[0], mesh, ans);
	l.solve(&LU[0], &LU1[0], &B[0]);

	fprintf(stderr, "laplace err=%.2le\n", mke_dist(&U[0], &LU[0], mesh));
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
