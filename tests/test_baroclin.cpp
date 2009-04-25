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
#include <math.h>
#include <stdlib.h>

#include "baroclin.h"
#include "util.h"

using namespace std;

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [mesh.txt|-]\n", name);
	exit(1);
}

double f1(double x, double y)
{
	return sin(x) * sin(y);
}

double f2(double x, double y)
{
	return cos(x) * cos(y);
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

void test_boclinic(const Mesh & m)
{
	int sz = m.ps.size();
	int os = m.outer.size();

	double tau = 0.001;
	int steps = 100000;
	Baroclin bc(m, tau, 1.6e-2, 8e-5);

	vector < double > u1(sz);
	vector < double > u2(sz);
	vector < double > bnd(std::max(os, 1));

	mke_proj(&u1[0], m, f1);
	mke_proj(&u2[0], m, f2);

	setbuf(stdout, 0);
	for (int i = 0; i < steps; ++i) {
		bc.calc(&u1[0], &u2[0], &u1[0], &u2[0], &bnd[0], (double)i * tau);

		fprintf(stderr, " === NORM1 = %le\n",
			mke_norm(&u1[0], m, sphere_scalar_cb));
		fprintf(stderr, " === NORM2 = %le\n",
			mke_norm(&u2[0], m, sphere_scalar_cb));

		print_function(stdout, &u1[0], m, x, y, z);
		print_function(stdout, &u2[0], m, x, y, z);
	}
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
	set_fpe_except();

	test_boclinic(mesh);
	return 0;
}

