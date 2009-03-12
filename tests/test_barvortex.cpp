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

#include "barvortex.h"
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

double an(double x, double y)
{
	return (-sin(x)*cos(y)*sin(x)*cos(y)+cos(x)*sin(y)*cos(x)*sin(y))/cos(x);
	//return (-sin(x)*cos(y)*sin(x)*cos(y)) / cos(x);
	//return (-cos(x)*sin(y)*cos(x)*sin(y))/cos(x);

	//return cos(x) * sin(y) / cos(x);	
	//return sin(x) * cos(y) / cos(x);
}

void test_jacobian(const Mesh & m)
{
	Jacobian j(m);
	vector < double > F1;
	vector < double > F2;
	vector < double > ans1;
	vector < double > rans1;
	vector < double > bnd;

	mke_proj(m, F1, f1);
	mke_proj(m, F2, f2);
	mke_proj(m, rans1, an);
	mke_proj_bnd(m, bnd, an);

	ans1.resize(m.ps.size());
	j.calc1(&ans1[0], &F1[0], &F2[0], &bnd[0]);

	fprintf(stderr, "jacobian  err=%.2le\n", 
		mke_dist(&ans1[0], &rans1[0], m, sphere_scalar_cb));

	//vector < double > p1(m.inner.size());
	//mke_u2p(&p1[0], &rans1[0], m);
	//vector_print(&p1[0], p1.size());
	//mke_u2p(&p1[0], &ans1[0], m);
	//vector_print(&p1[0], p1.size());

	//vector_print(&rans1[0], rans1.size());
	//vector_print(&ans1[0], ans1.size());
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

	test_jacobian(mesh);
	return 0;
}

