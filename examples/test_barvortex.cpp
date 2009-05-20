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
#include <time.h>

#include "barvortex.h"
#include "util.h"

using namespace std;

void usage (const char * name)
{
	fprintf (stderr, "usage: %s [mesh.txt|-]\n", name);
	exit (1);
}

double f1 (double x, double y)
{
	return sin (x) * sin (y) * cos(x);
}

double f2 (double x, double y)
{
	return 1000.0 * cos (x) * cos (y) * sin(x);
}

double f3 (double x, double y)
{
	return -(-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) ) / cos (x);
}

double an (double x, double y)
{
	return (-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) ) / cos (x);
	//return (-sin(x)*cos(y)*sin(x)*cos(y)) / cos(x);
	//return (-cos(x)*sin(y)*cos(x)*sin(y))/cos(x);

	//return cos(x) * sin(y) / cos(x);
	//return sin(x) * cos(y) / cos(x);
}

void test_jacobian (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int rs = (int)m.inner.size();
	int os = (int)m.outer.size();

	Jacobian j (m);
	vector < double > F1 (sz);
	vector < double > F2 (sz);
	vector < double > ans1 (sz);
	vector < double > rans1 (sz);
	vector < double > bnd (os);

	mke_proj (&F1[0], m, f1);
	mke_proj (&F2[0], m, f2);
	mke_proj (&rans1[0], m, an);
	mke_proj_bnd (&bnd[0], m, an);

	j.calc1 (&ans1[0], &F1[0], &F2[0], &bnd[0]);

	fprintf (stderr, "jacobian  err=%.2le\n",
	         mke_dist (&ans1[0], &rans1[0], m, sphere_scalar_cb) );

	//vector < double > p1(m.inner.size());
	//mke_u2p(&p1[0], &rans1[0], m);
	//vector_print(&p1[0], p1.size());
	//mke_u2p(&p1[0], &ans1[0], m);
	//vector_print(&p1[0], p1.size());

	//vector_print(&rans1[0], rans1.size());
	//vector_print(&ans1[0], ans1.size());
}

void rand_init(double * h, int n)
{
	for (int i = 0; i < n; ++i) {
		h[i] = (double)rand() / (double)RAND_MAX;
	}
}

void test_jacobian_T (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int rs = (int)m.inner.size();
	int os = (int)m.outer.size();
	double nr1, nr2;

	Jacobian j (m);
	vector < double > u (sz);
	vector < double > w (sz);
	vector < double > v (sz);
	vector < double > tmp(rs);
	vector < double > ans1 (sz);
	vector < double > ans2 (sz);

	srand((unsigned int)time(0));
	//mke_proj (&u[0], m, f1);
	//mke_proj (&w[0], m, f2);
	//mke_proj (&v[0], m, f3);
	rand_init(&u[0], (int)u.size()); mke_u2p(&tmp[0], &u[0], m); mke_p2u(&u[0], &tmp[0], 0, m); 
	rand_init(&v[0], (int)v.size()); mke_u2p(&tmp[0], &v[0], m); mke_p2u(&v[0], &tmp[0], 0, m); 
	rand_init(&w[0], (int)w.size()); mke_u2p(&tmp[0], &w[0], m); mke_p2u(&w[0], &tmp[0], 0, m); 

	j.calc1 (&ans1[0], &v[0], &w[0], 0);
	j.calc1t(&ans2[0], &u[0], &w[0], 0);

	nr1 = mke_scalar(&ans1[0], &u[0], m, sphere_scalar_cb);
	nr2 = mke_scalar(&ans2[0], &v[0], m, sphere_scalar_cb);

	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));

	j.calc1 (&ans1[0], &w[0], &v[0], 0);
	j.calc1t(&ans2[0], &w[0], &u[0], 0);

	nr1 = mke_scalar(&ans1[0], &u[0], m, sphere_scalar_cb);
	nr2 = mke_scalar(&ans2[0], &v[0], m, sphere_scalar_cb);

	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));
}

static double x (double u, double v)
{
	return cos (u) * cos (v);
}

static double y (double u, double v)
{
	return cos (u) * sin (v);
}

static double z (double u, double v)
{
	return sin (u);
}

double u0 (double x, double y)
{
	//return 0.0;
	return sin (x) * sin (y);
}

double z0 (double x, double y)
{
	return sin(x) * cos(x) / 25.0;
}

double rp (double x, double y, double t, double mu, double sigma)
{
	// double b = -6.0 * exp(t) * sin(y) * sin(2.0 * x);
	// return a - mu * b + sigma * a;
	// return 0.0;

	// lapl from -3.5 * sigma * ipow (sin (x), 3);
	return sigma*(2.*ipow(cos(x),2)-1.)*sin(x);
}

double zero_coriolis (double phi, double lambda)
{
	return 0.0;
}

double coriolis (double phi, double lambda)
{
	double omg = 0.0000727000000000;
	double l = omg * 2.0 * sin (phi);
	double h = cos (2.0 * lambda) * ipow (sin (2.0 * phi), 2);
	return l + h;
}

double an1 (double x, double y, double t)
{
	return x*sin(y+t)*ipow(cos(x),4);
}

double omg_an1 (double x, double y, double t)
{
	return -sin(y+t)*ipow(cos(x),2)*(9*sin(x)*cos(x)+20*x*ipow(cos(x),2)-15*x);
}

double rp1(double x, double y, double t, double mu, double sigma)
{
	return 390*mu*sin(y+t)*x*ipow(cos(x),2)
		+147*mu*sin(y+t)*sin(x)*cos(x)-
		400*mu*sin(y+t)*x*
		ipow(cos(x),4)-360*mu*sin(y+t)*
		sin(x)*ipow(cos(x),3)-20*sigma*sin(y+t)*
		ipow(cos(x),4)*x+15*sigma*sin(y+t)*
		ipow(cos(x),2)*x-9*sigma*sin(y+t)*
		ipow(cos(x),3)*sin(x)+30*cos(y+t)*
		ipow(cos(x),4)*sin(y+t)*x*x*sin(x)+9*cos(y+t)*
		ipow(cos(x),6)*sin(y+t)*sin(x)-45*mu*sin(y+t)*x-
		9*cos(y+t)*ipow(cos(x),5)*sin(y+t)*x-20*x*cos(y+t)*
		ipow(cos(x),4)-9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2);
}

double rp2(double x, double y, double t, double mu, double sigma)
{
	return -9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2)-20*x*cos(y+t)*
		ipow(cos(x),4)-9*sigma*sin(y+t)*
		ipow(cos(x),3)*sin(x)+15*sigma*sin(y+t)*
		ipow(cos(x),2)*x-20*sigma*sin(y+t)*
		ipow(cos(x),4)*x-360*mu*sin(y+t)*sin(x)*
		ipow(cos(x),3)+390*mu*sin(y+t)*x*
		ipow(cos(x),2)-400*mu*sin(y+t)*x*
		ipow(cos(x),4)+147*mu*sin(y+t)*sin(x)*cos(x)-45*mu*sin(y+t)*x;
}

#ifdef WIN32
#include <windows.h>
#undef max
#endif

void test_barvortex_L (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.05;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 8e-5;   //8e-5;
	double sigma = 1.6e-2; //1.6e-2;

	//BarVortex bv (m, rp1, zero_coriolis, tau, sigma, mu);
	BarVortex bv (m, rp, coriolis, tau, sigma, mu);

	vector < double > u (sz);
	vector < double > z (sz);
	vector < double > bnd (std::max (os, 1));
	vector < double > Ans(sz);
	vector < double > Ans1(sz);

	//mke_proj (&u[0], m, an1, 0);
	mke_proj (&u[0], m, u0);

	mke_proj (&z[0], m, z0);

	//if (!bnd.empty()) mke_proj_bnd(&bnd[0], m, f1);

	setbuf (stdout, 0);

	bv.L_step(&Ans[0], &u[0], &z[0]);
	bv.L_1_step(&Ans1[0], &Ans[0], &z[0]);

	fprintf (stderr, " ||L(L^1)|| = %le \n",
		mke_dist (&u[0], &Ans1[0], m, sphere_scalar_cb));
}

void test_barvortex_LT (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.05;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 8e-5;   //8e-5;
	double sigma = 1.6e-2; //1.6e-2;

	//BarVortex bv (m, rp1, zero_coriolis, tau, sigma, mu);
	BarVortex bv (m, rp, coriolis, tau, sigma, mu);

	vector < double > u  (sz);
	vector < double > v  (sz);
	vector < double > lv (sz);
	vector < double > ltu(sz);

	vector < double > z (sz);
	vector < double > bnd (std::max (os, 1));

	mke_proj (&u[0], m, f1);
	mke_proj (&v[0], m, f2);

	mke_proj (&z[0], m, z0);

	setbuf (stdout, 0);

	bv.L_step (&lv[0],  &v[0], &z[0]);
	bv.LT_step(&ltu[0], &u[0], &z[0]);

	double nr1 = mke_scalar(&lv[0],  &u[0], m, sphere_scalar_cb);
	double nr2 = mke_scalar(&ltu[0], &v[0], m, sphere_scalar_cb);

	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));
}

void test_barvortex (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.05;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 8e-5;   //8e-5;
	double sigma = 1.6e-2; //1.6e-2;

	//BarVortex bv (m, rp1, zero_coriolis, tau, sigma, mu);
	BarVortex bv (m, rp, coriolis, tau, sigma, mu);

	vector < double > u (sz);
	vector < double > bnd (std::max (os, 1));
	vector < double > Ans(sz);

	//mke_proj (&u[0], m, an1, 0);
	mke_proj (&u[0], m, u0);

	//if (!bnd.empty()) mke_proj_bnd(&bnd[0], m, f1);

	setbuf (stdout, 0);

	while (t < T)
	{
#if 1
		Timer tm;
		bv.calc (&u[0], &u[0], &bnd[0], t);
		if (i % 1 == 0) {
			fprintf (stderr, " === NORM = %le, STEP %lf of %lf: %lf\n",
			         mke_norm (&u[0], m, sphere_scalar_cb), t, T, tm.elapsed());
			// 3d print
			print_function (stdout, &u[0], m, x, y, z);
			// flat print
			// print_function (stdout, &u[0], m, 0, 0, 0);
		}
#endif

		i += 1;
		t += tau;
#if 0
		{
			mke_proj(&Ans[0], m, an1, t);
			fprintf(stderr, "time %lf/ norm %le\n", t, 
				mke_dist(&u[0], &Ans[0], m, sphere_scalar_cb));
//			print_function (stdout, &Ans[0], m, x, y, z);
		}
//		Sleep(500);
#endif
	}
}

int main (int argc, char *argv[])
{
	Mesh mesh;

	if (argc > 1)
	{
		FILE * f = (strcmp (argv[1], "-") == 0) ? stdin : fopen (argv[1], "rb");
		if (!f)
		{
			usage (argv[0]);
		}
		if (!mesh.load (f) )
		{
			usage (argv[0]);
		}
		fclose (f);
	}
	else
	{
		usage (argv[0]);
	}
	set_fpe_except();

	mesh.info();

	//test_jacobian(mesh);
	//test_jacobian_T(mesh);
	test_barvortex_L(mesh);
	//test_barvortex_LT(mesh);
	//test_barvortex (mesh);
	return 0;
}

