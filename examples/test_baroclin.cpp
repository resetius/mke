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

void usage (const char * name)
{
	fprintf (stderr, "usage: %s [mesh.txt|-]\n", name);
	exit (1);
}

double f1 (double x, double y)
{
	return sin (x) * sin (y);
}

double f2 (double x, double y)
{
	return cos (x) * cos (y);
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

double rp_f (double x, double y, double t, double sigma, 
			 double mu, double sigma1,
			double mu1, double alpha, double theta)
{
	return -sigma* (2.*ipow (cos (x), 2) - 1.) *sin (x);
}

double rp_g (double x, double y, double t, double sigma, 
			 double mu, double sigma1,
					double mu1, double alpha, double theta)
{
	return -sigma* (2.*ipow (cos (x), 2) - 1.) *sin (x);
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

double rp_f1(double x, double y, double t, double sigma, 
			 double mu, double sigma1,
					double mu1, double alpha, double theta)
{
	return -45*mu*sin(y+t)*x-(9./2.)*sigma*
		ipow(cos(x),3)*sin(y+t)*sin(x)+(9./2.)*sigma*
		ipow(cos(x),3)*sin(x)*cos(y+t)-10*sigma*
		ipow(cos(x),4)*x*sin(y+t)+10*sigma*
		ipow(cos(x),4)*x*cos(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*sin(y+t)-(15./2.)*sigma*
		ipow(cos(x),2)*x*cos(y+t)-360*mu*sin(y+t)*sin(x)*
		ipow(cos(x),3)+147*mu*sin(y+t)*sin(x)*cos(x)-400*mu*sin(y+t)*x*
		ipow(cos(x),4)+390*mu*sin(y+t)*x*
		ipow(cos(x),2)-20*x*cos(y+t)*
		ipow(cos(x),4)-9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2);
}

double rp_g1(double x, double y, double t, double sigma, 
			 double mu, double sigma1,
					double mu1, double alpha, double theta)
{
	double alpha2 = alpha * alpha;
	double r = 20*x*sin(y+t)*
		ipow(cos(x),4)+9*sin(y+t)*
		ipow(cos(x),3)*sin(x)-15*x*sin(y+t)*
		ipow(cos(x),2)+390*mu*cos(y+t)*x*
		ipow(cos(x),2)-10*sigma*
		ipow(cos(x),4)*x*cos(y+t)-(9./2.)*sigma*
		ipow(cos(x),3)*sin(y+t)*sin(x)-alpha2*
		ipow(cos(x),7)*x-45*mu*cos(y+t)*x+18*
		ipow(cos(x),6)*
		ipow(cos(y+t),2)*sin(x)-9*
		ipow(cos(x),6)*sin(x)+9*x*
		ipow(cos(x),5)-(9./2.)*sigma*
		ipow(cos(x),3)*sin(x)*cos(y+t)-30*x*x*
		ipow(cos(x),4)*sin(x)-18*x*
		ipow(cos(x),5)*
		ipow(cos(y+t),2)+alpha2*
		ipow(cos(x),4)*x*sin(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*sin(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*cos(y+t)-400*mu*cos(y+t)*x*
		ipow(cos(x),4)+147*mu*cos(y+t)*sin(x)*cos(x)+60*x*x*
		ipow(cos(x),4)*
		ipow(cos(y+t),2)*sin(x)-360*mu*cos(y+t)*sin(x)*
		ipow(cos(x),3)-10*sigma*
		ipow(cos(x),4)*x*sin(y+t)+4*alpha2*
		ipow(cos(x),6)*x*x*sin(x)-9*alpha2*
		ipow(cos(x),3)*mu1*cos(y+t)*sin(x)-20*alpha2*
		ipow(cos(x),4)*mu1*cos(y+t)*x+15*alpha2*
		ipow(cos(x),2)*mu1*cos(y+t)*x-alpha2*
		ipow(cos(x),4)*sigma1*x*cos(y+t);
	r /= alpha2;
	return r;
}

double u1_t (double x, double y, double t)
{
	return x*sin(y+t)*ipow(cos(x),4);
}

double u2_t (double x, double y, double t)
{
	return x*cos(y+t)*ipow(cos(x),4);
}

void test_boclinic (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.001;
	double t = 0.0;
	int steps = 100000;
	double sigma  = 1.6e-2;
	double mu     = 8e-2;
	double sigma1 = sigma;
	double mu1    = mu;
	double alpha  = 1.0;

//	Baroclin bc (m, rp_f, rp_g, coriolis, tau, 
//		sigma, mu, sigma1, mu1, alpha);

	Baroclin bc (m, rp_f1, rp_g1, zero_coriolis, tau, 
		sigma, mu, sigma1, mu1, alpha);

	vector < double > u1 (sz);
	vector < double > u2 (sz);
	vector < double > bnd (std::max (os, 1) );

	vector < double > u1r (sz);
	vector < double > u2r (sz);

//	mke_proj (&u1[0], m, f1);
//	mke_proj (&u2[0], m, f2);

	mke_proj (&u1[0], m, u1_t, 0);
	mke_proj (&u2[0], m, u2_t, 0);

	setbuf (stdout, 0);
	for (int i = 0; i < steps; ++i)
	{
		bc.calc (&u1[0], &u2[0], &u1[0], &u2[0], &bnd[0], t);

		t += tau;

		fprintf (stderr, " === NORM1 = %le\n", bc.norm(&u1[0]));
		fprintf (stderr, " === NORM2 = %le\n", bc.norm(&u2[0]));

		{
			mke_proj (&u1r[0], m, u1_t, t);
			mke_proj (&u2r[0], m, u2_t, t);

			fprintf (stderr, " === DIST1 = %le\n", bc.dist(&u1[0],&u1r[0]));
			fprintf (stderr, " === DIST2 = %le\n", bc.dist(&u2[0],&u2r[0]));
		}

//		print_function (stdout, &u1r[0], m, x, y, z);
//		print_function (stdout, &u2r[0], m, x, y, z);

//		print_function (stdout, &u1[0], m, x, y, z);
		print_function (stdout, &u2[0], m, x, y, z);
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

	test_boclinic (mesh);
	return 0;
}

