/* -*- charset: utf-8 -*- */

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "baroclin.h"
#include "util.h"

using namespace phelm;
using std::string;
using std::vector;

void usage (const char * name)
{
	fprintf (stderr, "usage: %s [-f|--file mesh.txt|-] [-t|--threads number] [--task task] [--verbose|-v number]\n", name);
	fprintf (stderr, "tasks:\n"
	         "test\n"
	         "kornev1\n"
	        );

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
             double mu1, double alpha)
{
	return -sigma* (2.*ipow (cos (x), 2) - 1.) *sin (x);
}

double rp_g (double x, double y, double t, double sigma,
             double mu, double sigma1,
             double mu1, double alpha)
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

double rp_f1 (double x, double y, double t, double sigma,
              double mu, double sigma1,
              double mu1, double alpha)
{
	return -45*mu*sin (y + t) *x - (9. / 2.) *sigma*
	       ipow (cos (x), 3) *sin (y + t) *sin (x) + (9. / 2.) *sigma*
	       ipow (cos (x), 3) *sin (x) *cos (y + t) - 10*sigma*
	       ipow (cos (x), 4) *x*sin (y + t) + 10*sigma*
	       ipow (cos (x), 4) *x*cos (y + t) + (15. / 2.) *sigma*
	       ipow (cos (x), 2) *x*sin (y + t) - (15. / 2.) *sigma*
	       ipow (cos (x), 2) *x*cos (y + t) - 360*mu*sin (y + t) *sin (x) *
	       ipow (cos (x), 3) + 147*mu*sin (y + t) *sin (x) *cos (x) - 400*mu*sin (y + t) *x*
	       ipow (cos (x), 4) + 390*mu*sin (y + t) *x*
	       ipow (cos (x), 2) - 20*x*cos (y + t) *
	       ipow (cos (x), 4) - 9*cos (y + t) *
	       ipow (cos (x), 3) *sin (x) + 15*x*cos (y + t) *
	       ipow (cos (x), 2);
}

double rp_g1 (double x, double y, double t, double sigma,
              double mu, double sigma1,
              double mu1, double alpha)
{
	double alpha2 = alpha * alpha;
	double r = 20 * x * sin (y + t) *
	           ipow (cos (x), 4) + 9 * sin (y + t) *
	           ipow (cos (x), 3) * sin (x) - 15 * x * sin (y + t) *
	           ipow (cos (x), 2) + 390 * mu * cos (y + t) * x *
	           ipow (cos (x), 2) - 10 * sigma *
	           ipow (cos (x), 4) * x * cos (y + t) - (9. / 2.) * sigma *
	           ipow (cos (x), 3) * sin (y + t) * sin (x) - alpha2 *
	           ipow (cos (x), 7) * x - 45 * mu * cos (y + t) * x + 18 *
	           ipow (cos (x), 6) *
	           ipow (cos (y + t), 2) * sin (x) - 9 *
	           ipow (cos (x), 6) * sin (x) + 9 * x *
	           ipow (cos (x), 5) - (9. / 2.) * sigma *
	           ipow (cos (x), 3) * sin (x) * cos (y + t) - 30 * x * x *
	           ipow (cos (x), 4) * sin (x) - 18 * x *
	           ipow (cos (x), 5) *
	           ipow (cos (y + t), 2) + alpha2 *
	           ipow (cos (x), 4) * x * sin (y + t) + (15. / 2.) * sigma *
	           ipow (cos (x), 2) * x * sin (y + t) + (15. / 2.) * sigma *
	           ipow (cos (x), 2) * x * cos (y + t) - 400 * mu * cos (y + t) * x *
	           ipow (cos (x), 4) + 147 * mu * cos (y + t) * sin (x) * cos (x) + 60 * x * x *
	           ipow (cos (x), 4) *
	           ipow (cos (y + t), 2) * sin (x) - 360 * mu * cos (y + t) * sin (x) *
	           ipow (cos (x), 3) - 10 * sigma *
	           ipow (cos (x), 4) * x * sin (y + t) + 4 * alpha2 *
	           ipow (cos (x), 6) * x * x * sin (x) - 9 * alpha2 *
	           ipow (cos (x), 3) * mu1 * cos (y + t) * sin (x) - 20 * alpha2 *
	           ipow (cos (x), 4) * mu1 * cos (y + t) * x + 15 * alpha2 *
	           ipow (cos (x), 2) * mu1 * cos (y + t) * x - alpha2 *
	           ipow (cos (x), 4) * sigma1 * x * cos (y + t);
	r /= alpha2;
	return r;
}

/**
 * Операторы:
 * proc (u1, u2) options operator, arrow; diff(L(u1), t)+(1/2)*sigma*L(u1-u2)-mu*L(L(u1)) end proc
 * proc (u1, u2) options operator, arrow; diff(L(u2), t)+(1/2)*sigma*L(u1+u2)-mu*L(L(u2)) end proc
 */
double rp_f2 (double x, double y, double t, double sigma,
              double mu, double sigma1,
              double mu1, double alpha)
{
	return -9*cos (y + t) *
	       ipow (cos (x), 3) *sin (x) - 20*x*cos (y + t) *
	       ipow (cos (x), 4) + 15*x*cos (y + t) *
	       ipow (cos (x), 2) - (9. / 2.) *sigma*
	       ipow (cos (x), 3) *sin (y + t) *sin (x) + (9. / 2.) *sigma*
	       ipow (cos (x), 3) *sin (x) *cos (y + t) - 10*sigma*
	       ipow (cos (x), 4) *x*sin (y + t) + 10*sigma*
	       ipow (cos (x), 4) *x*cos (y + t) + (15. / 2.) *sigma*
	       ipow (cos (x), 2) *x*sin (y + t) - (15. / 2.) *sigma*
	       ipow (cos (x), 2) *x*cos (y + t) - 360*mu*sin (y + t) *sin (x) *
	       ipow (cos (x), 3) + 147*mu*sin (y + t) *sin (x) *cos (x) - 400*mu*sin (y + t) *x*
	       ipow (cos (x), 4) + 390*mu*sin (y + t) *x*
	       ipow (cos (x), 2) - 45*mu*sin (y + t) *x;
}

double rp_g2 (double x, double y, double t, double sigma,
              double mu, double sigma1,
              double mu1, double alpha)
{
	return -15*x*sin (y + t) *
	       ipow (cos (x), 2) + 20*x*sin (y + t) *
	       ipow (cos (x), 4) + 9*sin (y + t) *
	       ipow (cos (x), 3) *sin (x) - (9. / 2.) *sigma*
	       ipow (cos (x), 3) *sin (x) *cos (y + t) - 10*sigma*
	       ipow (cos (x), 4) *x*sin (y + t) - 10*sigma*
	       ipow (cos (x), 4) *x*cos (y + t) - (9. / 2.) *sigma*
	       ipow (cos (x), 3) *sin (y + t) *sin (x) + (15. / 2.) *sigma*
	       ipow (cos (x), 2) *x*sin (y + t) + (15. / 2.) *sigma*
	       ipow (cos (x), 2) *x*cos (y + t) - 360*mu*cos (y + t) *sin (x) *
	       ipow (cos (x), 3) + 147*mu*cos (y + t) *sin (x) *cos (x) - 400*mu*cos (y + t) *x*
	       ipow (cos (x), 4) + 390*mu*cos (y + t) *x*
	       ipow (cos (x), 2) - 45*mu*cos (y + t) *x;
}

double u1_t (double x, double y, double t)
{
	return x*sin (y + t) *ipow (cos (x), 4);
}

double u2_t (double x, double y, double t)
{
	return x*cos (y + t) *ipow (cos (x), 4);
}

void test_boclinic (const Mesh & m)
{
	int sz = (int) m.ps.size();
	int os = (int) m.outer.size();

	double tau = 0.001;
	double t = 0.0;
	int steps = 100000;
	double sigma  = 1.6e-2;
	double mu     = 8e-2;
	double sigma1 = sigma;
	double mu1    = mu;
	double alpha  = 1.0;
	//double alpha  = 0.0;

//	Baroclin bc (m, rp_f, rp_g, coriolis, tau,
//		sigma, mu, sigma1, mu1, alpha);

	Baroclin bc (m, rp_f1, rp_g1, zero_coriolis, tau,
	             sigma, mu, sigma1, mu1, alpha);
	bc.info();

//	Baroclin bc (m, rp_f2, rp_g2, zero_coriolis, tau,
//		sigma, mu, sigma1, mu1, alpha);

	vector < double > u1 (sz);
	vector < double > u2 (sz);
	vector < double > bnd (std::max (os, 1) );

	vector < double > u1r (sz);
	vector < double > u2r (sz);

//	proj (&u1[0], m, f1);
//	proj (&u2[0], m, f2);

	proj (&u1[0], m, u1_t, 0);
	proj (&u2[0], m, u2_t, 0);

	setbuf (stdout, 0);
	for (int i = 0; i < steps; ++i)
	{
		bc.calc (&u1[0], &u2[0], &u1[0], &u2[0], &bnd[0], t);

		t += tau;

		fprintf (stderr, " === NORM1 = %le\n", bc.norm (&u1[0]) );
		fprintf (stderr, " === NORM2 = %le\n", bc.norm (&u2[0]) );

		{
			proj (&u1r[0], m, u1_t, t);
			proj (&u2r[0], m, u2_t, t);

			fprintf (stderr, " === DIST1 = %le\n", bc.dist (&u1[0], &u1r[0]) );
			fprintf (stderr, " === DIST2 = %le\n", bc.dist (&u2[0], &u2r[0]) );
		}

//		print_function (stdout, &u1r[0], m, x, y, z);
//		print_function (stdout, &u2r[0], m, x, y, z);

//		print_function (stdout, &u1[0], m, x, y, z);
//		print_function (stdout, &u2[0], m, x, y, z);
	}
}

double kornev1_rp_ (double phi, double lambda)
{
	double omg   = 2.0 * M_PI / 24. / 60. / 60.;
	double T0    = 1. / omg;
	double H     = 5000;
	double sigma = 1. / 20. / 2. / M_PI;
	double R     = 6.371e+6;
	double x   = phi;

	double pt1 = -0.5 * (sin (x) * M_PI * x - 2 * sin (x) * x * x);
	if (fabs (pt1) > 1e-14)
	{
		pt1 /= cos (x);
	}

	double pt2 = -0.5 * (-M_PI + 4 * x);
	return -T0 / R * 16.0 / M_PI / M_PI * 30.0 * (pt1 + pt2);
}

double kornev1_rp1 (double phi, double lambda,
                    double t, double sigma,
                    double mu, double sigma1,
                    double mu1, double alpha)
{
	return kornev1_rp_ (phi, lambda);
}

double kornev1_rp2 (double phi, double lambda, double t,
                    double sigma,
                    double mu, double sigma1,
                    double mu1, double alpha)
{
	return 0.1 * kornev1_rp_ (phi, lambda);
}

double kornev1_coriolis (double phi, double lambda)
{
	double omg  = 2.0 * M_PI / 24. / 60. / 60.;
	double H    = 1;
	return 2.*sin (phi) + // l
	       0.5 * cos (2*lambda) *sin (2*phi) *sin (2*phi); //h
}

double kornev1_u1 (double phi, double lambda)
{
	double omg = 2.*M_PI / 24. / 60. / 60.; // ?
	double T0  = 1. / omg;
	double R   = 6.371e+6;

	return -T0 / R * 16.0 / M_PI / M_PI * 30.0 *
	       (M_PI / 4 * phi * phi - phi * phi * phi / 3);
}

double kornev1_u2 (double phi, double lambda)
{
	return 0.1 * kornev1_u1 (phi, lambda);
}

void test_kornev1 (const Mesh & m)
{
	int sz = m.size;
	int os = m.outer_size;

	double days = 30;
	double T = days * 2.0 * M_PI;

	double tau = 0.001;
	double t = 0.0;
	int steps = 100000;
	double sigma  = 1.14e-2;
	double mu     = 6.77e-5;
	double sigma1 = sigma;
	double mu1    = mu;
	double alpha  = 1.0;
	double nr, elapsed;
	int i = 0;

	Baroclin bc (m, kornev1_rp1, kornev1_rp2, kornev1_coriolis, tau,
	             sigma, mu, sigma1, mu1, alpha);

	vector < double > u1 (sz);
	vector < double > u2 (sz);
	vector < double > bnd (std::max (os, 1) );

	proj (&u1[0], m, kornev1_u1);
	proj (&u2[0], m, kornev1_u2);

	setbuf (stdout, 0);

	while (t < T)
	{
		Timer tm;
		bc.calc (&u1[0], &u2[0], &u1[0], &u2[0], &bnd[0], t);

		elapsed = tm.elapsed();

		nr = bc.norm (&u1[0]);

		fprintf (stderr, "t=%le; nr=%le; min=%le; max=%le; work=%le;\n",
		         t, nr, vec_find_min (&u1[0], sz),
		         vec_find_max (&u1[0], sz), tm.elapsed() );

		print_function (stdout, &u1[0], m, x, y, z);

		if (isnan (nr) || isinf (nr) )
		{
			return;
		}

		nr = bc.norm (&u2[0]);

		fprintf (stderr, "t=%le; nr=%le; min=%le; max=%le; work=%le;\n",
		         t, nr, vec_find_min (&u2[0], sz),
		         vec_find_max (&u2[0], sz), tm.elapsed() );

		print_function (stdout, &u2[0], m, x, y, z);

		if (isnan (nr) || isinf (nr) )
		{
			return;
		}

		i += 1;
		t += tau;
	}
}

int main (int argc, char *argv[])
{
	Mesh mesh;
	string task;

	for (int i = 0; i < argc; ++i)
	{
		if (!strcmp (argv[i], "--help") || !strcmp (argv[i], "-h") )
		{
			usage (argv[0]);
		}
		else if (!strcmp (argv[i], "--file") || !strcmp (argv[i], "-f") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			FILE * f = (strcmp (argv[i + 1], "-") == 0) ? stdin : fopen (argv[i + 1], "rb");

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
		else if (!strcmp (argv[i], "--threads") || !strcmp (argv[i], "-t") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			int threads = atoi (argv[i + 1]);
			set_num_threads (threads);
		}
		else if (!strcmp (argv[i], "--task") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}
			task = argv[i + 1];
		}
	}

#if defined(WIN32)
	set_fpe_except();
#endif

	mesh.info();

	if (task == "test")
	{
		test_boclinic (mesh);
	}
	else if (task == "kornev1")
	{
		test_kornev1 (mesh);
	}
	else
	{
		usage (argv[0]);
	}
	return 0;
}

