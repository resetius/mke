/* -*- charset: utf-8 -*- */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vector>

#include "util.h"
#include "laplace.h"
#include "norm.h"

using namespace std;
using namespace phelm;
using namespace linal;

static void usage (const char * name)
{
	fprintf (stderr, "usage: %s [--order order(1|2|3, before -f)] [-f|--file mesh.txt|-] [-d|--double] [-t|--threads number]\n", name);
	exit (1);
}

static double func_lapl (double x, double y)
{
	//return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
	return 2.0 * y * y - 2.0 * y + 2.0 * x * x - 2.0 * x;// + 2;
}

static double func (double x, double y)
{
	//return sin(M_PI * x) * sin(M_PI * y) + 10.0;
	return x * (x - 1) * y * (y - 1);// + x * x;
}

static double bnd (double x, double y)
{
	return func (x, y);
	//return 0.0;
	//return 1.0;
}

static double nr2 (double * a, double * b, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; ++i)
	{
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt (sum);
}

template < typename T >
void test_invert (Mesh & mesh)
{
	typedef Array < T, Allocator < T > > vector;

	int sz = (int) mesh.ps.size();
	int rs = (int) mesh.inner.size();
	int os = (int) mesh.outer.size();
	FlatNorm < T > nr (mesh);

	std::vector < T > F (sz);
	std::vector < T > B (os);
	std::vector < T > func1 (sz);
	std::vector < T > rans (sz);

	vector cF (sz);
	vector cB (os);
	vector cAns (sz);
	vector crans (sz);

	proj (&F[0], mesh, func_lapl);
	proj (&rans[0], mesh, func);
	proj_bnd (&B[0], mesh, bnd);

	vec_copy_from_host (&cF[0], &F[0], sz);
	vec_copy_from_host (&crans[0], &rans[0], sz);
	vec_copy_from_host (&cB[0], &B[0], os);

	Timer t;
	Laplace < T > l (mesh);
	fprintf (stderr, "l -> %lf\n", t.elapsed() );
	l.solve (&cAns[0], &cF[0], &cB[0]);

	vec_copy_from_device (&func1[0], &cAns[0], sz);

	{
		FILE * f = fopen ("lu_1_real.txt", "w");
		print_function (f, &rans[0], mesh);
		fclose (f);
		f = fopen ("lu_1_calc.txt", "w");
		print_function (f, &func1[0], mesh);
		fclose (f);
	}

	fprintf (stdout, "L1: invert  err=%.2le\n",
	         nr.dist (&cAns[0], &crans[0]) );
}

template < typename T >
void test_laplace (Mesh & mesh)
{
	typedef Array < T, Allocator < T > > vector;
	int sz = (int) mesh.ps.size();
	int rs = (int) mesh.inner.size();
	int os = (int) mesh.outer.size();
	FlatNorm < T > nr (mesh);
	Timer t;

	std::vector < T > U (sz);
	std::vector < T > LU (sz);
	std::vector < T > LU1 (sz);
	std::vector < T > B1 (os);
	std::vector < T > B2 (os);

	vector cU (sz);
	vector cLU (sz);
	vector cLU1 (sz);
	vector cB1 (os);
	vector cB2 (os);

	fprintf (stderr, "prepare ... \n");

	proj (&U[0], mesh, func);
	proj (&LU[0], mesh, func_lapl);
	proj_bnd (&B1[0], mesh, func_lapl);
	proj_bnd (&B2[0], mesh, func);

	vec_copy_from_host (&cU[0], &U[0], sz);
	vec_copy_from_host (&cLU[0], &LU[0], sz);
	vec_copy_from_host (&cB1[0], &B1[0], os);
	vec_copy_from_host (&cB2[0], &B2[0], os);

	fprintf (stderr, "start ... \n");

	Laplace < T > l (mesh);

	t.restart();
	l.calc1 (&cLU1[0], &cU[0], &cB1[0]);
	fprintf (stderr, "calc1 time: %lf\n", t.elapsed() );

	vec_copy_from_device (&LU1[0], &cLU1[0], sz);

	fprintf (stdout, "L2: laplace err=%.2le\n",
	         nr.dist (&cLU[0], &cLU1[0]) );
	{
		FILE * f = fopen ("lu_real.txt", "w");
		print_inner_function (f, &LU[0], mesh);
		fclose (f);
		f = fopen ("lu_calc.txt", "w");
		print_inner_function (f, &LU1[0], mesh);
		fclose (f);

		f = fopen ("lu_diff.txt", "w");
		vector tmp1 (sz);
		std::vector < T > tmp2 (sz);

		vec_diff (&tmp1[0], &cLU[0], &cLU1[0], sz);
		vec_copy_from_device (&tmp2[0], &tmp1[0], sz);

		print_inner_function (f, &tmp2[0], mesh);
		fclose (f);
	}

	t.restart();
	l.solve (&cLU[0], &cLU1[0], &cB2[0]);
	fprintf (stderr, "solve time: %lf\n", t.elapsed() );

	vec_copy_from_device (&LU[0], &cLU[0], sz);

//	vector_print(&U[0], U.size());
//	vector_print(&LU[0], LU.size());

	fprintf (stdout, "L3: laplace err=%.2le\n", nr.dist (&cU[0], &cLU[0]) );
}

int test_laplace (int argc, char *argv[])
{
	Mesh mesh;
	bool use_double = false;
	int order = 1;

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
			if (!mesh.load (f, order) )
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
		else if (!strcmp (argv[i], "--double") || !strcmp (argv[i], "-d") )
		{
			use_double = true;
		}
		else if (!strcmp(argv[i], "--order"))
		{
			if (i == argc - 1)
			{
				usage(argv[0]);
			}

			order = atoi(argv[i + 1]);
		}
	}

	if (mesh.ps.empty() )
	{
		usage (argv[0]);
	}

	mesh.info();

	linal_init();

	Timer t;
	try
	{
		if (check_device_supports_double() && use_double)
		{
			fprintf (stderr, "using double\n");
			test_invert < double > (mesh);
			test_laplace < double > (mesh);
		}
		else
		{
			fprintf (stderr, "using float\n");
			test_invert < float > (mesh);
			test_laplace < float > (mesh);
		}
	}
	catch (const std::exception & e)
	{
		fprintf (stderr, "exception: %s\n", e.what() );
	}

	fprintf (stderr, "elapsed: %lf\n", t.elapsed() );

	linal_shutdown();
	return 0;
}

