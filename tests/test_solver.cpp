/* -*- charset: utf-8 -*- */

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "gmres.h"
#include "util.h"
#include "solver.h"
#include "phelm.h"
#include "laplace.h"
#include "slaplace.h"

using namespace phelm;
using namespace std;

void usage (const char * n)
{
	fprintf (stderr, "%s [--threads|-t n] [--double|-d] \
[--task] [--dim] [--iters] [--mesh] [--operator]\n", n);
	fprintf (stderr, "--threads n - sets the number of threads\n");
	fprintf (stderr, "--double - use double or float ? \n");
	fprintf (stderr, "--mesh - mesh for mult_sparse \n");
	fprintf (stderr, "--operator - operator for mult_sparse \n");
	fprintf (stderr, "--task - all\n"
	         "mult_dense\n"
	         "mult_sparse\n"
	         "invert_sparse\n");
	exit (-1);
}

bool check (float val)
{
	return fabs (val) < 1e-5;
}

bool check (double val)
{
	return fabs (val) < 1e-12;
}

template < typename Matrix >
void init_matrix2 (Matrix & m, int n)
{
	int i, j = 0;

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			m.add (i, j, 1. / (double) (i + j + 1) );
		}
	}
}

template < typename Matrix >
void init_matrix (Matrix & m, int n)
{
	int i, j = 0;

	/**
	 * 4 2 0 .... 0
	 * 1 4 2 .... 0
	 * 0 1 4 2 .. 0
	 * 0 0 1 4 .. 0
	 * ............
	 * ...... 1 4 2
	 * .........1 4
	 */

	for (i = 0; i < n; ++i)
	{
		if (i == 0)
		{
			/* 4 2 0 .... 0 */
			m.add (i, i,     4);
			m.add (i, i + 1, 2);
		}
		else if (i < n - 1)
		{
			/* 1 4 2 .... 0 */
			m.add (i, i - 1, 1);
			m.add (i, i,     4);
			m.add (i, i + 1, 2);
		}
		else
		{
			/* 0 .... 1 4 */
			m.add (i, i - 1, 1);
			m.add (i, i,     4);
		}
	}
}

template < typename T >
bool test_solve (int n, int iters)
{
	int i, j = 0;
	if (n <= 0)     n     = 320000;
	if (iters <= 0) iters = 1;

	bool ret = true;

	vector < T > b (n);
	vector < T > x1 (n);
	vector < T > x2 (n);
	vector < T > x3 (n);
	vector < T > v (n);
	fprintf (stderr, "n=%d, iters=%d\n", n, iters);

	/* матрицу записываем по строкам! */
	SparseSolver  < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > > M1 (n);
	init_matrix (M1, n);

#ifdef UMFPACK
	UmfPackSolver < T, StoreCSR < T , Allocator > > M2 (n);
	init_matrix (M2, n);
#endif

#ifdef SUPERLU
	SuperLUSolver < T, StoreCSR < T , Allocator > > M3 (n);
	init_matrix (M1, n);
#endif

	for (i = 0; i < n; ++i)
	{
		b[i] = 1;
	}

	Timer t;

	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M1.solve (&x1[0], &b[0]);
	}
	fprintf (stdout, "gmres solve: %lf\n", t.elapsed() );
#ifdef UMFPACK
	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M2.solve (&x2[0], &b[0]);
	}
	fprintf (stdout, "umfpack solve: %lf\n", t.elapsed() );
#endif
#ifdef SUPERLU
	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M3.solve (&x3[0], &b[0]);
	}
	fprintf (stdout, "superlu solve: %lf\n", t.elapsed() );
#endif
	T nr;
#if defined(UMFPACK) && defined(SUPERLU)
	vec_diff (&v[0], &x2[0], &x3[0], (int) x1.size() );
	nr = vec_norm2 (&v[0], (int) x2.size() );
	fprintf (stderr, "%.16le\n", (double) nr);
	ret &= check (nr);
#endif
	vec_diff (&v[0], &x1[0], &x2[0], (int) x1.size() );
	nr = vec_norm2 (&v[0], (int) x2.size() );
	fprintf (stderr, "%.16le\n", (double) nr);
//	ret &= check(nr);

	return ret;
}

template < typename M, typename V >
void mult_loop(V & x1, M & M1, V & b, int iters)
{
	Timer t;
	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M1.mult_vector (&x1[0], &b[0]);
	}
	fprintf (stdout, "M1 mult_vector: %lf\n", t.elapsed() );
}

template < typename T >
bool test_mult (int n, int iters, const Mesh & m, const string & operator_name)
{
	int i, j = 0;
	if (n <= 0)     n     = 320000;
	if (iters <= 0) iters = 10000;

	bool ret = true;

	vector < T > b (n);
	vector < T > x1 (n);
	vector < T > x2 (n);
	vector < T > x3 (n);
	vector < T > v (n);
	fprintf (stderr, "n=%d, iters=%d\n", n, iters);

	for (i = 0; i < n; ++i)
	{
		b[i] = 1;
	}

	if (operator_name.empty()) {
		/* матрицу записываем по строкам! */
		SparseSolver  < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > > M1 (n);
		init_matrix (M1, n);
		mult_loop(x1, M1, b, iters);
	} else {
		if (operator_name == "laplace") {
			Laplace < T > lapl(m);
			mult_loop(x1, lapl.laplace_, b, iters);
		} else if (operator_name == "slaplace") {
			SphereLaplace < T > lapl(m);
			mult_loop(x1, lapl.laplace_, b, iters);
		} else {
			fprintf(stderr, "unknown operator name %s\n", operator_name.c_str());
			usage(0);
		}
	}

	return ret;
}

template < typename T >
bool test_simple_mult (int n, int iters)
{
	int i, j = 0, k;
	double mx = 0;
	if (n <= 0)     n     = 10000;
	if (iters <= 0) iters = 100;

	bool ret = true;

	vector < T > b (n);
	vector < T > x1 (n);
	vector < T > x2 (n);
	vector < T > x3 (n);
	vector < T > v (n);

	fprintf (stderr, "n=%d, iters=%d\n", n, iters);

	/* матрицу записываем по строкам! */
	vector < T > M1 (n*n);

	srand(time(0));
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			//if (i > j) {
			//	M1[i * n + j] = 1;
			//} else if (i < j) {
			//	M1[i * n + j] = 2;
			//}
			//M1[i * n + j] = (double)rand() / (double)RAND_MAX;
			M1[i * n + j] = 1. / (double) (i + j + 1);
		}
	}

	for (i = 0; i < n; ++i)
	{
		b[i] = 1;
	}

	Timer t;

	t.restart();

	for (k = 0; k < iters; ++k)
	{
		mat_mult_vector (&x1[0], &M1[0], &b[0], n);
	}
	fprintf (stdout, "M1 (simple) mult_vector: %lf\n", t.elapsed() );
#if 0 
	mat_mult_vector_stupid(&x2[0], &M1[0], &b[0], n);

	for (i = 0; i < n; ++i) {
		double a = fabs(x2[i] - x1[i]);
		if (mx < a) {
			mx = a;
		}
		fprintf(stderr, "%.16lf %.16lf\n", x1[i], x2[i]);
	}
	fprintf(stderr, "diff = %.16le\n", mx);
#endif

	return ret;

}

int main (int argc, char * argv[])
{
	Mesh mesh;

	string operator_name;
	string mesh_file;

	bool result        = true;
	bool use_double    = false;
	bool mult_dense    = false;
	bool mult_sparse   = false;
	bool invert_sparse = false;

	int dim   = 0;
	int iters = 0;

	for (int i = 0; i < argc; ++i)
	{
		if (!strcmp (argv[i], "--threads") || !strcmp (argv[i], "-t") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			int threads = atoi (argv[i + 1]);

			fprintf(stderr, "threads=%d\n", threads);
			set_num_threads (threads);
		}
		if (!strcmp (argv[i], "--double") || !strcmp (argv[i], "-d") )
		{
			use_double    = true;
		}
		if (!strcmp (argv[i], "--task") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			if (!strcmp (argv[i + 1], "mult_dense") )
			{
				mult_dense    = true;
			}
			if (!strcmp (argv[i + 1], "mult_sparse") )
			{
				mult_sparse   = true;
			}
			if (!strcmp (argv[i + 1], "invert_sparse") )
			{
				invert_sparse = true;
			}
		}
		if (!strcmp (argv[i], "--dim") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			dim = atoi (argv[i + 1]);
		}
		if (!strcmp (argv[i], "--iters") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			iters = atoi (argv[i + 1]);
		}
		if (!strcmp (argv[i], "--operator"))
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			operator_name = argv[i + 1];
		}
		if (!strcmp (argv[i], "--mesh"))
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			FILE * f = fopen(argv[i + 1], "rb");
			if (!f) {
				usage(argv[0]);
			}

			mesh.load(f);

			fclose(f);

			if (mesh.ps.empty() )
			{
				usage (argv[0]);
			}
		}
		if (!strcmp (argv[i], "--help") || !strcmp (argv[i], "-h") )
		{
			usage (argv[0]);
		}
	}

	try
	{
		phelm_init();

		int has_double = check_device_supports_double();
		fprintf (stderr, "has double: %d\n", has_double);

		Timer t;
		burn (1.0);
		t.restart();
		if (use_double && has_double)
		{
			fprintf (stderr, "testing double:\n");

			if (invert_sparse)
			{
				t.restart();
				result &= test_solve < double > (dim, iters);
				fprintf (stderr, "test_solve < double > (): %lf, %d\n", t.elapsed(), (int) result);
			}

			if (mult_sparse)
			{
				t.restart();
				result &= test_mult < double > (dim, iters, mesh, operator_name);
				fprintf (stderr, "test_mult < double > (): %lf, %d\n", t.elapsed(), (int) result);
			}

			if (mult_dense)
			{
				t.restart();
				result &= test_simple_mult < double > (dim, iters);
				fprintf (stderr, "test_simple_mult < double > (): %lf, %d\n", t.elapsed(), (int) result);
			}
		}
		else
		{
			fprintf (stderr, "testing float:\n");

			if (invert_sparse)
			{
				t.restart();
				result &= test_solve < float > (dim, iters);
				fprintf (stderr, "test_solve < float > (): %lf, %d\n", t.elapsed(), (int) result);
			}

			if (mult_sparse)
			{
				t.restart();
				result &= test_mult < float > (dim, iters, mesh, operator_name);
				fprintf (stderr, "test_mult < float > (): %lf, %d\n", t.elapsed(), (int) result);
			}

			if (mult_dense)
			{
				t.restart();
				result &= test_simple_mult < float > (dim, iters);
				fprintf (stderr, "test_simple_mult < float > (): %lf, %d\n", t.elapsed(), (int) result);
			}
		}
		fprintf (stderr, "elapsed: %lf\n", t.elapsed() );

		phelm_shutdown();
	}
	catch (const std::exception & e)
	{
		fprintf (stderr, "exception: %s\n", e.what() );
	}
	return (int) (!result);
}

