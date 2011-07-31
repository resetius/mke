#include "deriv.h"
#include "linal.h"
#include "phelm.h"
#include "norm.h"

using namespace std;
using namespace phelm;

static void usage (const char * name)
{
	fprintf (stderr, "usage: %s [-f|--file mesh.txt|-] [-d|--double] [-t|--threads number]\n", name);
	exit (1);
}

static double func(double x, double y)
{
	return x * x + y * y;
}

static double func_x(double x, double y)
{
	return 2 * x + y * y;
}

static double func_y(double x, double y)
{
	return x * x + 2 * y;
}

template < typename T >
void test_deriv(const Mesh & mesh)
{
	int n = mesh.size;

	FlatNorm < T > norm(mesh);
	Deriv < T > deriv(mesh);
	ArrayHost < T > ff(n);
	ArrayHost < T > ff_x(n);
	ArrayHost < T > ff_y(n);

	ArrayHost < T > ff_rx(n);
	ArrayHost < T > ff_ry(n);

	proj(&ff[0], mesh, func);
	proj(&ff_rx[0], mesh, func_x);
	proj(&ff_ry[0], mesh, func_y);

	deriv.calc_x(&ff_x[0], &ff[0]);
	deriv.calc_y(&ff_y[0], &ff[0]);

	fprintf(stderr, "dist 1 = %.16le\n",
		norm.dist(&ff_x[0], &ff_rx[0]));
	fprintf(stderr, "dist 2 = %.16le\n",
		norm.dist(&ff_y[0], &ff_ry[0]));
}

int main(int argc, char * argv[])
{
	Mesh mesh;
	bool use_double = false;

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
		else if (!strcmp (argv[i], "--double") || !strcmp (argv[i], "-d") )
		{
			use_double = true;
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
			test_deriv < double > (mesh);
		}
		else
		{
			fprintf (stderr, "using float\n");
			test_deriv < float > (mesh);
		}
	}
	catch (const std::exception & e)
	{
		fprintf (stderr, "exception: %s\n", e.what() );
	}

	fprintf (stderr, "elapsed: %lf\n", t.elapsed() );

	linal_shutdown();
}
