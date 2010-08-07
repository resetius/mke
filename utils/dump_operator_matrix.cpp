/* -*- charset: utf-8 -*- */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "util.h"
#include "solver.h"
#include "phelm.h"
#include "laplace.h"
#include "slaplace.h"

using namespace phelm;
using namespace linal;
using namespace std;

void usage (const char * n)
{
	fprintf (stderr, "%s [--double|-d] [--mesh] [--operator]\n", n);
	fprintf (stderr, "--threads n - sets the number of threads\n");
	fprintf (stderr, "--double - use double or float ? \n");
	fprintf (stderr, "--mesh - mesh \n");
	fprintf (stderr, "--operator - operator \n");
	exit (-1);
}

template < typename Matrix >
void dump (const string & filename, Matrix & m)
{
	FILE * f = fopen (filename.c_str(), "wb");
	if (f)
	{
		m.dump (f);
		fclose (f);
	}
}

template < typename T >
void dump_operator_matrix (const string & filename, const Mesh & m, const string & operator_name)
{
	/* матрицу записываем по строкам! */

	if (operator_name == "laplace")
	{
		Laplace < T > lapl (m);
		dump (filename, lapl.laplace_);
	}
	else if (operator_name == "slaplace")
	{
		SphereLaplace < T > lapl (m);
		dump (filename, lapl.laplace_);
	}
	else
	{
		fprintf (stderr, "unknown operator name %s\n", operator_name.c_str() );
		usage (0);
	}
}

int main (int argc, char * argv[])
{
	Mesh mesh;

	bool use_double = false;

	string filename;
	string operator_name;
	string mesh_file;

	for (int i = 0; i < argc; ++i)
	{
		if (!strcmp (argv[i], "--double") || !strcmp (argv[i], "-d") )
		{
			use_double    = true;
		}
		if (!strcmp (argv[i], "--operator") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			operator_name = argv[i + 1];
		}
		if (!strcmp (argv[i], "--mesh") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			FILE * f = fopen (argv[i + 1], "rb");
			if (!f)
			{
				usage (argv[0]);
			}

			mesh.load (f);

			fclose (f);

			if (mesh.ps.empty() )
			{
				usage (argv[0]);
			}
		}
		if (!strcmp (argv[i], "-o") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}
			filename = argv[i + 1];
		}
		if (!strcmp (argv[i], "--help") || !strcmp (argv[i], "-h") )
		{
			usage (argv[0]);
		}
	}

	try
	{
		if (use_double)
		{
			dump_operator_matrix < double > (filename, mesh, operator_name);
		}
		else
		{
			dump_operator_matrix < float > (filename, mesh, operator_name);
		}
	}
	catch (const std::exception & e)
	{
		fprintf (stderr, "exception: %s\n", e.what() );
	}
	return 0;
}

