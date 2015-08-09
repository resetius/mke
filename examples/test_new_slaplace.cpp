#include "mesh.h"
#include "phelm.h"

using namespace phelm;

static void usage(const char * name)
{
	fprintf(stderr, "usage: %s [-f|--file mesh.txt|-]\n", name);
	exit(1);
}


static double id_cb(const Triangle::NewElem & phi_i,
					const Triangle::NewElem & phi_j,
					const Triangle & trk,
					int z,
					const Mesh & m,
					int point_i, int point_j,
					int, int,
					void *)
{
	printf("(%d, %d) -> \n", point_i, point_j);
	FuncPtr f = phi_i.f * phi_j.f;
	f->print();
	printf(" <- \n\n");
	exit(0);
	return 0.0;
}

void test_laplace(const Mesh & mesh)
{
	typedef phelm::Solver<double> Matrix;
	Matrix idt(mesh.inner_size);

	new_generate_matrix(idt, mesh, id_cb, (void*)0);
}

int test_new_slaplace(int argc, char * argv[])
{
	Mesh mesh;

	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
		{
			usage(argv[0]);
		}
		else if (!strcmp(argv[i], "--file") || !strcmp(argv[i], "-f"))
		{
			if (i == argc - 1)
			{
				usage(argv[0]);
			}

			FILE * f = (strcmp(argv[i + 1], "-") == 0) ? stdin : fopen(argv[i + 1], "rb");

			if (!f)
			{
				usage(argv[0]);
			}
			if (!mesh.load(f))
			{
				usage(argv[0]);
			}

			fclose(f);
		}
	}

	if (mesh.ps.empty())
	{
		usage(argv[0]);
	}

	mesh.info();
	Timer t;
	linal_init();

	try {
		test_laplace(mesh);
	}
	catch (const std::exception & e)
	{
		fprintf(stderr, "exception: %s\n", e.what());
	}

	fprintf(stderr, "elapsed: %lf\n", t.elapsed());

	linal_shutdown();

	return 0;
}
