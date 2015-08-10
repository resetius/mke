#include "mesh.h"
#include "phelm.h"
#include "norm.h"

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
	FuncPtr f = phi_i.f * phi_j.f;
#if 0
	printf("(%d, %d) -> \n", point_i, point_j);
	f->print();
	printf(" <- \n\n");
	/*
	printf(" conv h:\n");
	phi_i.h.hx->print(); printf("{:\n");
	phi_i.h.hy->print(); printf("{:\n");
	phi_i.h.hz->print(); printf("{:\n");
	*/

	f->bind_args({ "x1", "y1", "z1" });
	FuncPtr fh = f->apply({ phi_i.h.hx, phi_i.h.hy, phi_i.h.hz });
	fh->print();
	printf(" <<- \n\n");
	fh->bind_args({ "x", "y", "z" });
	FuncPtr gfh = fh->apply({ phi_i.g.gx, phi_i.g.gy, phi_i.g.gz });
	gfh->print();
	printf(" <<- \n\n");

	FuncPtr Dgfh = gfh->diff("u");
	Dgfh->print();
	printf(" <<<- \n\n");
#endif

	f->bind_args({ "x1", "y1" });
	return integrate_generic_new(trk, [f](double x, double y) {
		return f->apply({ x, y })->value();
	});
}

static double integrate_cb(const Triangle::NewElem & phi_i,
					const Triangle::NewElem & phi_j,
					const Triangle & trk,
					int z,
					const Mesh & m,
					int point_i, int point_j,
					int, int,
					void *)
{
	// FuncPtr f1 = 
	// x1 <-(u, v)
	// y1 <-(u, v)
	// h(g)
	
	FuncPtr hgx = phi_i.h.hx->apply({ phi_i.g.gx, phi_i.g.gy, phi_i.g.gz });
	FuncPtr hgy = phi_i.h.hy->apply({ phi_i.g.gx, phi_i.g.gy, phi_i.g.gz });
	FuncPtr hgz = phi_i.h.hz->apply({ phi_i.g.gx, phi_i.g.gy, phi_i.g.gz });

	FuncPtr f1 =
		(phi_i.f->diff("x1") * hgx->diff("u") +
		phi_i.f->diff("y1") * hgy->diff("u")) *
		(phi_j.f->diff("x1") * hgx->diff("u") +
		phi_j.f->diff("y1") * hgy->diff("u"));

	// u <- (x1, y1, z1)
	// v <- (x1, y1, z1)

	FuncPtr g1u = phi_i.g1.g1u->apply({ phi_i.h1.h1x, phi_i.h1.h1y, phi_i.h1.h1z });
	FuncPtr g1v = phi_i.g1.g1v->apply({ phi_i.h1.h1x, phi_i.h1.h1y, phi_i.h1.h1z });

	f1->bind_args({ "x1", "y1", "z1", "u", "v" });
	g1u->bind_args({ "x1", "y1", "z1" });
	g1v->bind_args({ "x1", "y1", "z1" });

/*	printf("->\n");
	g1u->print();
	printf("\n");
	g1v->print();
	printf("\n");*/

	
	double r1 = - integrate_generic_new(trk, [f1,g1u,g1v](double x1, double y1) {
		double u = g1u->apply({ x1, y1, 0.0 })->value();
		double v = g1v->apply({ x1, y1, 0.0 })->value();
/*		if (isnan(u) || isnan(v)) {
			printf("%d, %d, x1, y1 -> %le, %le\n", 
				   (int)isnan(u), (int)isnan(v),
				   x1, y1);
		}*/
		return f1->apply({ x1, y1, 0.0, u, v })->value();
	});

/*	printf("%le\n", r);
	
	if (isnan(r)) {
		exit(0);
	}*/


	FuncPtr f2 =
		(phi_i.f->diff("x1") * hgx->diff("v") +
		phi_i.f->diff("y1") * hgy->diff("v")) *
		(phi_j.f->diff("x1") * hgx->diff("v") +
		phi_j.f->diff("y1") * hgy->diff("v"));

	double r2 = -integrate_generic_new(trk, [f1, g1u, g1v](double x1, double y1) {
		double u = g1u->apply({ x1, y1, 0.0 })->value();
		double v = g1v->apply({ x1, y1, 0.0 })->value();
		return f1->apply({ x1, y1, 0.0, u, v })->value() / cos(u) / cos(u);
	});

	return r1+r2;
}

static double func_lapl(double x, double y)
{
	//return -6.0 * sin(y) * sin(2.0 * x);

	double sx = sin(x);
	double cx = cos(x);
	double sy = sin(y);
	double cy = cos(y);

	double ans = 0.0;

	ans += (sx * sx * (cx * cy - 1) * cy);
	ans += -(cx * cy - 1) * cy;
	ans += -2 * sx * (cx * cy - sx) * (-sx * cy - cx);
	ans += -2 * (cx * cy - sx) * cy;
	if (fabs(ans) > 1e-15)
	{
		ans /= cx;
	}

	ans += sx * sx * cy * cy;
	ans += -cx * (cx * cy - 1) * cy;
	ans += sy * sy;
	ans += 2 * ipow(-sx * cy - cx, 2);
	ans += 2 * (cx * cy - sx) * (-cx * cy + sx);
	ans += 2 * sy * sy;
	ans += -8 * (sx - 3.0) * sx + 4.0 * cx * cx;
	return ans;
}

static double func(double x, double y)
{
	//return sin(y) * sin(2.0 * x);

	double sx = sin(x);
	double cx = cos(x);
	//double sy = sin(y);
	double cy = cos(y);

	return 0.5 * ipow(cx * cy - 1.0, 2) +
		ipow(cx * cy - sx, 2) +
		2.0 * ipow(sx - 3, 2);
}

void test_laplace(const Mesh & mesh)
{
	typedef phelm::Solver<double> Matrix;
	Matrix idt(mesh.inner_size);
	Matrix laplace(mesh.inner_size);
	Matrix bnd3(mesh.inner_size);

	new_generate_matrix(idt, mesh, id_cb, (void*)0);
	new_generate_matrix(laplace, mesh, integrate_cb, (void*)0);
	new_generate_boundary_matrix(bnd3, mesh, integrate_cb, (void*)0);

//	idt.print(stdout);
//	laplace.print(stdout);

	std::vector<double> U(mesh.size);
	std::vector<double> LU(mesh.size);
	std::vector<double> LU1(mesh.size);
	std::vector<double> tLU1(mesh.inner_size);
	std::vector<double> B1(std::max(mesh.outer_size, 1));

	proj(&U[0], mesh, func);
	proj(&LU[0], mesh, func_lapl);
	proj_bnd(&B1[0], mesh, func_lapl);

	std::vector<double> in(mesh.inner_size);
	std::vector<double> out(mesh.inner_size);
	std::vector<double> tmp(mesh.outer_size);
	u2p(&in[0], &U[0], mesh);
	proj_bnd(&tmp[0], &U[0], mesh);
	laplace.mult_vector(&out[0], &in[0]);
	bnd3.mult_vector(&in[0], &tmp[0]);
	vec_sum(&out[0], &out[0], &in[0], (int)in.size());
	idt.solve(&tLU1[0], &out[0]);
	p2u(&LU1[0], &tLU1[0], &B1[0], mesh);

	SphereNorm<double> nr(mesh);

	fprintf(stdout, "L2: laplace err=%.2le\n", (double)nr.dist(&LU[0], &LU1[0]));

//	laplace.print(stdout);
//	idt.print(stdout);
//	bnd3.print(stdout);
//	printf("\n\n");
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

	set_num_threads(1);

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
