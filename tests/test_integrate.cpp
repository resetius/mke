/* -*- charset: utf-8 -*- */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>

#include "phelm.h"
#include "util.h"

using namespace std;
using namespace phelm;

static void check (double a, double b)
{
	if (fabs (a - b) > 1e-14)
	{
		fprintf (stderr, "%.16lf != %.16lf\n", a, b);
		exit (-1);
	}
}

#if 1
static void test_simple()
{
	double a;
	fprintf (stderr, "test_simple\n");
	MeshPoint p1 (0, 0), p2 (1, 0), p3 (1, 1);
	Mesh::points_t ps;
	ps.push_back (p1);
	ps.push_back (p2);
	ps.push_back (p3);
	Triangle t1 (0, 1, 2, ps);
	Polynom p = (P2X - 1) * (P2Y - 2);
	fprintf (stderr, "triangle : [0, 0]-[1, 0]-[1, 1]\n");
	fprintf (stderr, "(x - 1) (y - 2) = ");
	p.print();
	//-> 0.291(6)
	a = integrate (p, t1, 0);
	fprintf (stderr, "integral (x - 1) (y - 2) = %.16lf\n", a);
	check (a, 0.2916666666666666);

	fprintf (stderr, "done\n");
}

static void test_simple2()
{
	fprintf (stderr, "int 0.5->1 int x->1 dx dy\n");
	double a = trapezoid_integral (0, 0, 1, 0, 0, 1, 0.5, 1);
	fprintf (stderr, "int 0.5->1 int x->1 dx dy = %.16lf\n", a);
	check (a, 0.125);
}

static void test_simple3()
{
	double ans = 0.0;
	fprintf (stderr, "test_simple\n");
	MeshPoint p1 (0, 0), p2 (1, 0), p3 (1, 1);
	Mesh::points_t ps;
	ps.push_back (p1);
	ps.push_back (p2);
	ps.push_back (p3);
	Triangle t1 (0, 1, 2, ps);
	double r;

	Polynom p = (P2X - 1) * (P2Y - 2);

	fprintf (stderr, "triangle : [0, 0]-[1, 0]-[1, 1]\n");
	fprintf (stderr, "(x - 1) (y - 2) = ");
	p.print();

	r = integrate_cos(p, t1, 0);
	fprintf (stderr, "integral (x - 1) (y - 2) cos x = %.16lf\n", r );
	ans = (4.0 * sin (1.) - 9.0 * cos (1.) ) / 2.0 + 1.0;
	check(r, ans);
	r = integrate_generic_new(p1.pr, p2.pr, p3.pr, [p](double x, double y){
		return p.apply(x, y) * cos(x);
	}, 5);
	fprintf(stderr, "integral (x - 1) (y - 2) cos x = %.16lf\n", r);
	check(r, ans);
	fprintf (stderr, "                    double ans = %.16lf\n", ans);

	r = integrate_sin(p, t1, 0);
	fprintf (stderr, "integral (x - 1) (y - 2) sin x = %.16lf\n", r );
	ans = 5.0 - (9.0 * sin (1.) + 4.0 * cos (1.) ) / 2.0;
	check(r, ans);
	fprintf (stderr, "                    double ans = %.16lf\n", ans);

	fprintf (stderr, "integral (x - 1) (y - 2) / cos x = %.16lf\n", integrate_1_cos (p, t1, 0) );

	fprintf (stderr, "done\n");
}

static void print_poly (vector < Polynom > & p)
{
	for (size_t i = 0; i < p.size(); ++i)
	{
		p[i].print();
	}
}
#if 0
static void test_elem()
{
	Mesh m;
	FILE * f = fopen ("input_0.txt", "rb");
	if (!f) return;
	m.load (f);
	Triangle & t1 = m.tr[0];
	vector < Polynom > p1 = m.elem1 (t1);
	print_poly (p1);
	fprintf (stderr, "\n");

	std::swap (t1.p[1], t1.p[2]);
	vector < Polynom > p2 = m.elem1 (t1);
	print_poly (p2);
	fprintf (stderr, "\n");

	std::swap (t1.p[0], t1.p[2]);
	vector < Polynom > p3 = m.elem1 (t1);
	print_poly (p3);
	fprintf (stderr, "\n");

	fclose (f);
}
#endif
static void test_trapezoid_cos()
{
	double a, b;
	double step = M_PI / 8;

	a = trapezoid_integral_cos (0, 0, -1.0,
	                            0.78539816339744828, 0, 0, 0, 0.78539816339744828);
	fprintf (stderr, "a=%.16lf\n", a);

	a = trapezoid_integral_cos (0, 0, -1.0,
	                            0.78539816339744828 + step, 0, step, 0, 0.78539816339744828);
	fprintf (stderr, "a=%.16lf\n", a);

	b = trapezoid_integral_1_cos (1, 0, -1.0,
	                              0.78539816339744828, 0, 0, 0, M_PI / 2);
	fprintf (stderr, "b=%.16lf\n", b);
}

static void test_sphere_area()
{
	fprintf (stderr, "test_sphere_area\n");

	double pp[] = { 1.0 };
	Polynom p (0, 0, pp, 1); // f(x, y) = 1

	double S1, S2, S3, S4;

	// треугольник, состоящий из 4х треугольников:
	//   S2
	//   S4
	// S1  S3
	{
		MeshPoint p1 (0, 0), p2 (0, M_PI / 8.0), p3 (M_PI / 8.0, 0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S1 = integrate_cos (p, t1, 0);

		fprintf (stderr, "S1 = %.16lf\n", S1);
		check (S1, 0.0761204674887130);
	}

	{
		MeshPoint p1 (0, M_PI / 8.0), p2 (0, M_PI / 4.0), p3 (M_PI / 8, M_PI / 8.0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S3 = integrate_cos (p, t1, 0);

		fprintf (stderr, "S3 = %.16lf\n", S3);
		check (S3, 0.0761204674887130);
	}

	{
		MeshPoint p1 (0, M_PI / 8.0), p2 (M_PI / 8, 0), p3 (M_PI / 8, M_PI / 8.0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S4 = integrate_cos (p, t1, 0);

		fprintf (stderr, "S4 = %.16lf\n", S4);
		check (S4, 0.0741589649823731);
	}

	{
		MeshPoint p1 (M_PI / 8.0, 0), p2 (M_PI / 4, 0), p3 (M_PI / 8, M_PI / 8.0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S2 = integrate_cos (p, t1, 0);
		//S2 = 0.460075592;

		fprintf (stderr, "S2 = %.16lf\n", S2);
		check (S2, 0.0664933188536525);
	}

	fprintf (stderr, "S1 + S2 + S3 + S4 = %.16lf\n", S1 + S2 + S3 + S4);
	fprintf (stderr, "S1 +      S3      = %.16lf\n", S1 + S3);
	fprintf (stderr, "                  = %.16lf\n", M_PI / 2);

	fprintf (stderr, "... = %.16lf\n", M_PI / 4. * (sin (M_PI / 2.) - sin (M_PI / 4.) ) );


	// swap x<->y
	{
		MeshPoint p1 (0, 0), p2 (M_PI / 8.0, 0), p3 (0, M_PI / 8.0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S1 = integrate_cos (p, t1, 0);

		fprintf (stderr, "S1 = %.16lf\n", S1);
		check (S1, 0.0761204674887130);
	}

	{
		MeshPoint p1 (M_PI / 8.0, 0), p2 (M_PI / 4.0, 0), p3 (M_PI / 8, M_PI / 8.0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S3 = integrate_cos (p, t1, 0);

		fprintf (stderr, "S3 = %.16lf\n", S3);
		check (S3, 0.0664933188536525);
	}

	{
		MeshPoint p1 (M_PI / 8.0, 0), p2 (0, M_PI / 8), p3 (M_PI / 8, M_PI / 8.0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S4 = integrate_cos (p, t1, 0);

		fprintf (stderr, "S4 = %.16lf\n", S4);
		check (S4, 0.0741589649823731);
	}

	{
		MeshPoint p1 (0, M_PI / 8.0), p2 (0, M_PI / 4), p3 (M_PI / 8, M_PI / 8.0);
		Mesh::points_t ps;
		ps.push_back (p1);
		ps.push_back (p2);
		ps.push_back (p3);
		Triangle t1 (0, 1, 2, ps);
		S2 = integrate_cos (p, t1, 0);
		//S2 = 0.460075592;

		fprintf (stderr, "S2 = %.16lf\n", S2);
		check (S1, 0.0761204674887130);
	}

	fprintf (stderr, "done\n");
}

#endif

extern "C" int test_integrate (int argc, char * argv[])
{
	test_simple();
	test_simple2();
	//test_elem();
	test_simple3();
	test_sphere_area();
	//test_trapezoid_cos();
	return 0;
}

