/* -*- charset: utf-8 -*- */
 
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vector>

#include "phelm.h"
#include "util.h"
#include "laplace.h"

using namespace std;
using namespace phelm;

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [mesh.txt|-]\n", name);
	exit(1);
}

double u(double x, double y)
{
	return sin(M_PI * x) * sin(M_PI * y) + 1.0;
}

double v(double x, double y)
{
	return x * (x - 1) * y * (y - 1) + 1.0;
}

double f(double x, double y)
{
	return -2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) 
		+ v(x, y);
}

double g(double x, double y)
{
	return 2.0 * y * y - 2.0 * y + 2.0 * x * x - 2.0 * x
		+ u(x, y);
}

static elements_t
laplace_integrate_cb( const Polynom & phi_i,
                      const Polynom & phi_j, 
                      const Triangle & trk, 
                      const Mesh & m,
                      int point_i, int point_j,
                      int i, int j,
                      void * user_data)
{
	int rs = (int)m.inner.size();
	elements_t r;
	double a;

	a = laplace(phi_i, phi_j, trk, m.ps);
	// Delta u
	r.push_back(Element(i,      j,      a));
	// Delta v
	r.push_back(Element(i + rs, j + rs, a));

	a = integrate(phi_i * phi_j, trk, m.ps);
	// v
	r.push_back(Element(i,      j + rs, a));
	// u
	r.push_back(Element(i + rs, j,      a));
	return r;
}

struct laplace_right_part_cb_data
{
	double * F;
	double * G;

	double * BU;
	double * BV;
};

static elements_t 
laplace_right_part_cb( const Polynom & phi_i,
                       const Polynom & phi_j,
                       const Triangle & trk, 
                       const Mesh & m,
                       int point_i, int point_j,
                       int i, int j,
                       laplace_right_part_cb_data * d)
{
	int rs = (int)m.inner.size();
	elements_t r;
	const double * F = d->F;
	const double * G = d->G;

	if (m.ps_flags[point_j] == 1) {    
		int j0       = m.p2io[point_j];
		const double * BU = d->BU;
		const double * BV = d->BV;
		double a = laplace(phi_i, phi_j, trk, m.ps);
		double b = integrate(phi_i * phi_j, trk, m.ps);

		r.push_back(Element(i,      j, -(BU[j0] * a + BV[j0] * b)));
		r.push_back(Element(i + rs, j, -(BV[j0] * a + BU[j0] * b)));
	} else {
		double a = integrate(phi_i * phi_j, trk, m.ps);
		// F
		r.push_back(Element(i,      j, F[point_j] * a));
		// G
		r.push_back(Element(i + rs, j, G[point_j] * a));
	}

	return r;
}

/**
 * \Delta u + v = f(x, y)
 * u + \Delta v = g(x, y)
 */
void test_invert(Mesh & m)
{
	int sz = (int)m.ps.size();
	int rs = (int)m.inner.size();
	int os = (int)m.outer.size();

	// right part
	vector < double > F(sz);
	vector < double > G(sz);

	// real answer
	vector < double > RU(sz);
	vector < double > RV(sz);

	// answer
	vector < double > U(sz);
	vector < double > V(sz);

	// boundary
	vector < double > BU(os);
	vector < double > BV(os);

	// right part
	proj(&F[0], m, f);
	proj(&G[0], m, g);

	// real answer
	proj(&RU[0], m, u);
	proj(&RV[0], m, v);

	// boundary
	proj_bnd(&BU[0], m, u);
	proj_bnd(&BV[0], m, v);

	Matrix A(2 * rs);
	vector < double > RP(2 * rs);
	vector < double > Ans(2 * rs);

	generate_matrix(A, m, laplace_integrate_cb, (void*)0);

	laplace_right_part_cb_data data;
	data.F = &F[0];
	data.G = &G[0];
	data.BU = &BU[0];
	data.BV = &BV[0];
	generate_right_part(&RP[0], m, laplace_right_part_cb, &data);

	A.solve(&Ans[0], &RP[0]);
	p2u(&U[0], &Ans[0],  &BU[0], m);
	p2u(&V[0], &Ans[rs], &BV[0], m);

	fprintf(stdout, "answer nev: U = %le\n", dist(&U[0], &RU[0], m));
	fprintf(stdout, "answer nev: V = %le\n", dist(&V[0], &RV[0], m));
}

int main(int argc, char *argv[])
{
	Mesh mesh;
	if (argc > 1) {
		FILE * f = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "rb");
		if (!f) {
			usage(argv[0]);
		}
		if (!mesh.load(f)) {
			usage(argv[0]);
		}
		fclose(f);
	} else {
		usage(argv[0]);
	}

	mesh.info();
	test_invert(mesh);
	return 0;
}

