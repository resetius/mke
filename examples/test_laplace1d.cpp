#include <stdio.h>
#include <vector>

#include "linal.h"
#include "solver.h"
#include "polynom.h"
#include "util.h"

using namespace std;
using namespace linal;

static double func(double x)
{
	//return 4 - x * x;
	return x * x;
}

static double lapl_func(double x)
{
	//return -2.0;

	return 2.0;
}

static double fx(double x, double * d)
{
	//return (x - d->xl) * (-x + d->xc);
	// dx, dx
	return *d;
}

struct fx_data
{
	double xl;
	double xc;
	double xr;
	int type;
};

static double fx2(double x, struct fx_data * d)
{
	switch (d->type) {
	case 0:
		return (x - d->xl) * (x - d->xl);
	case 1:
		return (-x + d->xr) * (-x + d->xr);
	case 2:
		return (x - d->xl) * (-x + d->xc);
	case 3:
		return (-x + d->xr) * (x - d->xc);
	default:
		throw "error";
	}
	return 0.0;
}

static void calc_lapl()
{
	double x0  = -2.0;
	double x1  =  2.0;
	int points =  5; // including edges
	int rs     = points - 2; // inner points
	double dx  = (x1 - x0) / (points - 1);

	linal::SparseSolver < double, 
		linal::StoreCSR < double , linal::Allocator > , 
		linal::StoreCSR < double , linal::Allocator > > A(rs);
	vector < double > rp(rs);
	vector < double > ans(rs);
	vector < double > rans(rs);

	// inner points
	for (int i = 0; i < rs; ++i) {
		int p = i + 1;

		double xc = p * dx + x0;
		double xl = (p - 1) * dx + x0;
		double xr = (p + 1) * dx + x0;
		double a;
		double b;

		int p2 = p;
		int j  = p2 - 1;
		double k;

		fx_data data;
		data.xc = xc;
		data.xl = xl;
		data.xr = xr;

		rans[i] = lapl_func(xc);

		k  = 1.0;
		b  = gauss_kronrod15(xl, xc, (fx_t)fx, &k);
		b += gauss_kronrod15(xc, xr, (fx_t)fx, &k);
		rp[i] = -func(xc) * b;

		data.type = 0;
		a  = gauss_kronrod15(xl, xc, (fx_t)fx2, &data);
		A.add(i, j, a);
		data.type = 1;
		a  = gauss_kronrod15(xc, xr, (fx_t)fx2, &data);
		A.add(i, j, a);

		p2 = p - 1;
		j  = p2 - 1;
		k  = -1.0;

		b  = gauss_kronrod15(xl, xc, (fx_t)fx, &k);
		rp[i] += -func(xl) * b;

		data.type = 2;
		a  = gauss_kronrod15(xl, xc, (fx_t)fx2, &data);
		if (j >= 0) {
			A.add(i, j, a);
		} else {
		}

		p2 = p + 1;
		j  = p2 - 1;

		b = gauss_kronrod15(xc, xr, (fx_t)fx, &k);
		rp[i] += -func(xr) * b;

		data.type = 3;
		a  = gauss_kronrod15(xc, xr, (fx_t)fx2, &data);
		if (j < rs) {
			A.add(i, j, a);
		} else {
		}
	}

	A.print(stderr);
	A.solve(&ans[0], &rp[0]);

	fprintf(stderr, "done\n");
}

static void solve_lapl()
{
	double x0  = -2.0;
	double x1  =  2.0;
	int points =  5; // including edges
	int rs     = points - 2; // inner points
	double dx  = (x1 - x0) / (points - 1);

	linal::SparseSolver < double, 
		linal::StoreCSR < double , linal::Allocator > , 
		linal::StoreCSR < double , linal::Allocator > > A(rs);
	vector < double > rp(rs);
	vector < double > ans(rs);
	vector < double > rans(rs);

	// inner points
	for (int i = 0; i < rs; ++i) {
		int p = i + 1;

		double xc = p * dx + x0;
		double xl = (p - 1) * dx + x0;
		double xr = (p + 1) * dx + x0;
		double a;
		double b;

		fx_data data;
		data.xc = xc;
		data.xl = xl;
		data.xr = xr;

		rans[i] = func(xc);
		rp[i] = 0; // - lapl_func(xc);

		// left:   x - xl; /
		// right: -x + xr; \
		// left left:   -x + xc; \
		// right right:  x - xc; /

		int p2 = p;
		int j  = p2 - 1;
		double k;

		k = 1.0;
		a  = gauss_kronrod15(xl, xc, (fx_t)fx, &k);
		a += gauss_kronrod15(xc, xr, (fx_t)fx, &k);
		A.add(i, j, a);

		// b = int (x - xl) (x + xl)
		data.type = 0;
		b  = gauss_kronrod15(xl, xc, (fx_t)fx2, &data);
		rp[i] += -lapl_func(xc) * b;
		// b = int (-x + xr) (-x + xr)
		data.type = 1;
		b  = gauss_kronrod15(xc, xr, (fx_t)fx2, &data);
		rp[i] += -lapl_func(xc) * b;

		p2 = p - 1;
		j  = p2 - 1;
		k  = -1.0;

		a  = gauss_kronrod15(xl, xc, (fx_t)fx, &k);
		if (j >= 0) {
			A.add(i, j, a);
		} else {
			rp[i] -= func(xl) * a;
		}

		data.type = 2;
		b  = gauss_kronrod15(xl, xc, (fx_t)fx2, &data);
		rp[i] += -lapl_func(xl) * b;

		p2 = p + 1;
		j  = p2 - 1;

		a = gauss_kronrod15(xc, xr, (fx_t)fx, &k);
		if (j < rs) {
			A.add(i, j, a);
		} else {
			rp[i] -= func(xr) * a;
		}

		data.type = 3;
		b  = gauss_kronrod15(xc, xr, (fx_t)fx2, &data);
		rp[i] += -lapl_func(xr) * b;
	}

	A.solve(&ans[0], &rp[0]);
//	A.print(stderr);

	vec_diff(&rans[0], &rans[0], &ans[0], rs);
	double nr = vec_norm2(&rans[0], rs);

	fprintf(stderr, "norm = %.16le\n", nr);
	fprintf(stderr, "done\n");
}

int main(int argc, char ** argv)
{
	//solve_lapl();
	calc_lapl();
}
