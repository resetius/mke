#include <stdio.h>
#include <vector>

#include "linal.h"
#include "solver.h"
#include "polynom.h"
#include "util.h"

using namespace std;

static double func(double x)
{
	return 4 - x * x;
}

static double lapl_func(double x)
{
	return -2.0;
}

static void calc_lapl()
{
}

struct fx_data
{
	double xl;
	double xc;
};

static double fx(double x, struct fx_data * d)
{
	//return (x - d->xl) * (-x + d->xc);
	// dx, dx
	return -1.0;
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

	// inner points
	for (int i = 0; i < rs; ++i) {
		int p = i + 1;

		double xc = p * dx + x0;
		double xl = (p - 1) * dx + x0;
		double xr = (p + 1) * dx + x0;
		double a;

		rp[i] = - lapl_func(xc);

		// left:   x - xl; /
		// right: -x + xr; \
		// left left:   -x + xc; \
		// right right:  x - xc; /

		int p2 = p - 1;
		int j  = p2 - 1;

		struct fx_data data;
		data.xc = xc;
		data.xl = xl;

		a = gauss_kronrod15(xl, xc, (fx_t)fx, &data);

		if (j >= 0) {
			A.add(i, j, a);
		} else {
		//	rp[i] += func(xc) * a;
		}

		p2 = p + 1;
		j  = p2 - 1;

		data.xc = xr;
		data.xl = xc;
		a = gauss_kronrod15(xc, xr, (fx_t)fx, &data);
		if (j < rs) {
			A.add(i, j, a);
		} else {
		//	rp[i] += func(xc) * a;
		}
	}

	A.solve(&ans[0], &rp[0]);
	A.print(stderr);

	fprintf(stderr, "done\n");
}

int main(int argc, char ** argv)
{
	solve_lapl();
}
