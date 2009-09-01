/* -*- charset: utf-8 -*- */
 
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "barvortex.h"
#include "util.h"

using namespace std;
using namespace phelm;

void usage(const char * name)
{
	fprintf(stderr, "usage: %s [-f|--file mesh.txt|-] [-t|--threads number] [--task task] [--verbose|-v number]\n", name);
	fprintf(stderr, "tasks:\n"
			"jacobian\n"
			"jacobian_T\n"
			"test\n"
			"L\n"
			"L2\n"
			"LT\n"
			"dymnikov_196\n"
			"kornev1\n"
			"laplace_LT\n");

	exit(1);
}

double f1 (double x, double y)
{
	return sin (x) * sin (y) * cos(x);
}

double f2 (double x, double y)
{
	return 1000.0 * cos (x) * cos (y) * sin(x);
}

double f3 (double x, double y)
{
	return -(-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) ) / cos (x);
}

double an (double x, double y)
{
	return (-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) ) / cos (x);
	//return (-sin(x)*cos(y)*sin(x)*cos(y)) / cos(x);
	//return (-cos(x)*sin(y)*cos(x)*sin(y))/cos(x);

	//return cos(x) * sin(y) / cos(x);
	//return sin(x) * cos(y) / cos(x);
}

double test_jacobian_f1 (double x, double y)
{
	return sin (x) * sin (y);
}

double test_jacobian_f2 (double x, double y)
{
	return cos (x) * cos (y);
}

double test_jacobian_an (double x, double y)
{
	return (-sin (x) *cos (y) *sin (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) ) / cos (x);
}

void test_jacobian (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int rs = (int)m.inner.size();
	int os = (int)m.outer.size();

	SphereJacobian j (m);
	vector < double > F1 (sz);
	vector < double > F2 (sz);
	vector < double > ans1 (sz);
	vector < double > rans1 (sz);
	vector < double > bnd (os);

	proj (&F1[0], m, test_jacobian_f1);
	proj (&F2[0], m, test_jacobian_f2);
	proj (&rans1[0], m, test_jacobian_an);
	proj_bnd (&bnd[0], m, test_jacobian_an);

	j.calc1 (&ans1[0], &F1[0], &F2[0], &bnd[0]);

	fprintf (stdout, "jacobian err=%.2le\n",
	         dist (&ans1[0], &rans1[0], m, sphere_scalar_cb, (void*)0) );

	//vector < double > p1(m.inner.size());
	//u2p(&p1[0], &rans1[0], m);
	//vector_print(&p1[0], p1.size());
	//u2p(&p1[0], &ans1[0], m);
	//vector_print(&p1[0], p1.size());

	//vector_print(&rans1[0], rans1.size());
	//vector_print(&ans1[0], ans1.size());
}

void rand_init(double * h, int n)
{
	for (int i = 0; i < n; ++i) {
		h[i] = (double)rand() / (double)RAND_MAX;
	}
}

void test_jacobian_T (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int rs = (int)m.inner.size();
	int os = (int)m.outer.size();
	double nr1, nr2;

	SphereJacobian j (m);
	vector < double > u (sz);
	vector < double > w (sz);
	vector < double > v (sz);
	vector < double > tmp(rs);
	vector < double > ans1 (sz);
	vector < double > ans2 (sz);

	srand((unsigned int)time(0));
	//proj (&u[0], m, f1);
	//proj (&w[0], m, f2);
	//proj (&v[0], m, f3);
	rand_init(&u[0], (int)u.size()); u2p(&tmp[0], &u[0], m); p2u(&u[0], &tmp[0], 0, m); 
	rand_init(&v[0], (int)v.size()); u2p(&tmp[0], &v[0], m); p2u(&v[0], &tmp[0], 0, m); 
	rand_init(&w[0], (int)w.size()); u2p(&tmp[0], &w[0], m); p2u(&w[0], &tmp[0], 0, m); 

	j.calc1 (&ans1[0], &v[0], &w[0], 0);
	j.calc1t(&ans2[0], &u[0], &w[0], 0);

	nr1 = scalar(&ans1[0], &u[0], m, sphere_scalar_cb, (void*)0);
	nr2 = scalar(&ans2[0], &v[0], m, sphere_scalar_cb, (void*)0);

	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));

	j.calc1 (&ans1[0], &w[0], &v[0], 0);
	j.calc1t(&ans2[0], &w[0], &u[0], 0);

	nr1 = scalar(&ans1[0], &u[0], m, sphere_scalar_cb, (void*)0);
	nr2 = scalar(&ans2[0], &v[0], m, sphere_scalar_cb, (void*)0);

	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));
}

static double x (double u, double v)
{
	return cos (u) * cos (v);
}

static double y (double u, double v)
{
	return cos (u) * sin (v);
}

static double z (double u, double v)
{
	return sin (u);
}

double u0 (double x, double y)
{
	//return 0.0;
	return sin (x) * sin (y);
}

double z0 (double x, double y)
{
	return sin(x) * cos(x) / 25.0;
}

double rp (double x, double y, double t, double mu, double sigma)
{
	// double b = -6.0 * exp(t) * sin(y) * sin(2.0 * x);
	// return a - mu * b + sigma * a;
	// return 0.0;

	// lapl from -3.5 * sigma * ipow (sin (x), 3);
	return sigma*(2.*ipow(cos(x),2)-1.)*sin(x);
}

double zero_coriolis (double phi, double lambda)
{
	return 0.0;
}

double coriolis (double phi, double lambda)
{
	double omg = 0.0000727000000000;
	double l = omg * 2.0 * sin (phi);
	double h = cos (2.0 * lambda) * ipow (sin (2.0 * phi), 2);
	return l + h;
}

double dymnikov_196_coriolis(double phi, double lambda)
{
	double R    = 6371.3;
	double beta = 2e-11;
	double l0   = 9.3e-5;
	return l0 + beta * R * (M_PI / 10 + phi);
}

double dymnikov_196_rp(double phi, double lambda, double t, double mu, double sigma)
{
	double R    = 6371.3;
	double y    = R * (M_PI / 10 + phi);
	double k    = 10.0;
	double tau0 = 1.1;
	double tau  = 1.0; //?
	double rho  = 1000;
	double L    = 4000;
	double H    = 500;
	double ans = -k * 2.0 * M_PI * tau0 / rho / L / H * sin(2 * tau * y / L);
	return ans;
}

double an1 (double x, double y, double t)
{
	return x*sin(y+t)*ipow(cos(x),4);
}

double omg_an1 (double x, double y, double t)
{
	return -sin(y+t)*ipow(cos(x),2)*(9*sin(x)*cos(x)+20*x*ipow(cos(x),2)-15*x);
}

double rp1(double x, double y, double t, double mu, double sigma)
{
	return 390*mu*sin(y+t)*x*ipow(cos(x),2)
		+147*mu*sin(y+t)*sin(x)*cos(x)-
		400*mu*sin(y+t)*x*
		ipow(cos(x),4)-360*mu*sin(y+t)*
		sin(x)*ipow(cos(x),3)-20*sigma*sin(y+t)*
		ipow(cos(x),4)*x+15*sigma*sin(y+t)*
		ipow(cos(x),2)*x-9*sigma*sin(y+t)*
		ipow(cos(x),3)*sin(x)+30*cos(y+t)*
		ipow(cos(x),4)*sin(y+t)*x*x*sin(x)+9*cos(y+t)*
		ipow(cos(x),6)*sin(y+t)*sin(x)-45*mu*sin(y+t)*x-
		9*cos(y+t)*ipow(cos(x),5)*sin(y+t)*x-20*x*cos(y+t)*
		ipow(cos(x),4)-9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2);
}

double rp2(double x, double y, double t, double mu, double sigma)
{
	return -9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2)-20*x*cos(y+t)*
		ipow(cos(x),4)-9*sigma*sin(y+t)*
		ipow(cos(x),3)*sin(x)+15*sigma*sin(y+t)*
		ipow(cos(x),2)*x-20*sigma*sin(y+t)*
		ipow(cos(x),4)*x-360*mu*sin(y+t)*sin(x)*
		ipow(cos(x),3)+390*mu*sin(y+t)*x*
		ipow(cos(x),2)-400*mu*sin(y+t)*x*
		ipow(cos(x),4)+147*mu*sin(y+t)*sin(x)*cos(x)-45*mu*sin(y+t)*x;
}

void test_barvortex_L (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.05;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 8e-5;   //8e-5;
	double sigma = 1.6e-2; //1.6e-2;

	//BarVortex bv (m, rp1, zero_coriolis, tau, sigma, mu, 1.0, 1.0);
	BarVortex bv (m, rp, coriolis, tau, sigma, mu, 1.0, 1.0);

	vector < double > u (sz);
	vector < double > z (sz);
	vector < double > bnd (std::max (os, 1));
	vector < double > Ans(sz);
	vector < double > Ans1(sz);

	//proj (&u[0], m, an1, 0);
	proj (&u[0], m, u0);

	proj (&z[0], m, z0);

	//if (!bnd.empty()) proj_bnd(&bnd[0], m, f1);

	setbuf (stdout, 0);

	bv.L_step(&Ans[0], &u[0], &z[0]);
	bv.L_1_step(&Ans1[0], &Ans[0], &z[0]);

	fprintf (stderr, " ||L(L^1)|| = %le \n",
		dist (&u[0], &Ans1[0], m, sphere_scalar_cb, (void*)0));
}

void test_barvortex_LT (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.05;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 8e-5;   //8e-5;
	double sigma = 1.6e-2; //1.6e-2;

	//BarVortex bv (m, rp1, zero_coriolis, tau, sigma, mu, 1.0, 1.0);
	BarVortex bv (m, rp, coriolis, tau, sigma, mu, 1.0, 1.0);

	vector < double > u  (sz);
	vector < double > v  (sz);
	vector < double > lv (sz);
	vector < double > ltu(sz);

	vector < double > z (sz);
	vector < double > bnd (std::max (os, 1));

	proj (&u[0], m, f1);
	proj (&v[0], m, f2);

	proj (&z[0], m, z0);

	setbuf (stdout, 0);

	bv.L_step (&lv[0],  &v[0], &z[0]);
	bv.LT_step(&ltu[0], &u[0], &z[0]);

	double nr1 = scalar(&lv[0],  &u[0], m, sphere_scalar_cb, (void*)0);
	double nr2 = scalar(&ltu[0], &v[0], m, sphere_scalar_cb, (void*)0);

	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));
}

double rnd1(double x, double y)
{
	if (fabs(x) < 1e-8) {
		return 0;
	} else {
		return (double)rand() / (double)RAND_MAX;
	}
}

void test_laplace_LT (const Mesh & m)
{
	double nr1, nr2;
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	int i = 0;

	//srand(time(0));
	srand(0);
	SphereLaplace < double > l (m);
	SphereJacobian j(m);
	SphereNorm  < double >   s(m);

	vector < double > u  (sz);
	vector < double > v  (sz);
	vector < double > w  (sz);
	vector < double > w1 (sz);
	vector < double > lv (sz);
	vector < double > ltu(sz);

	vector < double > z (sz);
	vector < double > bnd (std::max (os, 1));

	proj (&u[0], m, rnd1);
	proj (&v[0], m, rnd1);
	proj (&w[0], m, f1);

	setbuf (stdout, 0);

	l.calc1 (&lv[0],  &v[0], &bnd[0]);
	l.calc1 (&ltu[0], &u[0], &bnd[0]);

	nr1 = scalar(&lv[0],  &u[0], m, sphere_scalar_cb, (void*)0);
	nr2 = scalar(&ltu[0], &v[0], m, sphere_scalar_cb, (void*)0);

	fprintf (stderr, "Laplace:\n");
	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));

//	lv = v;
//	ltu = u;
	vec_mult(&lv[0],  &v[0], &w[0], sz);
	vec_mult(&ltu[0], &u[0], &w[0], sz);

//	j.calc1 (&lv[0],  &v[0], &w[0], &bnd[0]);
//	j.calc1t(&ltu[0], &u[0], &w[0], &bnd[0]);
//	nr1 = scalar(&u[0],  &lv[0], m, sphere_scalar_cb, (void*)0);	
//	nr2 = scalar(&ltu[0], &v[0], m, sphere_scalar_cb, (void*)0);
	nr1 = s.scalar(&u[0],  &lv[0]);
	nr2 = s.scalar(&ltu[0], &v[0]);

//	nr1 = scalar(&u[0],  &v[0], m, sphere_scalar_cb, (void*)0);
//	nr2 = scalar(&v[0], &u[0], m, sphere_scalar_cb, (void*)0);

	fprintf (stderr, "Jacobian:\n");
	fprintf (stderr, "(Lv, u)  = %le \n", nr1);
	fprintf (stderr, "(v, LTu) = %le \n", nr2);

	fprintf (stderr, " |(Lv, u) - (v, LTu)| = %le \n", fabs(nr1 - nr2));
}

void test_barvortex (const Mesh & m, int verbose, double T)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.001;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 8e-5;   //8e-5;
	double sigma = 1.6e-2; //1.6e-2;

	BarVortex bv (m, rp1, zero_coriolis, tau, sigma, mu, 1.0, 1.0);
	//BarVortex bv (m, rp, coriolis, tau, sigma, mu, 1.0, 1.0);

	vector < double > u (sz);
	vector < double > bnd (std::max (os, 1));
	vector < double > Ans(sz);

	proj (&u[0], m, an1, 0);
	//proj (&u[0], m, u0);

	//if (!bnd.empty()) proj_bnd(&bnd[0], m, f1);

	setbuf (stdout, 0);

	while (t < T)
	{
		Timer tm;
		bv.calc (&u[0], &u[0], &bnd[0], t);

		if (i % 1 == 0) {
			fprintf (stderr, " === NORM = %le, STEP %lf of %lf: %lf\n",
					 bv.norm (&u[0]), t, T, tm.elapsed());
			if (verbose) {

				// 3d print
				print_function (stdout, &u[0], m, x, y, z);
				// flat print
				// print_function (stdout, &u[0], m, 0, 0, 0);
			}
		}

		i += 1;
		t += tau;
#if 1
		{
			proj(&Ans[0], m, an1, t);
			fprintf(stderr, "time %lf/ norm %le\n", t, 
				bv.dist(&u[0], &Ans[0]));
//			print_function (stdout, &Ans[0], m, x, y, z);
		}
#endif
	}

	fprintf(stdout, "bv: norm %le\n",  
			bv.dist(&u[0], &Ans[0]));
}

void test_barvortex_L2 (const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.05;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 8e-5;   //8e-5;
	double sigma = 1.6e-2; //1.6e-2;

	BarVortex bv (m, rp1, zero_coriolis, tau, sigma, mu, 1.0, 1.0);
	//BarVortex bv (m, rp, coriolis, tau, sigma, mu, 1.0, 1.0);

	vector < double > u (sz);
	vector < double > u1(sz);
	vector < double > z(sz);
	vector < double > bnd (std::max (os, 1));
	vector < double > Ans(sz);

	proj (&u[0], m, an1, 0);
	//proj (&u[0], m, u0);

	//if (!bnd.empty()) proj_bnd(&bnd[0], m, f1);

	setbuf (stdout, 0);

	while (t < T)
	{
		Timer tm;
		bv.calc_L (&u1[0], &u[0], &z[0], &bnd[0], t);
#if 0
		if (i % 1 == 0) {
			fprintf (stderr, " === NORM = %le, STEP %lf of %lf: %lf\n",
			         norm (&u[0], m, sphere_scalar_cb), t, T, tm.elapsed());
			// 3d print
			print_function (stdout, &u[0], m, x, y, z);
			// flat print
			// print_function (stdout, &u[0], m, 0, 0, 0);
		}
#endif

		i += 1;
		t += tau;
#if 1
		{
			proj(&Ans[0], m, an1, t);
			fprintf(stderr, "time %lf/ norm %le\n", t, 
				dist(&u1[0], &Ans[0], m, sphere_scalar_cb, (void*)0));
//			print_function (stdout, &Ans[0], m, x, y, z);
		}
//		Sleep(500);
#endif
		u1.swap(u);
	}
}

/**
 * Тест из книги Дымникова.
 * Устойчивость и предсказуемость крупномасштабных атмосферных процессов, Москва, 
 * ИВМ РАН, 2007, 283с
 *
 * Все данные со страницы 196.
 */
void test_dymnikov_196(const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.0001;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double mu    = 1250;
	double sigma = 5e-8; //1.6e-2;

	BarVortex bv (m, dymnikov_196_rp, dymnikov_196_coriolis, tau, sigma, mu, 1.0, 1.0);

	vector < double > u (sz);
	vector < double > bnd (std::max (os, 1));

	//proj (&u[0], m, an1, 0);
	//proj (&u[0], m, u0);

	//if (!bnd.empty()) proj_bnd(&bnd[0], m, f1);

	setbuf (stdout, 0);

	while (t < T)
	{
#if 1
		Timer tm;
		bv.calc (&u[0], &u[0], &bnd[0], t);
		if (i % 1 == 0) {
			fprintf (stderr, " === NORM = %le, STEP %lf of %lf: %lf\n",
			         bv.norm (&u[0]), t, T, tm.elapsed());
			// 3d print
			//print_function (stdout, &u[0], m, x, y, z);
			// flat print
			// print_function (stdout, &u[0], m, 0, 0, 0);
		}
#endif

		i += 1;
		t += tau;
	}
}

double kornev1_rp_(double phi, double lambda)
{
	double omg   = 2.0 * M_PI/24./60./60.;
	double T0    = 1./omg;
	double H     = 5000;
	double sigma = 1./20./2./M_PI;
	double R     = 6.371e+6;
	double f     = - sigma * 180/1.15 * (6*(2*cos(phi)*cos(phi)-1)*sin(phi));

	return f*T0*T0/R/R;
}

double kornev1_rp(double phi, double lambda, double t, double mu, double sigma)
{
	return kornev1_rp_(phi, lambda);
}

double kornev1_coriolis(double phi, double lambda)
{
	double omg  = 2.0 * M_PI/24./60./60.;
	double H    = 1;
	return 2 * omg * sin(phi) + 0.1 * H * cos(2*lambda)*ipow(sin(2*phi),2);
}

double kornev1_u0(double phi, double lambda)
{
	return - 180/1.15 * ipow(sin(phi),3);
}

void test_kornev1(const Mesh & m)
{
	int sz = (int)m.ps.size();
	int os = (int)m.outer.size();

	double tau = 0.00001;
	double t = 0;
	//double T = 0.1;
	double days = 30;
	double T = days * 2.0 * M_PI;
	double month = 30.0 * 2.0 * M_PI;
	int i = 0;

	double R     = 6.371e+6;
	double H     = 5000;
	double sigma = 1./20./2./M_PI;
	double mu    = sigma / 100.;
	double omg   = 2.0 * M_PI/24./60./60.;
	double T0    = 1./omg;
	double k1    = 1.0;
	double k2    = T0/R/R;

//	sigma = sigma * T0;
//	mu    = 10e-2*sigma;

	BarVortex bv (m, kornev1_rp, kornev1_coriolis, tau, sigma, mu, k1, k2);

	vector < double > u (sz);
	vector < double > bnd (std::max (os, 1));

	proj (&u[0], m, kornev1_u0);

	{
		vector < double > tmp(sz);
		proj(&tmp[0], m, kornev1_rp_);

		print_function("kornev1_rp.txt", &tmp[0], m, x, y, z);
		print_function("kornev1_u0.txt", &u[0], m, x, y, z);
	}

	setbuf (stdout, 0);

	while (t < T)
	{
#if 1
		Timer tm;
		bv.calc (&u[0], &u[0], &bnd[0], t);
//		print_function("kornev1_u1.txt", &u[0], m, x, y, z);
//		exit(1);
		if (i % 1 == 0) {
			fprintf (stderr, " === NORM = %le, STEP %lf of %lf: %lf\n",
			         bv.norm (&u[0]), t, T, tm.elapsed());
			// 3d print
//			print_function (stdout, &u[0], m, x, y, z);
			// flat print
			// print_function (stdout, &u[0], m, 0, 0, 0);
		}
#endif

		i += 1;
		t += tau;
	}
}

int main (int argc, char *argv[])
{
	Mesh mesh;
	string task;
	int verbose = 0;
	double time = 1.0;

	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
		{
			usage(argv[0]);
		} else if (!strcmp(argv[i], "--file") || !strcmp(argv[i], "-f")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			FILE * f = (strcmp(argv[i + 1], "-") == 0) ? stdin : fopen(argv[i + 1], "rb");

			if (!f) {
				usage(argv[0]);
			}
			if (!mesh.load(f)) {
				usage(argv[0]);
			}

			fclose(f);
		} else if (!strcmp(argv[i], "--threads") || !strcmp(argv[i], "-t")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			int threads = atoi(argv[i + 1]);
			set_num_threads(threads);
		} else if (!strcmp(argv[i], "--time") || !strcmp(argv[i], "-T")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			time = atof(argv[i + 1]);
		} else if (!strcmp(argv[i], "--task")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}
			task = argv[i + 1];
		} else if (!strcmp(argv[i], "--verbose") || !strcmp(argv[i], "-v")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}
			verbose = atoi(argv[i + 1]);
		}
	}

	if (mesh.ps.empty()) {
		usage(argv[0]);
	}

#if defined(WIN32)
	set_fpe_except();
#endif

	mesh.info();

	if (task == "jacobian") {
		test_jacobian(mesh);
	} else if (task == "jacobian_T") {
		test_jacobian_T(mesh);
	} else if (task == "test") {
		test_barvortex(mesh, verbose, time);
	} else if (task == "L") {
		test_barvortex_L(mesh);
	} else if (task == "L2") {
		test_barvortex_L2(mesh);
	} else if (task == "LT") {
		test_barvortex_LT(mesh);
	} else if (task == "dymnikov_196") {
		test_dymnikov_196 (mesh);
	} else if (task == "kornev1") {
		test_kornev1(mesh);
	} else if (task == "laplace_LT") {
		test_laplace_LT(mesh);
	} else {
		usage(argv[0]);
	}
	return 0;
}
