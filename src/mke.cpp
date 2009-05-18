/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (Алексей Озерицкий)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "mke.h"
#include "util.h"
#include "solver.h"

using namespace std;

bool Mesh::load(FILE * f)
{
	int size;
	int lineno = 1;

#define _BUF_SZ 32768
	char s[_BUF_SZ];

	fprintf(stderr, "loading mesh ... \n");

	fgets (s, _BUF_SZ - 1, f); lineno++;

	//skip comments
	do
	{
		if (*s != '#') {
			break;
		}
		lineno ++;
	}
	while (fgets(s, _BUF_SZ - 1, f));

	// points
	do
	{
		double x, y;
		const char * sep = ";";
		char * str;

		if (*s == '#')
			break;

		MeshPoint p;

		for (str = strtok(s, sep); str; str = strtok(0, sep))
		{
			x = 0; y = 0;
			sscanf (str, "%lf%lf", &x, &y);
			p.add(Point(x, y));
		}

		ps.push_back (p); lineno++;
	}
	while (fgets (s, _BUF_SZ - 1, f));

	size = (int) ps.size();

	ps_flags.resize(size);
	if (!fgets (s, _BUF_SZ - 1, f)) {
		goto bad;
	}
	lineno++;

	// triangles
	do
	{
		int n1, n2, n3, tid, z = 1;

		if (*s == '#') 
			break;

		if (sscanf (s, "%d%d%d ; %d", &n1, &n2, &n3, &z) < 3)
		{
			goto bad;
		}

		//так как индексы в файле с 1 а не с 0
		--n1;
		--n2;
		--n3;
		--z;
		if (n1 >= size || n2 >= size || n3 >= size)
		{
			goto bad;
		}

		if (n1 < 0 || n1 < 0 || n3 < 0) {
			goto bad;
		}

		Triangle t(n1, n2, n3, z);
		tid = tr.size();
		tr.push_back (t);
		if ((int)adj.size() <= n1) adj.resize(n1 + 1);
		if ((int)adj.size() <= n2) adj.resize(n2 + 1);
		if ((int)adj.size() <= n3) adj.resize(n3 + 1);
		adj[n1].push_back(tid);
		adj[n2].push_back(tid);
		adj[n3].push_back(tid); lineno++;
	}
	while (fgets (s, _BUF_SZ - 1, f) );

	if (!fgets (s, _BUF_SZ - 1, f)) {
		goto make_inner;
	}
	lineno++;

	// boundary
	do
	{
		int n;
		if (*s == '#')
			break;

		if (sscanf(s, "%d", &n) != 1) {
			goto bad;
		}

		--n;

		if (n >= size || n < 0) {
			goto bad;
		}

		ps_flags[n] = 1; lineno++;
	}
	while (fgets (s, _BUF_SZ - 1, f) );

make_inner:
	fprintf(stderr, "done ... \n");
	fprintf(stderr, "preparing mesh ... \n");
	inner.reserve(ps.size());
	outer.reserve(ps.size());
	p2io.resize(ps.size());

	for (uint i = 0; i < ps.size(); ++i) {
		if (ps_flags[i] == 0) {
			p2io[i] = inner.size();
			inner.push_back(i);
		} else if (ps_flags[i] == 1) {
			p2io[i] = outer.size();
			outer.push_back(i);
		}
	}

	prepare();
	fprintf(stderr, "done ... \n");

	return true;

bad:
	{
		fprintf(stderr, "bad file format, line = %d\n", lineno);
		return false;
	}
}

void Mesh::prepare()
{
	triangles_t::iterator it, b = tr.begin(), e = tr.end();
	for (it = b; it != e; ++it) {
		it->prepare(ps);
	}
}

void Mesh::info()
{
	fprintf(stderr, "points: %lu\n", ps.size());
	fprintf(stderr, "inner points: %lu\n", inner.size());
	fprintf(stderr, "outer points: %lu\n", outer.size());
	fprintf(stderr, "triangles: %lu\n", tr.size());
}

void print_inner_function(FILE * to, double * ans, const Mesh & m, 
					x_t x, x_t y, x_t z)
{
	fprintf(to, "# points %lu\n", m.inner.size());
	for (uint i = 0; i < m.inner.size(); ++i) {
		int p = m.inner[i];
		double u = m.ps[p].x();
		double v = m.ps[p].y();
		double f = ans[p];

		if (x) {
			fprintf(to, "%.16lf ", x(u, v));
		} else {
			fprintf(to, "%.16lf ", u);
		}

		if (y) {
			fprintf(to, "%.16lf ", y(u, v));
		} else {
			fprintf(to, "%.16lf ", v);
		}

		if (z) {
			fprintf(to, "%.16lf ", z(u, v));
		} else {
			fprintf(to, "%.16lf ", 0.0);
		}

		fprintf(to, "%.16lf %.16lf %.16lf \n", f, u, v);
	}

	fprintf(to, "# triangles %lu\n", m.tr.size());
	for (uint i = 0; i < m.tr.size(); ++i) {
		int p1 = m.tr[i].p[0];
		int p2 = m.tr[i].p[1];
		int p3 = m.tr[i].p[2];

		if (m.ps_flags[p1] == 0
			&& m.ps_flags[p2] == 0
			&& m.ps_flags[p3] == 0) 
		{
			fprintf(to, "%d %d %d\n", 
				m.p2io[p1] + 1, 
				m.p2io[p2] + 1, 
				m.p2io[p3] + 1);
		}
	}
	fprintf(to, "# end \n");
}

void print_inner_function(const char * to, double * ans, const Mesh & m, 
					x_t x, x_t y, x_t z)
{
	FILE * f = fopen(to, "wb");
	if (f) {
		print_inner_function(f, ans, m, x, y, z);
		fclose(f);
	}
}

void print_function(FILE * to, double * ans, const Mesh & m, 
					x_t x, x_t y, x_t z)
{
	fprintf(to, "# points %lu\n", m.ps.size());
	for (uint i = 0; i < m.ps.size(); ++i) {
		double u = m.ps[i].x();
		double v = m.ps[i].y();
		double f = ans[i];

		if (x) {
			fprintf(to, "%.16lf ", x(u, v));
		} else {
			fprintf(to, "%.16lf ", u);
		}

		if (y) {
			fprintf(to, "%.16lf ", y(u, v));
		} else {
			fprintf(to, "%.16lf ", v);
		}

		if (z) {
			fprintf(to, "%.16lf ", z(u, v));
		} else {
			fprintf(to, "%.16lf ", 0.0);
		}

		fprintf(to, "%.16lf %.16lf %.16lf \n", f, u, v);
	}

	fprintf(to, "# triangles %lu\n", m.tr.size());
	for (uint i = 0; i < m.tr.size(); ++i) {
		fprintf(to, "%d %d %d\n", 
			m.tr[i].p[0] + 1, 
			m.tr[i].p[1] + 1, 
			m.tr[i].p[2] + 1);
	}
	fprintf(to, "# end \n");
}

void print_function(const char * fname, double * ans, const Mesh & m, 
					x_t x, x_t y, x_t z)
{
	FILE * f = fopen(fname, "wb");
	if (f) {
		print_function(f, ans, m, x, y, z);
		fclose(f);
	}
}

void generate_matrix(Matrix & A, const Mesh & m, integrate_cb_t integrate_cb, void * user_data)
{
	int rs  = m.inner.size();     // размерность

	Timer t;
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) { // номер строки
		// по внутренним точкам
		int p = m.inner[i];
		;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				int j    = m.p2io[p2]; // номер внутренней точки
				                       // номер столбца
				if (m.ps_flags[p2] == 1) {
					; // граница
				} else {
					double a = integrate_cb(phi_i, phik[i0], trk, m, p, p2, user_data);
					A.add(i, j, a);
				}
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_matrix: %lf \n", t.elapsed()); 
#endif
}

void generate_full_matrix(Matrix & A, const Mesh & m, integrate_cb_t integrate_cb, void * user_data)
{
	int sz  = m.ps.size();

	Timer t;
//#pragma omp parallel for
	for (int p = 0; p < sz; ++p) {
		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				double a = integrate_cb(phi_i, phik[i0], trk, m, p, p2, user_data);
				A.add(p, p2, a);
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_full_matrix: %lf \n", t.elapsed()); 
#endif
}



void generate_right_part(double * b, const Mesh & m, right_part_cb_t right_part_cb, void * user_data)
{
	int rs  = m.inner.size();     // размерность

	Timer t;
#pragma omp parallel for
	for (int i = 0; i < rs; ++i)
	{
		// по внутренним точкам
		int p = m.inner[i];
		b[i]  = 0.0;
		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			// если граница не меняется по времени, то надо делать отдельный callback
			// для вычисления постоянной поправки к правой части?
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				double a = right_part_cb(phi_i, phik[i0], trk, m, p, p2, user_data);
				b[i]    += a;
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_right_part: %lf \n", t.elapsed());
#endif
}

void generate_full_right_part(double * b, const Mesh & m, right_part_cb_t right_part_cb, void * user_data)
{
	int sz  = m.ps.size();     // размерность
	vector < Polynom > phik;

	Timer t;
#pragma omp parallel for
	for (int p = 0; p < sz; ++p)
	{
		b[p]  = 0.0;
		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			// если граница не меняется по времени, то надо делать отдельный callback
			// для вычисления постоянной поправки к правой части?
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				double a = right_part_cb(phi_i, phik[i0], trk, m, p, p2, user_data);
				b[p]    += a;
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_right_part: %lf \n", t.elapsed());
#endif
}

/**
 * генерирует матрицу для интеграции краевых условий в правую часть
 * inner.size() x outer.size()
 * в cb передается phi_j где j точка границы
 * phi_i, где i внутренняя точка
 */
void generate_boundary_matrix(Matrix & A, const Mesh & m, right_part_cb_t right_part_cb, void * user_data)
{
	int os = m.outer.size(); // размер границы
	for (int j = 0; j < os; ++j) {
		// по внешним точкам
		int p2 = m.outer[j];

		for (uint tk = 0; tk < m.adj[p2].size(); ++tk) {
			// по треугольникам в точке
			int trk_j = m.adj[p2][tk];
			const Triangle & trk    = m.tr[trk_j];
			const vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_j           = trk.elem1(p2);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p    = m.tr[trk_j].p[i0];

				if (m.ps_flags[p] == 1) {
					;
				} else {
					// p - внутренняя точка
					int i = m.p2io[p];
					double a = right_part_cb(phik[i0], phi_j, trk, m, p, p2, user_data);
					A.add(i, j, a);
				}
			}
		}
	}
}

void mke_solve(double * Ans, const double * bnd, double * b, Matrix & A, const Mesh & m)
{
	int rs  = m.inner.size();     // размерность
	vector < double > x(rs);      // ответ

	mke_solve2(&x[0], b, A, m);
	mke_p2u(Ans, &x[0], bnd, m);
}

void mke_solve2(double * Ans, double * b, Matrix & A, const Mesh & m)
{
	int sz  = m.ps.size();
	int rs  = m.inner.size();     // размерность
	vector < double > x(rs);      // ответ

#ifdef _DEBUG
	Timer t;
	fprintf(stderr, "solve %dx%d: \n", rs, rs);
#endif
	A.solve(Ans, &b[0]);
#ifdef _DEBUG
	fprintf(stderr, "mke_solve: %lf \n", t.elapsed());
#endif
}

/* добавляем краевые условия */
void mke_p2u(double * u, const double * p, const double * bnd, const Mesh & m)
{
	int sz  = m.ps.size();

#pragma omp parallel for
	for (int i = 0; i < sz; ++i) {
		if (m.ps_flags[i] == 1) {
			//внешняя
			if (bnd) {
				u[i] = bnd[m.p2io[i]];
			} else {
				u[i] = 0;
			}
		} else {
			//внутренняя
			u[i] = p[m.p2io[i]];
		}
	}
}

/* убираем краевые условия */
void mke_u2p(double * p, const double * u, const Mesh & m)
{
	int j = 0;
	int rs = (int)m.inner.size();
	for (int i = 0; i < rs; ++i) {
		p[j++] = u[m.inner[i]];
	}
}

void mke_set_bnd(double *u, const double * bnd, const Mesh & m)
{
	int sz  = m.ps.size();
	for (int i = 0; i < sz; ++i) {
		if (m.ps_flags[i] == 1) {
			if (bnd) {
				u[i] = bnd[m.p2io[i]];
			} else {
				u[i] = 0;
			}
		}
	}
}

double generic_scalar_cb(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, const Mesh & m, int, int, void * )
{
	return integrate(phi_i * phi_j, trk, m.ps);
}

double sphere_scalar_cb(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, const Mesh & m, int, int, void * user_data)
{
	return integrate_cos(phi_i * phi_j, trk, m.ps);
}

void convolution(double * ans, const double * u, const double * v, const Mesh & m, scalar_cb_t cb, void * user_data)
{
	int sz  = m.ps.size(); // размерность

#pragma omp parallel for
	for (int i = 0; i < sz; ++i)
	{
		// по всем точкам
		int p = i;
		ans[i] = 0.0;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольникам в точке
			int trk_i               = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int j  = trk.p[i0];
				ans[i] += u[i] * v[j] * cb(phi_i, phik[i0], trk, m, i, j, user_data);
			}
		}
	}
}

double mke_scalar(const double * u, const double * v, const Mesh & m, scalar_cb_t cb, void * user_data)
{
	int sz  = m.ps.size();
	double s = 0.0;
	vector < double > nr(sz);
	convolution(&nr[0], u, v, m, cb, user_data);
#pragma omp parallel for reduction(+:s)
	for (int i = 0; i < sz; ++i) {
		s = s + nr[i];
	}
	return s;
}

void generate_scalar_matrix(Matrix & mat, const Mesh & m, scalar_cb_t cb, void * user_data)
{
	generate_full_matrix(mat, m, cb, user_data);
}

double mke_fast_scalar(const double * u, const double * v, const Mesh & m, Matrix & A)
{
	int sz = m.ps.size();
	vector < double > tmp(sz);
	A.mult_vector(&tmp[0], v);
	return vec_scalar2(u, &tmp[0], sz);
}

double mke_fast_norm(const double * u, const Mesh & m, Matrix & A)
{
	return sqrt(mke_fast_scalar(u, u, m, A));
}

double mke_fast_dist(const double * u, const double * v, const Mesh & m, Matrix & A)
{
	int sz  = m.ps.size(); // размерность
	vector < double > diff(sz);
	vec_diff(&diff[0], u, v, sz);
	return mke_fast_norm(&diff[0], m, A);
}

double mke_norm(const double * u, const Mesh & m, scalar_cb_t cb, void * user_data)
{
	return sqrt(mke_scalar(u, u, m, cb, user_data));
}

double mke_dist(const double * u, const double * v, const Mesh & m, scalar_cb_t cb, void * user_data)
{
	int sz  = m.ps.size(); // размерность
	vector < double > diff(sz);
	vec_diff(&diff[0], u, v, sz);
	return mke_norm(&diff[0], m, cb, user_data);
}

void mke_proj(double * F, const Mesh & mesh, f_xy_t f)
{
	size_t sz = mesh.ps.size();
	for (size_t i = 0; i < sz; ++i)
	{
		const Point & p = mesh.ps[i].p[0];
		F[i] = f(p.x, p.y);
	}
}

void mke_proj_bnd(double * F, const Mesh & m, f_xy_t f)
{
	for (size_t i = 0; i < m.outer.size(); ++i) {
		int p0 = m.outer[i];
		const Point & p = m.ps[p0].p[0];
		F[i] = f(p.x, p.y);
	}
}

void mke_proj_bnd(double * F, const double * F1, const Mesh & m)
{
#pragma omp parallel for
	for (int i = 0; i < (int)m.outer.size(); ++i) {
		int p0 = m.outer[i];
		F[i] = F1[p0];
	}
}

void mke_proj(double * F, const Mesh & mesh, f_xyt_t f, double t)
{
	for (size_t i = 0; i < mesh.ps.size(); ++i)
	{
		F[i] = f(mesh.ps[i].x(), mesh.ps[i].y(), t);
	}
}

void mke_proj_bnd(double * F, const Mesh & m, f_xyt_t f, double t)
{
	for (size_t i = 0; i < m.outer.size(); ++i) {
		int p0 = m.outer[i];
		const Point & p = m.ps[p0].p[0];
		F[i] = f(p.x, p.y, t);
	}
}

