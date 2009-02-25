/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (������� ���������)
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

vector < Polynom > Mesh::elem1(const Triangle & t) const
{
	vector < Polynom > r;

	// p0
	r.push_back((P2X - t.x(1, ps)) * (t.y(2, ps) - t.y(1, ps)) 
		- (P2Y - t.y(1, ps)) * (t.x(2, ps) - t.x(1, ps)));
	// p1
	r.push_back((P2X - t.x(0, ps)) * (t.y(2, ps) - t.y(0, ps)) 
		- (P2Y - t.y(0, ps)) * (t.x(2, ps) - t.x(0, ps)));
	// p2
	r.push_back((P2X - t.x(0, ps)) * (t.y(1, ps) - t.y(0, ps)) 
		- (P2Y - t.y(0, ps)) * (t.x(1, ps) - t.x(0, ps)));

	for (uint i = 0; i < 3; ++i)
	{
		r[i] /= r[i].apply(t.x(i, ps), t.y(i, ps));
	}

	return r;
}

vector < Polynom > Mesh::elem1_inner(const Triangle & t) const
{
	vector < Polynom > r;
	// p0
	if (ps_flags[t.p[0]] == 0) {
		Polynom p = ((P2X - t.x(1, ps)) * (t.y(2, ps) - t.y(1, ps))
				- (P2Y - t.y(1, ps)) * (t.x(2, ps) - t.x(1, ps)));
		p /= p.apply(t.x(0, ps), t.y(0, ps));
		r.push_back(p);
	}

	// p1
	if (ps_flags[t.p[1]] == 0) {
		Polynom p = ((P2X - t.x(0, ps)) * (t.y(2, ps) - t.y(0, ps))
				- (P2Y - t.y(0, ps)) * (t.x(2, ps) - t.x(0, ps)));
		p /= p.apply(t.x(1, ps), t.y(1, ps));
		r.push_back(p);
	}

	// p2
	if (ps_flags[t.p[2]] == 0) {
		Polynom p = ((P2X - t.x(0, ps)) * (t.y(1, ps) - t.y(0, ps))
				- (P2Y - t.y(0, ps)) * (t.x(1, ps) - t.x(0, ps)));
		p /= p.apply(t.x(2, ps), t.y(2, ps));
		r.push_back(p);
	}
	return r;
}

Polynom Mesh::elem1(const Triangle & t, int p) const
{
	if (p == t.p[0]) {
		Polynom p = (P2X - t.x(1, ps)) * (t.y(2, ps) - t.y(1, ps)) 
			- (P2Y - t.y(1, ps)) * (t.x(2, ps) - t.x(1, ps));
		p /= p.apply(t.x(0, ps), t.y(0, ps));
		return p;
	} else if (p == t.p[1]) {
		Polynom p = (P2X - t.x(0, ps)) * (t.y(2, ps) - t.y(0, ps)) 
			- (P2Y - t.y(0, ps)) * (t.x(2, ps) - t.x(0, ps));
		p /= p.apply(t.x(1, ps), t.y(1, ps));
		return p;
	} else if (p == t.p[2]) {
		Polynom p =(P2X - t.x(0, ps)) * (t.y(1, ps) - t.y(0, ps)) 
			- (P2Y - t.y(0, ps)) * (t.x(1, ps) - t.x(0, ps));
		p /= p.apply(t.x(2, ps), t.y(2, ps));
		return p;
	} else {
		assert(0);
	}
	return Polynom(1, 2);
}

void Mesh::load(FILE * f)
{
	int size;
	int a;
	int lineno = 1;

#define _BUF_SZ 32768
	char s[_BUF_SZ];

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
		double x, y, z;

		if (*s == '#')
			break;

		a = sscanf (s, "%lf%lf%lf", &x, &y, &z);
		if (a != 3 && a != 2)
		{
			goto bad;
		}

		ps.push_back (Point (x, y)); lineno++;
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
		int n1, n2, n3, tid;

		if (*s == '#') 
			break;

		if (sscanf (s, "%d%d%d", &n1, &n2, &n3) != 3)
		{
			goto bad;
		}

		//��� ��� ������� � ����� � 1 � �� � 0
		--n1;
		--n2;
		--n3;
		if (n1 >= size || n2 >= size || n3 >= size)
		{
			goto bad;
		}

		if (n1 < 0 || n1 < 0 || n3 < 0) {
			goto bad;
		}

		Triangle t(n1, n2, n3);
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

	return;

bad:
	{
		fprintf(stderr, "bad file format, line = %d\n", lineno);
		exit(1);
	}
}

void print_function(FILE * to, double * ans, const Mesh & m, 
					x_t x, x_t y, x_t z)
{
	fprintf(to, "# points %lu\n", m.ps.size());
	for (uint i = 0; i < m.ps.size(); ++i) {
		double u = m.ps[i].x;
		double v = m.ps[i].y;
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
	fprintf(to, "# end %lu\n");
}

void generate_matrix(Matrix & A, const Mesh & m, integrate_cb_t integrate_cb, void * user_data)
{
	int rs  = m.inner.size();     // �����������

	Timer t;
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		// �� ���������� ������
		int p = m.inner[i];
		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// �� ������������ � �����
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			Polynom phi_i           = m.elem1(trk, p);
			vector < Polynom > phik = m.elem1(trk);
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				int j    = m.p2io[p2]; //����� ���������� �����
				if (m.ps_flags[p2] == 1) {
					; // �������
				} else {
					double a = integrate_cb(phi_i, phik[i0], p2, trk_i, m, user_data);
					A.add(i, j, a);
				}
			}
		}
	}
	fprintf(stderr, "generate_matrix: %lf \n", t.elapsed()); 
}

void generate_right_part(double * b, const Mesh & m, right_part_cb_t right_part_cb, void * user_data)
{
	int rs  = m.inner.size();     // �����������

	Timer t;
#pragma omp parallel for
	for (int i = 0; i < rs; ++i)
	{
		// �� ���������� ������
		int p = m.inner[i];
		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// �� ������������ � �����
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			Polynom phi_i           = m.elem1(trk, p);
			//vector < Polynom > phik = m.elem1_inner(trk);
			
			// ���� ������� �� �������� �� �������, �� ���� ������ ��������� callback
			// ��� ���������� ���������� �������� � ������ �����?
			vector < Polynom > phik = m.elem1(trk);
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				double a = right_part_cb(phi_i, phik[i0], p2, trk_i, m, user_data);
				b[i]    += a;
			}
		}
	}
	fprintf(stderr, "generate_right_part: %lf \n", t.elapsed());
}

void mke_solve(double * Ans, const double * bnd, double * b, Matrix & A, const Mesh & m)
{
	int sz  = m.ps.size();
	int rs  = m.inner.size();     // �����������
	vector < double > x(rs);      // �����

	Timer t;
	fprintf(stderr, "solve %dx%d: \n", rs, rs);
	A.solve(&b[0], &x[0]);
	fprintf(stderr, "mke_solve: %lf \n", t.elapsed());
	mke_p2u(Ans, &x[0], bnd, m);
}

/* ��������� ������� ������� */
void mke_p2u(double * p, const double * u, const double * bnd, const Mesh & m)
{
	int sz  = m.ps.size();

	for (int i = 0; i < sz; ++i) {
		if (m.ps_flags[i] == 1) {
			//�������
			p[i] = bnd[m.p2io[i]];
		} else {
			//����������
			p[i] = u[m.p2io[i]];
		}
	}
}

/* ������� ������� ������� */
void mke_u2p(double * u, const double * p, const Mesh & m)
{
	int j = 0;
	for (int i = 0; i < (int)m.inner.size(); ++i) {
		u[j++] = p[m.inner[i]];
	}
}
