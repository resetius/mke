/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2015 Alexey Ozeritsky
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
 * 3. Redistributions in any form must be accompanied by information on
 *    how to obtain complete source code for the Phelm software and any
 *    accompanying software that uses the Phelm software.  The source code
 *    must either be included in the distribution or be available for no
 *    more than the cost of distribution plus a nominal fee, and must be
 *    freely redistributable under reasonable conditions.  For an
 *    executable file, complete source code means the source code for all
 *    modules it contains.  It does not include source code for modules or
 *    files that typically accompany the major components of the operating
 *    system on which the executable file runs.
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


/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 */

#include <map>
#include <string>
#include "mesh.h"

using namespace phelm;
using namespace std;

static map<string, CoordConv> std_convs;
convs_t phelm::std_id_convs;

static bool load_std_convs() {
	{
		CoordConv id;
		id.name = "id";
		id.g.gx = FuncPtr(new Symb("u"));
		id.g.gy = FuncPtr(new Symb("v"));
		id.g.gz = FuncPtr(new Const(0));

		id.g1.g1u = FuncPtr(new Symb("x"));
		id.g1.g1v = FuncPtr(new Symb("y"));
		std_convs[id.name] = id;

		std_id_convs.resize(1);
		std_id_convs[0] = id;
	}

	{
		CoordConv stdsphere1;

		stdsphere1.name = "stdsphere1";

		FuncPtr x(new Symb("x"));
		FuncPtr y(new Symb("y"));
		FuncPtr z(new Symb("z"));

		FuncPtr u(new Symb("u"));
		FuncPtr v(new Symb("v"));
		stdsphere1.g.gx = FuncPtr(new Cos(u))*FuncPtr(new Cos(v));
		stdsphere1.g.gy = FuncPtr(new Cos(u))*FuncPtr(new Sin(v));
		stdsphere1.g.gz = FuncPtr(new Sin(u));

		stdsphere1.g1.g1u = FuncPtr(new ASin(z));
		stdsphere1.g1.g1v = FuncPtr(new ATan2(y, x));

		std_convs[stdsphere1.name] = stdsphere1;
	}

	{
		CoordConv stdsphere2;

		stdsphere2.name = "stdsphere2";

		FuncPtr x(new Symb("x"));
		FuncPtr y(new Symb("y"));
		FuncPtr z(new Symb("z"));

		FuncPtr u(new Symb("u"));
		FuncPtr v(new Symb("v"));
		stdsphere2.g.gx = FuncPtr(new Cos(u))*FuncPtr(new Cos(v))*-1.0;
		stdsphere2.g.gy = FuncPtr(new Cos(u))*FuncPtr(new Sin(v))*-1.0;
		stdsphere2.g.gz = FuncPtr(new Sin(u));

		stdsphere2.g1.g1u = FuncPtr(new ASin(z));
		stdsphere2.g1.g1v = FuncPtr(new ATan2(y, x));

		std_convs[stdsphere2.name] = stdsphere2;
	}

	{
		CoordConv stdsphere3;

		stdsphere3.name = "stdsphere3";

		FuncPtr x(new Symb("x"));
		FuncPtr y(new Symb("y"));
		FuncPtr z(new Symb("z"));

		FuncPtr u(new Symb("u"));
		FuncPtr v(new Symb("v"));
		stdsphere3.g.gx = FuncPtr(new Cos(u))*FuncPtr(new Cos(v));
		stdsphere3.g.gy = FuncPtr(new Cos(u))*FuncPtr(new Sin(v));
		stdsphere3.g.gz = FuncPtr(new Sin(u))*-1.0;

		stdsphere3.g1.g1u = FuncPtr(new ASin(z))*-1.0;
		stdsphere3.g1.g1v = FuncPtr(new ATan2(y, x));

		std_convs[stdsphere3.name] = stdsphere3;
	}

	{
		CoordConv stdsphere4;

		stdsphere4.name = "stdsphere4";

		FuncPtr x(new Symb("x"));
		FuncPtr y(new Symb("y"));
		FuncPtr z(new Symb("z"));

		FuncPtr u(new Symb("u"));
		FuncPtr v(new Symb("v"));
		stdsphere4.g.gx = FuncPtr(new Cos(u))*FuncPtr(new Cos(v));
		stdsphere4.g.gy = FuncPtr(new Cos(u))*FuncPtr(new Sin(v))*-1.0;
		stdsphere4.g.gz = FuncPtr(new Sin(u));

		stdsphere4.g1.g1u = FuncPtr(new ASin(z));
		stdsphere4.g1.g1v = FuncPtr(new ATan2(y, x))*-1.0;

		std_convs[stdsphere4.name] = stdsphere4;
	}

	return true;
}

static bool loaded = load_std_convs();

bool Mesh::load(FILE * f)
{
	int size;
	int lineno = 1;

#define _BUF_SZ 32768
	char s[_BUF_SZ];

	fprintf(stderr, "#loading mesh ... \n");

	fgets(s, _BUF_SZ - 1, f);
	lineno++;

	//skip comments
	const char * predef = "#predef: zone";
	do
	{
		if (strstr(s, predef) == s) {
			int zone;
			char name[32768];
			if (sscanf(s + strlen(predef), "%d %s", &zone, name) == 2) {
				fprintf(stderr, "#zone %d conv: '%s'\n", zone, name);
				if (convs.size() < zone) {
					convs.resize(zone);
				}
				if (std_convs.find(name) == std_convs.end()) {
					fprintf(stderr, "#conv '%s' not found\n", name);
					exit(-1);
				}
				convs[zone - 1] = std_convs[name];
			}
		}
		if (*s != '#')
		{
			break;
		}
		lineno++;
	} while (fgets(s, _BUF_SZ - 1, f));

	// points
	do
	{
		double x, y, z;
		const char * sep = ";";
		char * str;

		if (*s == '#')
			break;

		MeshPoint p;
		bool pr_isset = false;

		for (str = strtok(s, sep); str; str = strtok(0, sep))
		{
			x = 0;
			y = 0;

			int i = sscanf(str, "%lf%lf%lf", &x, &y, &z);
			switch (i) {
			case 3:
				p.pr = Point(x, y, z);
				pr_isset = true;
				break;
			case 2:
				p.add(Point(x, y));
				break;
			default:
				abort();
			}
		}

		if (!pr_isset) {
			p.pr = Point(p.x(0), p.y(0), 0);
		}

		ps.push_back(p);
		lineno++;
	} while (fgets(s, _BUF_SZ - 1, f));

	size = (int)ps.size();

	if (!fgets(s, _BUF_SZ - 1, f))
	{
		goto bad;
	}
	lineno++;

	// triangles
	do
	{
		int n1, n2, n3, tid, z = 1;

		if (*s == '#')
			break;

		if (sscanf(s, "%d%d%d ; %d", &n1, &n2, &n3, &z) < 3)
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

		if (n1 < 0 || n1 < 0 || n3 < 0)
		{
			goto bad;
		}

		Triangle t(n1, n2, n3, ps, convs, z);
		tid = (int)tr.size();
		tr.push_back(t);
		if ((int)adj.size() <= n1) adj.resize(n1 + 1);
		if ((int)adj.size() <= n2) adj.resize(n2 + 1);
		if ((int)adj.size() <= n3) adj.resize(n3 + 1);
		adj[n1].push_back(tid);
		adj[n2].push_back(tid);
		adj[n3].push_back(tid);
		lineno++;
	} while (fgets(s, _BUF_SZ - 1, f));

	if (!fgets(s, _BUF_SZ - 1, f))
	{
		goto make_inner;
	}
	lineno++;

	// boundary
	do
	{
		int n;
		if (*s == '#')
			break;

		if (sscanf(s, "%d", &n) != 1)
		{
			goto bad;
		}

		--n;

		if (n >= size || n < 0)
		{
			goto bad;
		}

		ps[n].flags = MeshPoint::POINT_BOUNDARY;
		lineno++;
	} while (fgets(s, _BUF_SZ - 1, f));

make_inner:
	fprintf(stderr, "#done ... \n");
	fprintf(stderr, "#preparing mesh ... \n");
	inner.reserve(ps.size());
	outer.reserve(ps.size());
	p2io.resize(ps.size());

	for (uint i = 0; i < ps.size(); ++i)
	{
		if (ps[i].is_regular())
		{
			p2io[i] = (int)inner.size();
			inner.push_back(i);
		}
		else if (ps[i].is_boundary())
		{
			p2io[i] = (int)outer.size();
			outer.push_back(i);
		}
		else
		{
			abort();
		}
	}

	prepare();
	fprintf(stderr, "#done ... \n");

	return true;

bad:
	{
		fprintf(stderr, "#bad file format, line = %d\n", lineno);
		return false;
	}
}

int Triangle::point_number(int p1) const
{
	if (p1 == p[0])
	{
		return 0;
	}
	else if (p1 == p[1])
	{
		return 1;
	}
	else if (p1 == p[2])
	{
		return 2;
	}
	abort();
	return -1;
}


const Triangle::basis_t & Triangle::elem1(int zone) const
{
	if (zone < 0) {
		zone = z;
	}

	if ((int)phik.size() < zone + 1) {
		phik.resize((size_t)(zone + 1));
	}

	if (phik[zone].empty()) {
		phik[zone] = prepare_basis(zone);
	}

	return phik[zone];
}

const Polynom & Triangle::elem1(int p1, int zone) const
{
	const basis_t & phik = elem1(zone);

	if (p1 == p[0])
	{
		return phik[0];
	} 
	else if (p1 == p[1])
	{
		return phik[1];
	} 
	else if (p1 == p[2])
	{
		return phik[2];
	}
	abort();
	return *((Polynom*)0);
}

Triangle::basis_t Triangle::prepare_basis(int z) const
{
	basis_t r;
	r.reserve(3);

	// p0
	r.push_back((P2X - x(1, z)) * (y(2, z) - y(1, z))
		- (P2Y - y(1, z)) * (x(2, z) - x(1, z)));
	// p1
	r.push_back((P2X - x(0, z)) * (y(2, z) - y(0, z))
		- (P2Y - y(0, z)) * (x(2, z) - x(0, z)));
	// p2
	r.push_back((P2X - x(0, z)) * (y(1, z) - y(0, z))
		- (P2Y - y(0, z)) * (x(1, z) - x(0, z)));

	for (uint i = 0; i < 3; ++i)
	{
		r[i] /= r[i].apply(x(i, z), y(i, z));
	}

	/*
	x[0] = X (0, ps);
	y[0] = Y (0, ps);
	x[1] = X (1, ps);
	y[1] = Y (1, ps);
	x[2] = X (2, ps);
	y[2] = Y (2, ps);
	*/
	return r;
}

const Triangle::NewElem & Triangle::new_elem1(int p1, int zone) const
{
	auto & phik = new_elem1(zone);

	if (p1 == p[0])
	{
		return phik[0];
	}
	else if (p1 == p[1])
	{
		return phik[1];
	}
	else if (p1 == p[2])
	{
		return phik[2];
	}
	abort();
	return *((NewElem*)0);
}

const std::vector<Triangle::NewElem> & Triangle::new_elem1(int zone) const
{
	if (zone < 0) {
		zone = z;
	}

	if ((int)newphi.size() < zone + 1) {
		newphi.resize((size_t)(zone + 1));
	}

	if (newphi[zone].empty()) {
		newphi[zone] = prepare_new_basis(zone);
	}

	return newphi[zone];
}

static double sign(double a) {
	if (a < 0) {
		return -1;
	}
	else {
		return 1;
	}
}

std::vector<Triangle::NewElem> Triangle::prepare_new_basis(int z) const
{
	std::vector<NewElem> r;
	FuncPtr X1(new Symb("x1"));
	FuncPtr Y1(new Symb("y1"));
	FuncPtr Z1(new Symb("z1"));

	FuncPtr X(new Symb("x"));
	FuncPtr Y(new Symb("y"));
	FuncPtr Z(new Symb("z"));

	NewElem e0, e1, e2;
	Point p0 = ps[p[0]].pr;
	Point p1 = ps[p[1]].pr;
	Point p2 = ps[p[2]].pr;

	Point n = (p0 - p1) * (p0 - p2);
	n = n / n.len();
	double cosa = fabs(n.x) / sqrt(n.x*n.x + n.y*n.y);
	double cosb = fabs(n.z);
	double a = -sign(n.x*n.y)*acos(cosa);
	n = n.rotate_z(a);
	double b =  sign(n.x*n.z)*acos(cosb);
	n = n.rotate_y(b);

	Matrix m;
	m.rotate_z(a);
	m.rotate_y(b);

	Matrix m1;
	m1.rotate_y(-b);
	m1.rotate_z(-a);

	p0 = p0.apply(m);
	p1 = p1.apply(m);
	p2 = p2.apply(m);

	double zdiff = - p0.z;
	
	e0.f =
		(X1 - p1.x) * (p2.y - p1.y) -
		(Y1 - p1.y) * (p2.x - p1.x);
	
	e1.f = 
		(X1 - p0.x) * (p2.y - p0.y) - 
		(Y1 - p0.y) * (p2.x - p0.x);

	e2.f = 
		(X1 - p0.x) * (p1.y - p0.y) - 
		(Y1 - p0.y) * (p1.x - p0.x);

	e0.h1.h1x = m.m[0][0] * X + m.m[0][1] * Y + m.m[0][2] * Z;
	e0.h1.h1y = m.m[1][0] * X + m.m[1][1] * Y + m.m[1][2] * Z;
	e0.h1.h1z = m.m[2][0] * X + m.m[2][1] * Y + m.m[2][2] * Z + zdiff;
	e2.h1 = e1.h1 = e0.h1;

	e0.h.hx = m1.m[0][0] * X1 + m1.m[0][1] * Y1 + m1.m[0][2] * Z1;
	e0.h.hy = m1.m[1][0] * X1 + m1.m[1][1] * Y1 + m1.m[1][2] * Z1;
	e0.h.hz = m1.m[2][0] * X1 + m1.m[2][1] * Y1 + m1.m[2][2] * Z1 - zdiff;
	e2.h = e1.h = e0.h;

	e0.g.gx = e1.g.gx = e2.g.gx = convs[z].g.gx;
	e0.g.gy = e1.g.gy = e2.g.gy = convs[z].g.gy;
	e0.g.gz = e1.g.gz = e2.g.gz = convs[z].g.gz;

	e0.g1.g1u = e1.g1.g1u = e2.g1.g1u = convs[z].g1.g1u;
	e0.g1.g1v = e1.g1.g1v = e2.g1.g1v = convs[z].g1.g1v;

	r.push_back(e0);
	r.push_back(e1);
	r.push_back(e2);
	return r;
}
