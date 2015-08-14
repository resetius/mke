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
#include <algorithm>
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

template<typename T>
int Mesh::find_point2(T & seen, int pn1, int pn2) {
	if (pn1 > pn2) {
		std::swap(pn1, pn2);
	}
	auto it = seen.find({ pn1, pn2 });
	if (it != seen.end()) {
		return it->second;
	}
	else {
		MeshPoint & p1 = ps[pn1];
		MeshPoint & p2 = ps[pn2];
		MeshPoint p3;
		if (p1.is_boundary() && p2.is_boundary()) {
			p3.flags = MeshPoint::POINT_BOUNDARY;
		}
		p3.pr = (p1.pr + p2.pr) / 2;
		for (int i = 0; i < (int)p1.p.size(); ++i) {
			p3.add((p1.p[i] + p2.p[i]) / 2);
		}
		int pn = (int)ps.size();
		ps.push_back(p3);
		seen[{pn1, pn2}] = pn;
		return pn;
	}
}

template<typename T>
int Mesh::find_point3(T & seen, int pn1, int pn2, int pn3_or_n, int flag) {
	int tn1 = pn1, tn2 = pn2, tn3 = pn3_or_n;
	if (flag == 0) {
		if (tn1 > tn2) {
			std::swap(tn1, tn2);
			if (tn3 == 2) {
				tn3 = 1;
			}
			else {
				tn3 = 2;
			}
		}
	}
	else {
		if (tn1 > tn2) {
			std::swap(tn1, tn2);
		}
		if (tn1 > tn3) {
			std::swap(tn1, tn3);
		}
		if (tn2 > tn3) {
			std::swap(tn2, tn3);
		}
	}
	auto it = seen.find({ tn1, tn2, tn3, flag });
	if (it != seen.end()) {
		return it->second;
	}
	else {
		MeshPoint p4;
		if (flag) {
			MeshPoint & p1 = ps[pn1];
			MeshPoint & p2 = ps[pn2];
			MeshPoint & p3 = ps[pn3_or_n];			
			p4.pr = (p1.pr + p2.pr + p3.pr) / 3;
			for (int i = 0; i < (int)p1.p.size(); ++i) {
				p4.add((p1.p[i] + p2.p[i] + p3.p[i]) / 3);
			}
		}
		else {
			MeshPoint & p1 = ps[pn1];
			MeshPoint & p2 = ps[pn2];
			if (p1.is_boundary() && p2.is_boundary()) {
				p4.flags = MeshPoint::POINT_BOUNDARY;
			}
			double s = pn3_or_n;
			p4.pr = (p1.pr * (3.0 - s) + p2.pr * s) / 3;
			for (int i = 0; i < (int)p1.p.size(); ++i) {
				p4.add((p1.p[i] * (3.0 - s) + p2.p[i] * s) / 3);
			}
		}
		int pn = (int)ps.size();
		ps.push_back(p4);
		seen[{ tn1, tn2, tn3, flag }] = pn;
		return pn;
	}
}

void Mesh::prepare_additional_points() {
	struct Comparator {
		bool operator()(const std::vector<int> & a, const std::vector<int> & b)
		{
			return memcmp(&a[0], &b[0], a.size() * sizeof(int)) < 0;
		}
	};
	std::map<std::vector<int>, int, Comparator> seen;

	//for (int tid = 0; tid < (int)tr.size(); ++tid) {
	//	std::sort(tr[tid].p.begin(), tr[tid].p.end());
	//}

	switch (order)
	{
	case 3:
		for (int tid = 0; tid < (int)tr.size(); ++tid) {
			int pn1 = tr[tid].p[0];
			int pn2 = tr[tid].p[1];
			int pn3 = tr[tid].p[2];

			int pn4 = find_point3(seen, pn1, pn2, 1);
			int pn5 = find_point3(seen, pn1, pn2, 2);
			int pn6 = find_point3(seen, pn2, pn3, 1);
			int pn7 = find_point3(seen, pn2, pn3, 2);
			int pn8 = find_point3(seen, pn1, pn3, 1);
			int pn9 = find_point3(seen, pn1, pn3, 2);

			int pn10 = find_point3(seen, pn1, pn2, pn3, 1);

			assert(ps[pn10].pr == (ps[pn4].pr + ps[pn7].pr) / 2);
			assert(ps[pn6].pr == (ps[pn2].pr * 2 + ps[pn3].pr) / 3);
			assert(ps[pn8].pr == (ps[pn1].pr * 2 + ps[pn3].pr) / 3);
			assert(ps[pn10].pr == (ps[pn6].pr + ps[pn8].pr) / 2);
			assert(ps[pn10].pr == (ps[pn5].pr + ps[pn9].pr) / 2);

			tr[tid].p.push_back(pn4);
			tr[tid].p.push_back(pn5);
			tr[tid].p.push_back(pn6);
			tr[tid].p.push_back(pn7);
			tr[tid].p.push_back(pn8);
			tr[tid].p.push_back(pn9);
			tr[tid].p.push_back(pn10);
		}
		break;
	case 2:
		for (int tid = 0; tid < (int)tr.size(); ++tid) {
			int pn1 = tr[tid].p[0];
			int pn2 = tr[tid].p[1];
			int pn3 = tr[tid].p[2];			
			int pn4 = find_point2(seen, pn1, pn2);
			int pn5 = find_point2(seen, pn2, pn3);
			int pn6 = find_point2(seen, pn3, pn1);

			tr[tid].p.push_back(pn4);
			tr[tid].p.push_back(pn5);
			tr[tid].p.push_back(pn6);
		}
		break;
	case 1:
		break;
	default:
		break;
	}

	size = (int)ps.size();

	adj.resize(size);
	for (int tid = 0; tid < (int)tr.size(); ++tid) {
		for (int i = 0; i < (int)tr[tid].p.size(); ++i) {
			int p = tr[tid].p[i];
			adj[p].push_back(tid);
		}
	}
}

bool Mesh::load(FILE * f, int ord)
{
	int size;
	int lineno = 1;

	order = ord;

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

		Triangle t(n1, n2, n3, ps, convs, z, order);
		tid = (int)tr.size();
		tr.push_back(t);
		//if ((int)adj.size() <= n1) adj.resize(n1 + 1);
		//if ((int)adj.size() <= n2) adj.resize(n2 + 1);
		//if ((int)adj.size() <= n3) adj.resize(n3 + 1);
		//adj[n1].push_back(tid);
		//adj[n2].push_back(tid);
		//adj[n3].push_back(tid);
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

	prepare_additional_points();

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
	for (int i = 0; i < (int)p.size(); ++i) {
		if (p1 == p[i]) {
			return i;
		}
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

	for (int i = 0; i < (int)p.size(); ++i) {
		if (p1 == p[i]) {
			return phik[i];
		}
	}
	abort();
	return *((Polynom*)0);
}

static Polynom L(const Point & p1, const Point & p2) {
	return (P2X - p1.x) * (p2.y - p1.y) - (P2Y - p1.y) * (p2.x - p1.x);
}

Triangle::basis_t Triangle::prepare_basis(int z) const
{
	switch (order) {
	case 3:
		return prepare_basis3(z);
	case 2:
		return prepare_basis2(z);
	default:
		return prepare_basis1(z);
	}
}

Triangle::basis_t Triangle::prepare_basis1(int z) const
{
	basis_t r;
	r.reserve(3);

	// p0
	//r.push_back(
	//	(P2X - x(1, z)) * (y(2, z) - y(1, z)) - 
	//	(P2Y - y(1, z)) * (x(2, z) - x(1, z)));
	//// p1
	//r.push_back(
	//	(P2X - x(0, z)) * (y(2, z) - y(0, z)) - 
	//	(P2Y - y(0, z)) * (x(2, z) - x(0, z)));
	//// p2
	//r.push_back(
	//	(P2X - x(0, z)) * (y(1, z) - y(0, z)) - 
	//	(P2Y - y(0, z)) * (x(1, z) - x(0, z)));

	r.push_back(L(ps[p[1]].p[z], ps[p[2]].p[z]));
	r.push_back(L(ps[p[0]].p[z], ps[p[2]].p[z]));
	r.push_back(L(ps[p[0]].p[z], ps[p[1]].p[z]));

	for (uint i = 0; i < 3; ++i)
	{
		r[i] /= r[i].apply(x(i, z), y(i, z));
	}

	return r;
}

Triangle::basis_t Triangle::prepare_basis2(int z) const
{
	basis_t r;
	basis_t psi;
	psi.push_back(L(ps[p[1]].p[z], ps[p[2]].p[z]));
	psi.push_back(L(ps[p[0]].p[z], ps[p[2]].p[z]));
	psi.push_back(L(ps[p[0]].p[z], ps[p[1]].p[z]));
	
	psi.push_back(L(ps[p[3]].p[z], ps[p[5]].p[z]));
	psi.push_back(L(ps[p[3]].p[z], ps[p[4]].p[z]));
	psi.push_back(L(ps[p[4]].p[z], ps[p[5]].p[z]));	

	r.push_back(psi[0] * psi[3]);
	r.push_back(psi[1] * psi[4]);
	r.push_back(psi[2] * psi[5]);

	r.push_back(psi[0] * psi[1]);
	r.push_back(psi[1] * psi[2]);
	r.push_back(psi[0] * psi[2]);

	for (int i = 0; i < (int)r.size(); ++i) {
		r[i] /= r[i].apply(x(i, z), y(i, z));
	}

	for (int i = 0; i < (int)r.size(); ++i) {
		for (int j = 0; j < (int)p.size(); ++j) {
			double delta = (i == j) ? 1 : 0;
			double v = r[i].apply(x(j, z), y(j, z));
			assert(fabs(v - delta) < 1e-10);
		}
	}
	return r;
}

Triangle::basis_t Triangle::prepare_basis3(int z) const
{
	basis_t r;
	basis_t psi;

	psi.push_back(L(ps[p[1]].p[z], ps[p[2]].p[z]));
	psi.push_back(L(ps[p[0]].p[z], ps[p[2]].p[z]));
	psi.push_back(L(ps[p[0]].p[z], ps[p[1]].p[z]));

	psi.push_back(L(ps[p[3]].p[z], ps[p[7]].p[z]));
	psi.push_back(L(ps[p[4]].p[z], ps[p[8]].p[z]));
	psi.push_back(L(ps[p[4]].p[z], ps[p[5]].p[z]));

	psi.push_back(L(ps[p[3]].p[z], ps[p[6]].p[z]));
	psi.push_back(L(ps[p[6]].p[z], ps[p[8]].p[z]));
	psi.push_back(L(ps[p[5]].p[z], ps[p[7]].p[z]));

	r.push_back(psi[0] * psi[3] * psi[4]);
	r.push_back(psi[1] * psi[5] * psi[6]);
	r.push_back(psi[2] * psi[7] * psi[8]);

	r.push_back(psi[0] * psi[1] * psi[4]);
	r.push_back(psi[0] * psi[1] * psi[6]);
	r.push_back(psi[1] * psi[2] * psi[6]);

	r.push_back(psi[1] * psi[2] * psi[8]);
	r.push_back(psi[0] * psi[2] * psi[4]);
	r.push_back(psi[0] * psi[2] * psi[8]);

	r.push_back(psi[0] * psi[1] * psi[2]);

	for (int i = 0; i < (int)r.size(); ++i) {
		r[i] /= r[i].apply(x(i, z), y(i, z));
	}

	assert(ps[p[9]].pr == (ps[p[0]].pr + ps[p[1]].pr + ps[p[2]].pr) / 3);

	for (int i = 0; i < (int)r.size(); ++i) {
		for (int j = 0; j < (int)p.size(); ++j) {
			double delta = (i == j) ? 1 : 0;
			double v = r[i].apply(x(j, z), y(j, z));
			assert(fabs(v - delta) < 1e-10);
		}
	}

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

void Triangle::proj() {
	for (int i = 0; i < 3; ++i) {
		pp[i] = ps[p[i]].pr;
	}

	Point n = (pp[0] - pp[1]) * (pp[0] - pp[2]);

	if (n.len() < 1e-14) {
		return;
	}
	if (sqrt(n.x*n.x + n.y*n.y) < 1e-14) {
		return;
	}

	n = n / n.len();
	double cosa = fabs(n.x) / sqrt(n.x*n.x + n.y*n.y);
	double cosb = fabs(n.z);
	double a = -sign(n.x*n.y)*acos(cosa);
	n = n.rotate_z(a);
	double b = sign(n.x*n.z)*acos(cosb);
	n = n.rotate_y(b);

	m = Matrix();
	m.rotate_z(a);
	m.rotate_y(b);

	m1 = Matrix();
	m1.rotate_y(-b);
	m1.rotate_z(-a);

	for (int i = 0; i < 3; ++i) {
		pp[i] = pp[i].apply(m);
	}

	assert(fabs(pp[0].z - pp[1].z) < 1e-15 && fabs(pp[1].z - pp[2].z) < 1e-15);
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

	double zdiff = - pp[0].z;
	
	e0.f =
		(X1 - pp[1].x) * (pp[2].y - pp[1].y) -
		(Y1 - pp[1].y) * (pp[2].x - pp[1].x);
	
	e1.f = 
		(X1 - pp[0].x) * (pp[2].y - pp[0].y) - 
		(Y1 - pp[0].y) * (pp[2].x - pp[0].x);

	e2.f = 
		(X1 - pp[0].x) * (pp[1].y - pp[0].y) - 
		(Y1 - pp[0].y) * (pp[1].x - pp[0].x);

	e0.h.hx = m.m[0][0] * X + m.m[0][1] * Y + m.m[0][2] * Z;
	e0.h.hy = m.m[1][0] * X + m.m[1][1] * Y + m.m[1][2] * Z;
	e0.h.hz = m.m[2][0] * X + m.m[2][1] * Y + m.m[2][2] * Z + zdiff;
	e0.h.hx->bind_args({ "x", "y", "z" });
	e0.h.hy->bind_args({ "x", "y", "z" });
	e0.h.hz->bind_args({ "x", "y", "z" });
	e2.h = e1.h = e0.h;

	// TODO: check this transforms
	e0.h1.h1x = m1.m[0][0] * X1 + m1.m[0][1] * Y1 + m1.m[0][2] * (Z1 - zdiff);
	e0.h1.h1y = m1.m[1][0] * X1 + m1.m[1][1] * Y1 + m1.m[1][2] * (Z1 - zdiff);
	e0.h1.h1z = m1.m[2][0] * X1 + m1.m[2][1] * Y1 + m1.m[2][2] * (Z1 - zdiff);
	e0.h1.h1x->bind_args({ "x1", "y1", "z1" });
	e0.h1.h1y->bind_args({ "x1", "y1", "z1" });
	e0.h1.h1z->bind_args({ "x1", "y1", "z1" });
	e2.h1 = e1.h1 = e0.h1;

	e0.g.gx = e1.g.gx = e2.g.gx = convs[z].g.gx;
	e0.g.gy = e1.g.gy = e2.g.gy = convs[z].g.gy;
	e0.g.gz = e1.g.gz = e2.g.gz = convs[z].g.gz;

	convs[z].g1.g1u->bind_args({ "x", "y", "z" });
	convs[z].g1.g1v->bind_args({ "x", "y", "z" });

	e0.g1.g1u = e1.g1.g1u = e2.g1.g1u = convs[z].g1.g1u;
	e0.g1.g1v = e1.g1.g1v = e2.g1.g1v = convs[z].g1.g1v;

	r.push_back(e0);
	r.push_back(e1);
	r.push_back(e2);

	for (int i = 0; i < 3; ++i) {
		r[i].f->bind_args({ "x1", "y1" });
		r[i].f = r[i].f * (1.0/r[i].f->apply({ pp[i].x, pp[i].y })->value());
	}

	return r;
}
