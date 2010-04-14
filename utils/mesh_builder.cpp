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

#include <math.h>
#include <assert.h>

#include <map>
#include <set>
#include <list>

#include "mesh_builder.h"

using namespace std;

void sphere_orto_projector(Point & p)
{
	// ортогональная проекция на сферу
	double & x = p.x;
	double & y = p.y;
	double & z = p.z;
	double l = sqrt(x * x + y * y + z * z);
	x /= l;
	y /= l;
	z /= l;
}

void sphere_z_projector(Point & p)
{
	// неортогональная проекция (с сохранением z)
	double & x = p.x;
	double & y = p.y;
	double & z = p.z;
	double l = sqrt(1.0 - z * z) / sqrt(x * x + y * y);
	x *= l;
	y *= l;
}

void flat_projector(Point & p)
{
	// do nothing
}

static void divide(vector < Triangle > & mesh, vector < Point > & points, surface_projector project)
{
	vector < Triangle > new_mesh;
	vector < Point   > new_points;

	int old_size = (int)points.size();
	vector < map < int, int > > already(old_size);

	new_points.insert(new_points.end(), points.begin(), points.end());

	for (vector < Triangle >::iterator it = mesh.begin(); it != mesh.end(); ++it)
	{
		Point a = points[it->v1];
		Point b = points[it->v2];
		Point c = points[it->v3];

		int p1, p2, p3;

		if (already[it->v1].find(it->v2) != already[it->v1].end()) {
			p1 = already[it->v1][it->v2];
		} else {
			p1 = (int)new_points.size();
			Point v1 = (a + b) / 2;
			project(v1);
			new_points.push_back(v1);

			already[it->v1][it->v2] = already[it->v2][it->v1] = p1;
		}

		if (already[it->v1].find(it->v3) != already[it->v1].end()) {
			p2 = already[it->v1][it->v3];
		} else {
			p2 = (int)new_points.size();
			Point v2 = (a + c) / 2;
			project(v2);
			new_points.push_back(v2);

			already[it->v1][it->v3] = already[it->v3][it->v1] = p2;
		}

		if (already[it->v2].find(it->v3) != already[it->v2].end()) {
			p3 = already[it->v2][it->v3];
		} else {
			p3 = (int)new_points.size();
			Point v3 = (b + c) / 2;
			project(v3);
			new_points.push_back(v3);

			already[it->v2][it->v3] = already[it->v3][it->v2] = p3;
		}

		new_mesh.push_back(Triangle(it->v1, p1, p2));
		new_mesh.push_back(Triangle(it->v3, p2, p3));
		new_mesh.push_back(Triangle(it->v2, p3, p1));
		new_mesh.push_back(Triangle(p1, p2, p3));
	}
	mesh.swap(new_mesh);
	points.swap(new_points);
}

void iterate_mesh(vector < Triangle > & mesh, vector < Point > & points, int iterations, surface_projector project)
{
	for (int i = 0; i < iterations; ++i) {
		divide(mesh, points, project);
	}
}

void normalize_mesh(vector < Triangle > & mesh, vector < Point > & points, surface_projector project)
{
	for (vector < Point >::iterator it = points.begin();
		it != points.end(); ++it)
	{
		project(*it);
	}
}

void filter_mesh(vector < Triangle > & mesh,
				 vector < Point > & points,
				 vector < int > & boundary,
				 filter_condition ok)
{
	if (!ok) return;

	vector < Point > new_points;
	vector < Triangle > new_mesh;
	vector < int > nums(points.size());

	for (size_t i = 0; i < points.size(); ++i) {
		nums[i] = -1;
	}

	set < int > bnd;
	for (size_t i = 0; i < mesh.size(); ++i) {
		Triangle & t = mesh[i];

		Point & p1 = points[t.v1];
		Point & p2 = points[t.v2];
		Point & p3 = points[t.v3];

		if (ok(p1) && nums[t.v1] < 0) {
			nums[t.v1] = (int)new_points.size();
			new_points.push_back(p1);
		}

		if (ok(p2) && nums[t.v2] < 0) {
			nums[t.v2] = (int)new_points.size();
			new_points.push_back(p2);
		}

		if (ok(p3) && nums[t.v3] < 0) {
			nums[t.v3] = (int)new_points.size();
			new_points.push_back(p3);
		}

		if (ok(p1) && ok(p2) && ok(p3)) {
			t.v1 = nums[t.v1];
			t.v2 = nums[t.v2];
			t.v3 = nums[t.v3];
			new_mesh.push_back(t);
		} else {
			if (ok(p1)) {
				bnd.insert(nums[t.v1]);
			}

			if (ok(p2)) {
				bnd.insert(nums[t.v2]);
			}

			if (ok(p3)) {
				bnd.insert(nums[t.v3]);
			}
		}
	}

	boundary.insert(boundary.end(), bnd.begin(), bnd.end());

	new_mesh.swap(mesh);
	new_points.swap(points);
}

class Graph
{
	vector < set < int > > data;
	int n;

public:
	Graph(): n(0) {}

	void resize(int n1)
	{
		n = n1;
		data.resize(n);
	}

	void adj(int i, int j)
	{
		data[i].insert(j);
	}

	void symmetric_test()
	{
		for (int i = 0; i < n; ++i)
		{
			for (set < int >::iterator n = data[i].begin(); n != data[i].end(); ++n)
			{
				if (data[*n].find(i) == data[*n].end()) {
					assert(0);
				}
			}
		}
	}

	void rcm()
	{
		vector < multimap < int, int > > sorted1; // node->[neighbors], neighbor: degree->node
		vector < list < int > > sorted;
		vector < int > visited;
		vector < int > order; // old_id -> new_id
		sorted1.resize(n);
		sorted.resize(n);
		order.resize(n);

		symmetric_test();

		for (int i = 0; i < n; ++i)
		{
			for (set < int >::iterator n = data[i].begin(); n != data[i].end(); ++n)
			{
				sorted1[i].insert(make_pair((int)data[*n].size(), *n));
			}

			order[i] = -1;
		}

		for (int i = 0; i < n; ++i)
		{
			for (multimap < int, int >::iterator n = sorted1[i].begin(); 
				n != sorted1[i].end(); ++n)
			{
				sorted[i].push_back(n->second);
			}
		}

		int v = 0;
		int k = 0;
		list < int > vertices;

again:
		vertices.push_back(v);

		do
		{
			// print node

			v = vertices.front();
			vertices.pop_front();

			if (order[v] != -1) {
				continue;
			}

			order[v]   = k++;
			list < int > & neighbours = sorted[v];
			list < int > unvisited;

			for (list < int >::iterator it = neighbours.begin();
					it != neighbours.end(); ++it)
			{
				if (order[*it] < 0) {
					unvisited.push_back(*it);
				}
			}

			unvisited.swap(neighbours);
			vertices.insert(vertices.end(), neighbours.begin(), neighbours.end());
		} while (!vertices.empty());

		for (int i = 0; i < n; ++i)
		{
			if (order[i] == -1) {
				v = i;
				goto again;
			}
		}

		vector < set < int > > new_data;
		new_data.resize(n);

		for (int i = 0; i < n; ++i) {
			for (set < int > ::iterator it = data[i].begin(); it != data[i].end(); ++it)
			{
				new_data[order[i]].insert(order[*it]);
			}
		}

		new_data.swap(data);
	}

	void print(FILE * f, int w, int h)
	{
		int blocks_x = n / ((n + w - 1) / w);
		int blocks_y = n / ((n + h - 1) / h);
		vector < bool > p(blocks_x * blocks_y);

		fprintf(f, "P1\n");
		fprintf(f, "%d %d\n", blocks_x, blocks_y);

#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (set < int > ::iterator it = data[i].begin();
					it != data[i].end(); ++it)
			{
				int y = i * blocks_y / n;
				int x = *it * blocks_x / n;
				p[y * blocks_x + x] = true;
			}
		}

		for (int y = 0; y < blocks_y; ++y) {
			for (int x = 0; x < blocks_x; ++x) {
				fprintf(f, "%d ", (int)p[y * blocks_x + x]);
			}
			fprintf(f, "\n");
		}
	}
};

struct Mesh
{
	typedef std::vector < Triangle > triangles_t;///<triangles container
	typedef std::vector < Point >     points_t;  ///<points container
	typedef std::vector < int > points_flags_t;  ///<points properties container

	triangles_t tr; ///<triangles array
	points_t ps;    ///<points array

	/**
	 * Properties array/
	 *  - 0 - inner point
	 *  - 1 - boundary point
	 */
	points_flags_t ps_flags;

	/**
	 * mapping: point -> triangle in point.
	 */
	std::vector < std::vector < int > > adj;
	/**
	 * Sequence numbers of inner points.
	 */
	std::vector < int > inner;
	/**
	 * Sequence numbers of boundary points.
	 */
	std::vector < int > outer;
	/**
	 * mapping: global point number -> inner point number or outer point number.
	 */
	std::vector < int > p2io;

	void load(std::vector < Triangle > & tri,
		vector < Point > & points,
		vector < int > & boundary)
	{
		tr = tri;
		ps = points;

		adj.resize(ps.size());
		ps_flags.resize(ps.size());
		p2io.resize(ps.size());

		for (int tid = 0; tid < (int)tri.size(); ++tid) {
			Triangle & t = tri[tid];
			adj[t.v1].push_back(tid);
			adj[t.v2].push_back(tid);
			adj[t.v3].push_back(tid);
		}

		for (int i = 0; i < (int)boundary.size(); ++i) {
			ps_flags[boundary[i]] = 1;
		}

		for (int i = 0; i < (int)ps.size(); ++i)
		{
			if (ps_flags[i] == 0)
			{
				p2io[i] = (int) inner.size();
				inner.push_back (i);
			}
			else if (ps_flags[i] == 1)
			{
				p2io[i] = (int) outer.size();
				outer.push_back (i);
			}
		}
	}

	void generate_graph(Graph & g)
	{
		int rs  = (int) inner.size();    // размерность
		g.resize(rs);

		for (int i = 0; i < rs; ++i)   // номер строки
		{
			// по внутренним точкам
			int p = inner[i];

			for (int tk = 0; tk < (int)adj[p].size(); ++tk)
			{
				// по треугольникам в точке
				int trk_i = adj[p][tk];
				const Triangle & t = tr[trk_i];

				if (ps_flags[t.v1] == 0) {
					g.adj(i, p2io[t.v1]);
				}

				if (ps_flags[t.v2] == 0) {
					g.adj(i, p2io[t.v2]);
				}

				if (ps_flags[t.v2] == 0) {
					g.adj(i, p2io[t.v3]);
				}
			}
		}
	}
};

void reorder_mesh(std::vector < Triangle > & tri,
	vector < Point > & points,
	vector < int > & boundary)
{
	Mesh mesh;
	mesh.load(tri, points, boundary);
}

void vizualize_adj(vector < Triangle > & tri,
	vector < Point > & points,
	vector < int > & boundary, 
	FILE * f,
	int w,
	int h)
{
	if (!f) return;

	Mesh mesh;
	Graph g;
	mesh.load(tri, points, boundary);
	mesh.generate_graph(g);
	g.rcm();
	g.print(f, w, h);
}

