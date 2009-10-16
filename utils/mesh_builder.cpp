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
#include <map>
#include <set>

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
