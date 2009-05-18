/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky 
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>

using namespace std;

struct Point {
	double x;
	double y;

	int bnd;

	Point(double x1, double y1) : x(x1), y(y1), bnd(0) {}
	Point(double * x1) : x(x1[0]), y(x1[1]), bnd(0) {}

	Point & operator /= (double k) {
		x /= k;
		y /= k;
		return *this;
	}
};

Point operator + (const Point & a, const Point & b)
{
	return Point (a.x + b.x, a.y + b.y);
}

bool operator < (const Point &a, const Point & b)
{
	if (a.x < b.x) {
		return true;
	} else if (fabs(a.x - b.x) < 1e-15) {
		if (a.y < b.y) {
			return true;
		}
	}
	return false;
}

struct Triangle {
	int v1;
	int v2;
	int v3;

	Triangle(int v1_, int v2_, int v3_) : v1(v1_), v2(v2_), v3(v3_) 
	{
	}
};

void build_rectangle(double x, double y, double w, double h, vector < Triangle > & r, vector < Point > & p)
{
	double points[][2] = {
		{x, y},
		{x, y + h},
		{x + w, y + h},
		{x + w, y},
	};

	int triangles[][3] = {
		{0, 1, 2},
		{0, 2, 3},
	};

	int v = 4;
	int f = 2;
	int i;

	for (i = 0; i < f; ++i) {
		Triangle t1(triangles[i][0], triangles[i][1], triangles[i][2]);
		r.push_back(t1);
	}

	for (i = 0; i < v; ++i) {
		p.push_back(points[i]);
	}
}

void divide(vector < Triangle > & mesh, vector < Point > & points)
{
	vector < Triangle > new_mesh;
	vector < Point    > new_points;
	int old_size = points.size();
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
			p1 = new_points.size();
			Point v1 = a + b; v1 /= 2.0;
			new_points.push_back(v1);

			already[it->v1][it->v2] = already[it->v2][it->v1] = p1;
		}

		if (already[it->v1].find(it->v3) != already[it->v1].end()) {
			p2 = already[it->v1][it->v3];
		} else {
			p2 = new_points.size();
			Point v2 = a + c; v2 /= 2.0;
			new_points.push_back(v2);

			already[it->v1][it->v3] = already[it->v3][it->v1] = p2;
		}

		if (already[it->v2].find(it->v3) != already[it->v2].end()) {
			p3 = already[it->v2][it->v3];
		} else {
			p3 = new_points.size();
			Point v3 = b + c; v3 /= 2.0;
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

void iterate_mesh(vector < Triangle > & mesh, vector < Point > & points, int iterations)
{
	for (int i = 0; i < iterations; ++i) {
		divide(mesh, points);
	}
}

void print_mesh(double x, double y, double w, double h, const vector < Triangle > & mesh, vector < Point > & points)
{
	for (vector < Point >::const_iterator it = points.begin(); 
		it != points.end(); ++it)
	{
		fprintf(stdout, "%.16lf %.16lf\n", it->x, it->y);
	}
	fprintf(stdout, "# triangles\n");
	
	for (vector < Triangle >::const_iterator it = mesh.begin(); 
		it != mesh.end(); ++it)
	{
		fprintf(stdout, "%d %d %d\n", it->v1 + 1, it->v2 + 1, it->v3 + 1);
	}
	fprintf(stdout, "# boundary\n");

	for (size_t i = 0; i < points.size(); ++i) {
		if (fabs(points[i].x - x) < 1e-15 || fabs(points[i].y - y) < 1e-15 
				|| fabs(points[i].x - x - w) < 1e-15 || fabs(points[i].y - y - h) < 1e-15)
		{
			fprintf(stdout, "%lu \n", i + 1);
		}
	}
}

void usage(const char * name)
{
	fprintf(stderr, "usage: %s x y w h iters\n", name);
	exit(1);
}

int main(int argc, char * argv[])
{
	if (argc < 6) {
		usage(argv[0]);
	}

	vector < Triangle > mesh;
	vector < Point    >  points;

	double x  = atof(argv[1]);
	double y  = atof(argv[2]);
	double w  = atof(argv[3]);
	double h  = atof(argv[4]);
	int iters = atoi(argv[5]);
	fprintf(stderr, "x, y, w, h = %lf, %lf, %lf, %lf\n", x, y, w, h);
	fprintf(stderr, "iterations = %d\n", iters);

	build_rectangle(x, y, w, h, mesh, points);
	iterate_mesh(mesh, points, iters);
	print_mesh(x, y, w, h, mesh, points);
}

