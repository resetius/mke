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

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <vector>
#include <map>

using namespace std;

struct Vector {
	double x;
	double y;
	double z;

	Vector(double x1, double y1, double z1) : x(x1), y(y1), z(z1) {}
	Vector(double * x1) : x(x1[0]), y(x1[1]), z(x1[2]) {}

	void normalize(bool save_z) {
		if (save_z) {
			double l = sqrt(1.0 - z * z) / sqrt(x * x + y * y);
			x *= l;
			y *= l;
		} else {
			double l = sqrt(x * x + y * y + z * z);
			x /= l;
			y /= l;
			z /= l;
		}
	}
};

Vector operator + (const Vector & a, const Vector & b)
{
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector operator / (const Vector & a, double k)
{
	return Vector(a.x / k, a.y / k, a.z / k);
}

struct Triangle {
	int v1;
	int v2;
	int v3;

	Triangle(int v1_, int v2_, int v3_) : v1(v1_), v2(v2_), v3(v3_) 
	{
	}
};

void build_icosahedron(vector < Triangle > & r, vector < Vector > & p)
{
	double t = (1.0 + sqrt(5.0)) / 2.0;

	double points[][3] = {
		{t, 1.0, 0.0},
		{-t, 1.0, 0.0},
		{t, -1.0, 0.0},

		{-t, -1.0, 0.0},
		{1.0, 0.0, t},
		{1.0, 0.0, -t},

		{-1.0, 0.0, t},
		{-1.0, 0.0, -t},
		{0.0, t, 1.0},

		{0.0, -t, 1.0},
		{0.0, t, -1.0},
		{0.0, -t, -1.0},
	};

	int triangles[][3] = {
		{0, 8,  4},
		{1, 10, 7},
		{2, 9,  11},
		{7, 3,  1},

		{0, 5,  10},
		{3, 9,  6},
		{3, 11, 9},
		{8, 6,  4},

		{2, 4, 9},
		{3, 7, 11},
		{4, 2, 0},
		{9, 4, 6},

		{2,  11, 5},
		{0,  10, 8},
		{5,  0,  2},
		{10, 5,  7},

		{1,  6, 8},
		{1,  8, 10},
		{6,  1, 3},
		{11, 7, 5},
	};

	int v = 12;
	int f = 20;
	int i;

	for (i = 0; i < f; ++i) {
		Triangle t1(triangles[i][0], triangles[i][1], triangles[i][2]);
		r.push_back(t1);
	}

	for (i = 0; i < v; ++i) {
		p.push_back(points[i]);
	}
}

void build_hemisphere(vector < Triangle > & r, vector < Vector > & p)
{
	double points[][3] = {
		{0, 1, 0},
		{0.866025403784438, -0.5, 0},
		{-0.866025403784438, -0.5, 0},
		{0, 0.707106781186548, 0.707106781186547},
		{0.612372435695794, -0.353553390593274, 0.707106781186547},
		{-0.612372435695794, -0.353553390593274, 0.707106781186547},
		{0, -0.5, 0},
		{-0.433012701892219, 0.25, 0},
		{0.433012701892219, 0.25, 0},
	};

	int triangles[][3] = {
		{2, 7, 5},
		{7, 3, 0},
		{3, 0, 8},
		{4, 8, 1},
		{4, 6, 1},
		{2, 5, 6},
		{7, 5, 3},
		{3, 8, 4},
		{4, 6, 5},
		{3, 4, 5},
	};
	
	int v = 9;
	int f = 10;
	int i;
	
	for (i = 0; i < f; ++i) {
		Triangle t1(triangles[i][0], triangles[i][1], triangles[i][2]);
		r.push_back(t1);
	}
	
	for (i = 0; i < v; ++i) {
		p.push_back(points[i]);
	}
}

void build_test(vector < Triangle > & r, vector < Vector > & p)
{
	double z = sin(M_PI / 4.0);
	double points[][3] = {
		{-1  , 0  , 0},
		{-0.5, 0.5, 0},
		{ 0  , 1  , 0},
		{ 0.5, 0.5, 0},
		{ 1  ,  0  , 0},

		{-1  , 0  , z},
		{ 0  , 1  , z},
		{ 1  , 0  , z},
	};

	int triangles[][3] = {
		{0, 5, 1},
		{1, 5, 6},
		{1, 6, 2},
		{2, 6, 3},
		{3, 6, 7},
		{3, 7, 4},
	};

	int v = 8;
	int f = 6;
	int i;

	for (i = 0; i < f; ++i) {
		Triangle t1(triangles[i][0], triangles[i][1], triangles[i][2]);
		r.push_back(t1);
	}

	for (i = 0; i < v; ++i) {
		p.push_back(points[i]);
	}
}

void divide(vector < Triangle > & mesh, vector < Vector > & points, bool save_z)
{
	vector < Triangle > new_mesh;
	vector < Vector   > new_points;
	
	int old_size = points.size();
	vector < map < int, int > > already(old_size);

	new_points.insert(new_points.end(), points.begin(), points.end());

	for (vector < Triangle >::iterator it = mesh.begin(); it != mesh.end(); ++it)
	{
		Vector a = points[it->v1];
		Vector b = points[it->v2];
		Vector c = points[it->v3];

		int p1, p2, p3;

		if (already[it->v1].find(it->v2) != already[it->v1].end()) {
			p1 = already[it->v1][it->v2];
		} else {
			p1 = new_points.size();
			Vector v1 = (a + b) / 2; v1.normalize(save_z);
			new_points.push_back(v1);
			
			already[it->v1][it->v2] = already[it->v2][it->v1] = p1;
		}
		
		if (already[it->v1].find(it->v3) != already[it->v1].end()) {
			p2 = already[it->v1][it->v3];
		} else {
			p2 = new_points.size();
			Vector v2 = (a + c) / 2; v2.normalize(save_z);
			new_points.push_back(v2);
			
			already[it->v1][it->v3] = already[it->v3][it->v1] = p2;
		}
		
		if (already[it->v2].find(it->v3) != already[it->v2].end()) {
			p3 = already[it->v2][it->v3];
		} else {
			p3 = new_points.size();
			Vector v3 = (b + c) / 2; v3.normalize(save_z);
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

void iterate_mesh(vector < Triangle > & mesh, vector < Vector > & points, int iterations, bool save_z)
{
	for (int i = 0; i < iterations; ++i) {
		divide(mesh, points, save_z);
	}
}

void normalize_mesh(vector < Triangle > & mesh, vector < Vector > & points, bool save_z)
{
	for (vector < Vector >::iterator it = points.begin(); 
		it != points.end(); ++it)
	{
		it->normalize(save_z);
	}
}

double u(double x, double y, double z)
{
	return asin(z);
}

double v(double x, double y, double z)
{
	double r  = sqrt(x * x + y * y);

	if (r < 1e-15) return 0.0;

	double v1 = acos(x / r);

	if (y < 0) {
		v1 = 2.0 * M_PI - v1;
	}

	return v1;
}

void print_mesh(const vector < Triangle > & mesh, vector < Vector > & points, int type, bool local)
{
	fprintf(stdout, "# points %lu\n", points.size());
	fprintf(stdout, "# zone1 ; zone2 ; zone 3; ... \n");
	if (local) {
		fprintf(stdout, "# u v \n");
	} else {
		fprintf(stdout, "# x y z \n");
	}

	for (vector < Vector >::const_iterator it = points.begin(); 
		it != points.end(); ++it)
	{
		double x = it->x;
		double y = it->y;
		double z = it->z;

		/**
		 * v = lambda [0, 2pi]
		 * u = theta  [-pi/2, pi/2]
		 * x = cos u cos v
		 * y = cos u sin v
		 * z = sin u
		 */
		if (local) {
			fprintf(stdout, "%.16lf %.16lf \n", u(x, y, z), v(x, y, z));
		} else {
			fprintf(stdout, "%.16lf %.16lf %.16lf\n", x, y, z);
		}
	}
	fprintf(stdout, "# triangles %lu \n", mesh.size());
	
	for (vector < Triangle >::const_iterator it = mesh.begin(); 
		it != mesh.end(); ++it)
	{
		fprintf(stdout, "%lu %lu %lu\n", it->v1 + 1, it->v2 + 1, it->v3 + 1);
	}
	
	if (type == 1 || type == 2) {
		fprintf(stdout, "# boundary \n");
		for (size_t i = 0; i < points.size(); ++i) {
			double x = points[i].x;
			double y = points[i].y;
			double z = points[i].z;

			if (fabs(points[i].z) < 1e-15) {
				fprintf(stdout, "%lu \n", i + 1);
			} else if (type == 2) {
				// u <- 0,pi/4
				// v <- 0,pi
				if (fabs(u(x, y, z) - M_PI / 4.0) < 1e-15) {
					fprintf(stdout, "%lu \n", i + 1);
				} else if (fabs(v(x, y, z)) < 1e-15) {
					fprintf(stdout, "%lu \n", i + 1);
				} else if (fabs(v(x, y, z) - M_PI) < 1e-15) {
					fprintf(stdout, "%lu \n", i + 1);
				}
			}
		}
	}
}

void usage(const char * name)
{
	fprintf(stderr, "usage: %s --type [full|half] --coord [local|global] --iter [number]\n", name);
	fprintf(stderr, "--type\n"
			"   full -- fulsphere\n"
			"   half -- hemisphere\n"
			"--coord\n"
			"   local -- (u,v)\n"
			"   global -- (x,y,z)\n"
			"--iter\n"
			"   number -- number of iterations\n");
	exit(1);
}

int main(int argc, char * argv[])
{
	int iters = 0;

	/**
	 * 0 - full sphere
	 * 1 - hemisphere
	 * 2 - test
	 */
	int type  = 0;
	bool local = false;

	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "--h"))
		{
			usage(argv[0]);
		} else if (!strcmp(argv[i], "--type")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			if (!strcmp(argv[i + 1], "full")) {
				type = 0;
			} else if (!strcmp(argv[i + 1], "half")) {
				type = 1;
			} else if (!strcmp(argv[i + 1], "test")) {
				type = 2;
			} else {
				usage(argv[0]);
			}
		} else if (!strcmp(argv[i], "--coord")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			if (!strcmp(argv[i + 1], "local")) {
				local = true;
			} else if (!strcmp(argv[i + 1], "global")) {
				local = false;
			} else {
				usage(argv[0]);
			}
		} else if (!strcmp(argv[i], "--iter")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			iters = atoi(argv[i + 1]);
		}
	}

	fprintf(stderr, "iterations = %d\n", iters);
	fprintf(stderr, "local coords ? %d\n", (int)local);
	fprintf(stderr, "type = %d\n", type);

	vector < Triangle > mesh;
	vector < Vector >  points;

	switch (type) {
	case 1:
		build_hemisphere(mesh, points);
		break;
	case 0:
		build_icosahedron(mesh, points);
		break;
	case 2:
		build_test(mesh, points);
		break;
	default:
		assert(0);
	}

	normalize_mesh(mesh, points, type == 2);
	iterate_mesh(mesh, points, iters, type == 2);
	print_mesh(mesh, points, type, local);
	return 0;
}

