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

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <vector>
#include <map>
#include <set>
#include <string>

#include "mesh_builder.h"

using namespace std;

double x_ (double u, double v)
{
	return cos (u) * cos (v);
}

double y_ (double u, double v)
{
	return cos (u) * sin (v);
}

double z_ (double u, double v)
{
	return sin (u);
}

void build_icosahedron(vector < Triangle > & r, vector < Point > & p)
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

void build_test2(vector < Triangle > & r, vector < Point > & p)
{
	double t = (1.0 + sqrt(5.0)) / 2.0;

	double points[][3] = {
		{t, 1.0, 0.0},   // 0
		{-t, 1.0, 0.0},  // 1
		{t, -1.0, 0.0},  // 2

		{-t, -1.0, 0.0}, // 3
		{1.0, 0.0, t},   // 4
		{1.0, 0.0, -t},  // 5

		{-1.0, 0.0, t},  // 6
		{-1.0, 0.0, -t}, // 7
	};

	int triangles[][3] = {
		{7, 3,  1},
		{4, 2, 0},
		{5,  0,  2},
		{6,  1, 3},
	};

	int v = 8;
	int f = 4;
	int i;

	for (i = 0; i < f; ++i) {
		Triangle t1(triangles[i][0], triangles[i][1], triangles[i][2]);
		r.push_back(t1);
	}

	for (i = 0; i < v; ++i) {
		p.push_back(points[i]);
	}
}

void build_hemisphere(vector < Triangle > & r, vector < Point > & p)
{
	double points[][3] = {
		{0, 1, 0},
		{0.866025403784438, -0.5, 0},
		{-0.866025403784438, -0.5, 0},
		{0, 0.707106781186548,                   z_(M_PI / 3., 0)},
		{0.612372435695794, -0.353553390593274,  z_(M_PI / 3., 0)},
		{-0.612372435695794, -0.353553390593274, z_(M_PI / 3., 0)},
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

void build_test(vector < Triangle > & r, vector < Point > & p)
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

void build_test5(vector < Triangle > & r, vector < Point > & p)
{
	double ploc[][2] = {
		{-M_PI / 10, 0},
		{ M_PI / 10, 0},
		{ M_PI / 10, M_PI / 5},
		{-M_PI / 10, M_PI / 5},
	};

	double points[][3] = {
		{ x_(ploc[0][0], ploc[0][1]), y_(ploc[0][0], ploc[0][1]), z_(ploc[0][0], ploc[0][1]) },
		{ x_(ploc[1][0], ploc[1][1]), y_(ploc[1][0], ploc[1][1]), z_(ploc[1][0], ploc[1][1]) },
		{ x_(ploc[2][0], ploc[2][1]), y_(ploc[2][0], ploc[2][1]), z_(ploc[2][0], ploc[2][1]) },
		{ x_(ploc[3][0], ploc[3][1]), y_(ploc[3][0], ploc[3][1]), z_(ploc[3][0], ploc[3][1]) },
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

bool ok1(Point & v)
{
	if ((fabs(v.z - 1.0) < 1e-10) || (fabs(v.z + 1.0) < 1e-10)) {
		return false;
	}

	return true;
}

bool ok2(Point & v)
{
	if (v.z < -0.9) {
		return false;
	}

	return true;
}

bool ok3(Point & v)
{
	double z = z_(M_PI / 4, 0);
	if (v.z > z) {
		return false;
	}

	return true;
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

void print_mesh(const vector < Triangle > & mesh,
				vector < Point > & points,
				vector < int > & boundary,
				int type, bool local)
{
	fprintf(stdout, "#points:%lu\n", points.size());
	fprintf(stdout, "#triangles:%lu\n", mesh.size());
	fprintf(stdout, "# zone1 ; zone2 ; zone 3; ... \n");
	if (local) {
		fprintf(stdout, "# u = theta  [-pi/2, pi/2]\n"
				"# v = lambda [0, 2pi]\n"
				"# x = cos u cos v \n"
				"# y = cos u sin v \n"
				"# z = sin u \n"
				"# u v \n");
	} else {
		fprintf(stdout, "# x y z \n");
	}

	for (vector < Point >::const_iterator it = points.begin();
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
			double u_ = u(x, y, z), v_ = v(x, y, z);
			// zone 1
			fprintf(stdout, "%.16lf %.16lf ", u_, v_);
			// zone 2
			u_ = u(-x, -y, z); v_ = v(-x, -y, z);
			fprintf(stdout, " ; %.16lf %.16lf ", u_, v_);

			if (type == 0 || type == 1 || type == 4) {
				// zone 3
				u_ = u(x, -z, y);
				v_ = v(x, -z, y);
				fprintf(stdout, " ; %.16lf %.16lf ", u_, v_);

				// zone 4
				u_ = u(x, z, -y);
				v_ = v(x, z, -y);
				fprintf(stdout, " ; %.16lf %.16lf ", u_, v_);
			}

			fprintf(stdout, "\n");
		} else {
			fprintf(stdout, "%.16lf %.16lf %.16lf\n", x, y, z);
		}
	}
	fprintf(stdout, "# triangles %lu \n", mesh.size());

	for (vector < Triangle >::const_iterator it = mesh.begin();
		it != mesh.end(); ++it)
	{
		fprintf(stdout, "%d %d %d ", it->v1 + 1, it->v2 + 1, it->v3 + 1);
		if (type != 5) {
			double x1, y1, z1, x2, y2, z2, x3, y3, z3;
			double u1, v1, u2, v2, u3, v3;

			x1 = points[it->v1].x;
			y1 = points[it->v1].y;
			z1 = points[it->v1].z;

			x2 = points[it->v2].x;
			y2 = points[it->v2].y;
			z2 = points[it->v2].z;

			x3 = points[it->v3].x;
			y3 = points[it->v3].y;
			z3 = points[it->v3].z;

			u1 = u(x1, y1, z1);
			v1 = v(x1, y1, z1);

			u2 = u(x2, y2, z2);
			v2 = v(x2, y2, z2);

			u3 = u(x3, y3, z3);
			v3 = v(x3, y3, z3);

			if (type == 0 || type == 1 || type == 4) {
				if (z1 > 0.9 || z2 > 0.9 || z3 > 0.9) {
					// in zone 4
					fprintf(stdout, " ; 4 ");
				} else if (z1 < -0.9 || z2 < -0.9 || z3 < -0.9) {
					// in zone 3
					fprintf(stdout, " ; 3 ");
				} else
				if (   v1 < M_PI_2 || v2  < M_PI_2 || v2 < M_PI_2
					|| v1 > 1.5 * M_PI || v2  > 1.5 * M_PI || v2 > 1.5 * M_PI)
				{
					// in zone 2
					fprintf(stdout, " ; 2 ");
				}
			} else if (type != 2) {
				if (   v1 < M_PI_2 || v2  < M_PI_2 || v2 < M_PI_2
					|| v1 > 1.5 * M_PI || v2  > 1.5 * M_PI || v2 > 1.5 * M_PI)
				{
					// in zone 2
					fprintf(stdout, " ; 2 ");
				}
			}
		}

		fprintf(stdout, "\n");
	}

	if (type == 1 || type == 2 || type == 0 || type == 5)
	{
		fprintf(stdout, "# boundary \n");
		for (size_t i = 0; i < points.size(); ++i) {
			double x = points[i].x;
			double y = points[i].y;
			double z = points[i].z;

			if (type == 0 &&
				fabs(z + 1.0) < 1e-14)
			{
				// south pole
				fprintf(stdout, "%lu \n", i + 1);
			}

			if (type == 1 || type == 2) {
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
			} else if (type == 5) {
				if (fabs(u(x, y, z) - M_PI / 10) < 1e-14 ||
					fabs(u(x, y, z) + M_PI / 10) < 1e-14)
				{
					fprintf(stdout, "%lu \n", i + 1);
				}
			}
		}
	} else if (!boundary.empty()) {
		fprintf(stdout, "# boundary %lu \n", boundary.size());
		for (size_t i = 0; i < boundary.size(); ++i) {
			fprintf(stdout, "%d \n", boundary[i] + 1);
		}
	}
}

void usage(const char * name)
{
	fprintf(stderr, "usage: %s --type [full|half] --coord [local|global] --iter [number]\n", name);
	fprintf(stderr, "--type\n"
			"   full -- fulsphere\n"
			"   half -- hemisphere\n"
			"   test -- \n"
			"--coord\n"
			"   local -- (u,v)\n"
			"   global -- (x,y,z)\n"
			"--iter\n"
			"   number -- number of iterations\n"
			" \n\n"
			" You can use Vizualizer 3D v3.2 or higher to vizualize this data\n"
			" You must remove boundary section from the output \n"
			" before sending the data to Vizualizer 3D \n"
			" ");
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
				// половинка полусферы с отрезанной крышкой у полюса
				// две области
				type = 2;
			} else if (!strcmp(argv[i + 1], "test2")) {
				// сфера с вырезанными маленькими кружочками у полюсов
				// две области
				type = 3;
			} else if (!strcmp(argv[i + 1], "test3")) {
				// сфера с вырезанной шапкой у одного из полюсов
				// три области
				type = 4;
			} else if (!strcmp(argv[i + 1], "test4")) {
				// квадрат на сфере
				// lambda: [0,      pi/5 ]
				// phi:    [-pi/10, pi/10]
				type = 5;
			} else if (!strcmp(argv[i + 1], "test5")) {
				// полусфера с вырезанной шапкой у полюса
				// две области
				type = 6;
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

	fprintf(stdout, "#cmd:");
	for (int i = 0; i < argc; ++i) {
		fprintf(stdout, "%s ", argv[i]);
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "#build:$Id$\n");

	vector < Triangle > mesh;
	vector < Point >  points;
	vector < int >  boundary;

	switch (type) {
	case 6:
	case 1:
		build_hemisphere(mesh, points);
		break;
	case 0:
		build_icosahedron(mesh, points);
		break;
	case 2:
		build_test(mesh, points);
		break;
	case 3:
	case 4:
		build_icosahedron(mesh, points);
		break;
	case 5:
		build_test5(mesh, points);
		break;
	default:
		assert(0);
	}

	surface_projector project = sphere_orto_projector;
	if (type == 2 || type == 5 || type == 6) project = sphere_z_projector;
	filter_condition ok = 0;
	if (type < 3 || type == 5) {
		ok = 0; // alwais
	} else if (type == 3) {
		ok = ok1;
	} else if (type == 6) {
		ok = ok3;
	} else {
		ok = ok2;
	}

	normalize_mesh(mesh, points, project);
	iterate_mesh(mesh, points, iters, project);
	filter_mesh(mesh, points, boundary, ok);
	print_mesh(mesh, points, boundary, type, local);
	return 0;
}

