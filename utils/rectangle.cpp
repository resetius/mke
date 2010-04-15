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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <map>

#include "mesh_builder.h"

using namespace std;

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

void make_boundary(double x, double y, double w, double h, vector < Point > & points, vector < int > & boundary)
{
	for (size_t i = 0; i < points.size(); ++i) {
		if (fabs(points[i].x - x) < 1e-15 || fabs(points[i].y - y) < 1e-15
				|| fabs(points[i].x - x - w) < 1e-15 || fabs(points[i].y - y - h) < 1e-15)
		{
			boundary.push_back(i);
		}
	}
}


void print_mesh(FILE * f, double x, double y, double w, double h, 
		const vector < Triangle > & mesh, vector < Point > & points, vector < int > & boundary)
{
	for (vector < Point >::const_iterator it = points.begin();
		it != points.end(); ++it)
	{
		fprintf(f, "%.16lf %.16lf\n", it->x, it->y);
	}
	fprintf(f, "# triangles\n");

	for (vector < Triangle >::const_iterator it = mesh.begin();
		it != mesh.end(); ++it)
	{
		fprintf(f, "%d %d %d\n", it->v1 + 1, it->v2 + 1, it->v3 + 1);
	}
	fprintf(f, "# boundary\n");

	for (size_t i = 0; i < boundary.size(); ++i) {
		fprintf(f, "%lu \n", boundary[i] + 1);
	}
}

void usage(const char * name)
{
	fprintf(stderr, "usage: %s x y w h iters [-o output.txt] [-pgm output.pgm] [-nosort]\n", name);
	exit(1);
}

int main(int argc, char * argv[])
{
	if (argc < 6) {
		usage(argv[0]);
	}

	vector < Triangle > mesh;
	vector < Point    >  points;
	vector < int > boundary;

	double x  = atof(argv[1]);
	double y  = atof(argv[2]);
	double w  = atof(argv[3]);
	double h  = atof(argv[4]);
	int iters = atoi(argv[5]);
	bool sort = true;

	fprintf(stderr, "x, y, w, h = %lf, %lf, %lf, %lf\n", x, y, w, h);
	fprintf(stderr, "iterations = %d\n", iters);

	FILE * f1 = stdout;
	FILE * f2 = 0;

	for (int i = 6; i < argc; ++i) {
		if (!strcmp(argv[i], "-o")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			f1 = fopen(argv[i + 1], "wb");
			if (!f1) {
				fprintf(stderr, "cannot open %s\n", argv[6]);
				usage(argv[0]);
			}
		} else if (!strcmp(argv[i], "-pgm")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			f2 = fopen(argv[i + 1], "wb");
			if (!f2) {
				fprintf(stderr, "cannot open %s\n", argv[6]);
				usage(argv[0]);
			}
		} else if (!strcmp(argv[i], "-nosort")) {
			sort = false;
		}
	}

	build_rectangle(x, y, w, h, mesh, points);
	iterate_mesh(mesh, points, iters, flat_projector);
	make_boundary(x, y, w, h, points, boundary);
	if (sort) {
		reorder_mesh(mesh, points, boundary);
	}
	print_mesh(f1, x, y, w, h, mesh, points, boundary);

	fclose(f1);

	if (f2) {
		vizualize_adj(mesh, points, boundary, f2, 600, 600);
		fclose(f2);
	}
	return 0;
}

