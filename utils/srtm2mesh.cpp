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
#include <stdio.h>
#include <math.h>

#include <vector>
#include <string>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#ifdef WIN32
#include <windows.h>
#define htons(x) \
	(((((unsigned short)(x)) >>8) & 0xff) | \
                ((((unsigned short)(x)) & 0xff)<<8))
#else
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <arpa/inet.h>
#include <stdint.h>
#endif

using namespace std;

#ifdef WIN32
	typedef __int16 int16_t;
	typedef HANDLE file_t;
#else
	typedef int file_t;
#endif

struct Data {
	int16_t * d_;

#ifdef WIN32
	HANDLE fd_;
	HANDLE map_;
#else
	int fd_;
	size_t sz_;
#endif

	Data(): d_(0), 
		fd_((file_t)-1)
	{
	}

	~Data()
	{
#ifdef WIN32
		if (fd_ >= (file_t)-1) CloseHandle(fd_);
		if (map_) CloseHandle(map_);
		if (d_) UnmapViewOfFile(d_);
		map_ = 0;
#else
		if (d_) munmap(d_, sz_);
		if (fd_ >= 0) close(fd_);
#endif
		d_  = 0;
		fd_ = (file_t)-1;
	}

	bool load(const char * path) {
#ifdef WIN32
		fd_ = CreateFile(path, GENERIC_READ, NULL, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

		if (fd_ < 0) {
			return false;
		}

		map_ = CreateFileMapping(fd_, NULL, PAGE_READONLY, NULL, NULL, path);
		if (!map_) {
			CloseHandle(fd_);
			return false;
		}

		d_ = (int16_t *)MapViewOfFile(map_, FILE_MAP_READ, NULL, NULL, NULL);
		if (!d_) {
			CloseHandle(fd_);
			CloseHandle(map_);
			return false;
		}

		return true;
#else
		struct stat buf;
		fd_ = open(path, O_RDONLY);
		if (fd_ < 0) {
			return false;
		}

		if (stat(path, &buf) == 0) {
			sz_ = buf.st_size;
		} else {
			close(fd_); fd_ = -1;
			return false;
		}

		assert(sz_ >= 2 * 43200 * 21600);

		d_ = (int16_t*)mmap(0, sz_, PROT_READ, MAP_SHARED, fd_, 0);
		if (d_ == (void*)-1) {
			close(fd_); fd_ = -1;
			return false;
		}

		return true;
#endif
	}

	double get(double u, double v) {
		int i, j;
		assert(-M_PI / 2 <= u && u <= M_PI / 2);
		assert(0 <= v && v<= 2 * M_PI);
		j = (int)((v * (43200.0 - 1.0)) / 2.0 / M_PI);
		i = (int)((u + M_PI / 2) * (21600.0 - 1.0) / M_PI);
		i = 21599 - i;

		return (double)((int16_t)htons(d_[i * 43200 + j]));
	}
};

struct Point {
	double x;
	double y;
	double z;

	Point(): x(0), y(0), z(0) {}
	Point(double x1, double y1, double z1): x(x1), y(y1), z(z1) {}
};

/**
 * format:
 * # comments
 * # comments
 * x y z
 * x y z
 * .....
 * # not necessary triangles data 
 * n1 n2 n3
 * ....
 * # not necessary boundary data
 * n1
 * n1
 * ....
 */
struct SimpleMesh {
	vector < Point > ps;
	vector < string > triangles;

	SimpleMesh() {}

	bool load(FILE * f) {
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

			if (sscanf (s, "%lf%lf%lf", &x, &y, &z) != 3) {
				return false;
			}

			ps.push_back (Point(x, y, z)); lineno++;
		}
		while (fgets (s, _BUF_SZ - 1, f));

		// not necessary triangles data 
		// for vizualizer 3d compatibility

		triangles.push_back(s); lineno++;

		if (!fgets (s, _BUF_SZ - 1, f)) { 
			return true;
		}

		do
		{
			if (*s == '#')
				break;
			triangles.push_back(s);
			lineno++;
		}
		while (fgets (s, _BUF_SZ - 1, f));

		return true;
	}
};

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

void make_relief(vector < double > & relief, SimpleMesh & mesh, Data & data)
{
	size_t sz = mesh.ps.size();
	relief.resize(sz);

	for (size_t i = 0; i < sz; ++i) {
		Point & p = mesh.ps[i];
		double u1 = u(p.x, p.y, p.z); // latitude
		double v1 = v(p.x, p.y, p.z); // longitude

		relief[i] = data.get(u1, v1);
	}
}

void make_relief(vector < double > & relief, int n_phi, int n_la, Data & data)
{
	relief.resize(n_phi * n_la);

	double d_la  = 2.0*M_PI/(n_la); // dx on new grid
	double d_phi = M_PI/(2.0*(n_phi-1.0)+1.0); // dy on new grid
	double omg   = 2.*M_PI/24./60./60.; // ?
	double T     = 1./omg;

	for (int i = 0; i < n_phi; ++i) 
	{
        for (int j = 0; j < n_la; ++j) 
		{
			double u1 = i * d_phi; // latitude
			double v1 = j * d_la; // longitude

			relief[i * n_la + j] = T * data.get(u1, v1);
		}
	}
}

void print_relief(vector < double > & relief, int n_phi, int n_la, int type)
{
	for (int i = 0; i < n_phi; ++i) 
	{
        for (int j = 0; j < n_la; ++j) 
		{
			fprintf(stdout, "%.16le ", relief[i * n_la + j]);
		}
		fprintf(stdout, "\n");
	}
}

void print_relief(vector < double > & relief, SimpleMesh & mesh, int type)
{
	size_t sz = mesh.ps.size();

	if (type == 0) {
		for (size_t i = 0; i < sz; ++i) {
			fprintf(stdout, "%lf ", relief[i]);
		}
	} else {
		fprintf(stdout, "#\n");

		for (size_t i = 0; i < sz; ++i) {
			Point & p = mesh.ps[i];
			fprintf(stdout, "%.16lf %.16lf %.16lf %.16lf\n", p.x, p.y, p.z, relief[i]);
		}

		sz = mesh.triangles.size();
		for (size_t i = 0; i < sz; ++i) {
			fprintf(stdout, "%s", mesh.triangles[i].c_str());
		}
		fprintf(stdout, "#");
	}
}

void usage(const char * name)
{
	fprintf(stderr, "usage: %s --srtm [srtmpath] --type [f|tri] --mesh [mesh.txt]\n", name);
	fprintf(stderr, "--srtm [srtmpath] -- path to srtm directory (. is default)\n"
		"--type [f|tri] -- type of output\n"
		"    f    - only f data (elevation)\n"
		"    tri  - full data (for vizualizer 3d)\n"
		"--mesh [mesh.txt] - mesh generated by 'sphere', use global coordinates !!\n");
	exit(1);
}

int main(int argc, char * argv[])
{
	vector < double > relief;
	const char * srtm_path = ".";
	const char * mesh_file = 0;
	int output_type   = 0; //0 - f only, 1 - full
	SimpleMesh mesh;
	Data data;

	FILE * f;

	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "--h"))
		{
			usage(argv[0]);
		} else if (!strcmp(argv[i], "--type")) {
			if (i == argc - 1) {
				fprintf(stderr, "--type requires an argument\n");
				usage(argv[0]);
			}

			if (!strcmp(argv[i + 1], "f")) {
				output_type = 0;
			} else if (!strcmp(argv[i + 1], "tri")) {
				output_type = 1;
			}
		} else if (!strcmp(argv[i], "--mesh")) {
			if (i == argc - 1) {
				fprintf(stderr, "--mesh requires an argument\n");
				usage(argv[0]);
			}

			mesh_file = argv[i + 1];
		} else if (!strcmp(argv[i], "--srtm")) {
			if (i == argc - 1) {
				fprintf(stderr, "--srtm requires an argument\n");
				usage(argv[0]);
			}

			srtm_path = argv[i + 1];
		}
	}

	if (!data.load(srtm_path)) {
		fprintf(stderr, "cannot load data %s\n", srtm_path);
		usage(argv[0]);
	}

	if (mesh_file) {
		f = fopen(mesh_file, "rb");
		if (!f) {
			fprintf(stderr, "cannot open %s\n", mesh_file);
			usage(argv[0]);
		}

		if (!mesh.load(f)) {
			fprintf(stderr, "cannot load %s\n", mesh_file);
			usage(argv[0]);
		}

		make_relief(relief, mesh, data);
		print_relief(relief, mesh, output_type);
	} else {
		make_relief(relief, 24, 32, data);
		print_relief(relief, 24, 32, output_type);
	}

	return 0;
}

