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
 
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

template < typename A >
class Point
{
public:
	A x;
	A y;
	A z;
	A u;
	A v;

	template < typename T >
	Point(T x1, T y1, T z1, T u1, T v1): 
		x(x1), y(y1), z(z1), u(u1), v(v1)
	{
	}

	Point() : x(0), y(0), z(0), u(0), v(0)
	{
	}
};

typedef Point < double > Pointd;

void load_file(std::vector < Pointd > & points,  //!<координаты точек
			   std::vector < double > & colors,  //!<номера узлов
			   std::vector < std::vector < int > > & tri, 
			   FILE * f)
{
	int size;
	int lineno = 0;

#define _BUF_SZ 32768
	char s[_BUF_SZ];

	if (!fgets (s, _BUF_SZ - 1, f))
		goto bad;

	lineno ++;
	do
	{
		if (*s != '#')
			break;

		lineno ++;
	}
	while (fgets (s, _BUF_SZ - 1, f) );

	do
	{
		double x = 0.0, y = 0.0, z = 0.0, u = 0.0, v = 0.0;
		double f = 0.0;
		int m;

		if (*s == '#')
			break;

		m = sscanf (s, "%lf%lf%lf%lf%lf%lf", &x, &y, &z, &f, &u, &v);
		if (m < 2)
		{
			goto bad;
		}

		points.push_back(Pointd (x, y, z, u, v));
		colors.push_back(f);
		lineno ++;
	}
	while (fgets (s, _BUF_SZ - 1, f) );

	size = (int) points.size();

	fgets (s, _BUF_SZ - 1, f); lineno ++;
	do
	{
		int n1, n2, n3;
		int zone = 1;
		const char * sep = ";";
		char * str = s;


		if (*s == '#') {
			break;
		}

		str = strtok(s, sep);
		if (sscanf (str, "%d%d%d", &n1, &n2, &n3) != 3)
		{
			goto bad;
		}
		str = strtok(0, sep);
		if (str) {
			sscanf(str, "%d", &zone);
		}

		//так как индексы в файле с 1 а не с 0
		--n1;
		--n2;
		--n3;
		--zone;
		if (n1 >= size || n2 >= size || n3 >= size)
		{
			goto bad;
		}

		vector < int > tr(4);
		tr.push_back(n1); tr.push_back(n2); tr.push_back(n3); tr.push_back(zone);

		tri.push_back(tr);

		lineno ++;
	}
	while (fgets (s, _BUF_SZ - 1, f));

	return;

bad:
	{
		fclose(f);
		ostringstream ss;
		ss << "bad file format, lineno: " << lineno;
		string out = ss.str();
		throw logic_error (out.c_str());
	}
}

void print_func(const char * fn,
				std::vector < Pointd > & points,  
			    std::vector < double > & colors,  
			    std::vector < std::vector < int > > & tri)
{
	FILE * f = fopen(fn, "w");
	if (!f) return;

	fprintf(f, "#points %lu\n", points.size()); 
	for (size_t i = 0; i < points.size(); ++i)
	{
		Pointd & p = points[i];
		double & F = colors[i];
		fprintf(f, "%.16lf %16lf %.16lf %16lf %.16lf %16lf \n", 
			p.x, p.y, p.z, F, p.u, p.v);
	}

	fprintf(f, "# triangles %lu\n", tri.size());
	for (size_t i = 0; i < tri.size(); ++i)
	{
		fprintf(f, "%d %d %d %d\n", tri[i][0] + 1, tri[i][1] + 1, 
			tri[i][2] + 1,  tri[i][3] + 1);
	}
	fprintf(f, "# end \n");
	fclose(f);
}

void analize(FILE * fp)
{
	std::vector < Pointd >  points;  //!<координаты точек
	std::vector < double >  f;       //!<номера узлов
	std::vector < double >  m;
	std::vector < double >  m2;
	std::vector < double >  d;
	std::vector < std::vector < int > > tri; 
	int n = 0;

	while (!feof(fp)) {
		points.clear(); f.clear(); tri.clear();
		try {
			load_file(points, f, tri, fp);
		} catch (exception & e) {
			fprintf(stderr, "%s\n", e.what());
			break;
		}

		if (m.empty()) {
			m.resize(f.size());
			m2.resize(f.size());
		}

		double k1 = (double)n / (double)(n + 1);
		double k2 = 1.0 / (double)(n + 1);

		for (int i = 0; i < (int)m.size(); ++i) {
			m[i]  = k1 * m[i]  + k2 * f[i];
			m2[i] = k1 * m2[i] + k2 * f[i] * f[i];
		}
	}

	print_func("average.txt", points, m, tri);
	print_func("dispersion.txt", points, m2, tri);
}

int main(int argc, char * argv[])
{
	FILE * f = stdin;

	if (argc > 1) {
		f = fopen(argv[1], "r");
		if (!f) {
			fprintf(stderr, "cannot open %s\n", argv[1]);
			return -1;
		}
	}

	analize(f);

	fclose(f);
}
