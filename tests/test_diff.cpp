
#include <stdlib.h>
#include <math.h>
#include "mesh.h"

using namespace phelm;
using namespace std;

static bool cmp (double a, double b)
{
	if (fabs (a - b) > 1e-14)
	{
		return false;
	}
	return true;
}

static bool do_all (Mesh & m)
{
	int sz = m.size;
	double v;
	bool r = true;
	for (int i = 0; i < sz; ++i)
	{
		vector < double > vals;
		int zone = m.tr[m.adj[i][0]].z;
		for (uint tk = 0; tk < m.adj[i].size(); ++tk)
		{
			int trk_i = m.adj[i][tk];
			const Triangle & trk = m.tr[trk_i];
			const Polynom & p = trk.elem1 (i, zone);
			Polynom d = diff (p, 1);
			int z = trk.z;
			v = d.apply (m.ps[i].x (z), m.ps[i].y (z) );
			vals.push_back (v);
		}

		for (int i = 0; i < (int) vals.size(); ++i)
		{
			if (!cmp (v, vals[i]) )
			{
				fprintf (stderr, "%.16lf != %.16lf\n", v, vals[i]);
				r = false;
			}
		}
	}
	return r;
}

extern "C" int test_diff (int argc, char ** argv)
{
	Mesh m;
	bool r;
	FILE * f = fopen (argv[1], "r");
	m.load (f);
	fclose (f);
	r = do_all (m);
	fprintf (stderr, "result = %d\n", (int) r);
	return !r;
}

