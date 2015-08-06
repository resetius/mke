#undef NDEBUG

#include <assert.h>
#include "mesh_builder.h"

using namespace std;

static double sign(double a) {
	if (a < 0) {
		return -1;
	}
	else {
		return 1;
	}
}

void rotate_tr(const Triangle &tr, vector<Point> & ps)
{
	// единичная нормаль
	Point n = Point(
		ps[tr.v1].x + ps[tr.v2].x + ps[tr.v3].x,
		ps[tr.v1].y + ps[tr.v2].y + ps[tr.v3].y,
		ps[tr.v1].z + ps[tr.v2].z + ps[tr.v3].z
	);

	n = n / 3;
	n = n / n.abs();

	// углы в плоскости (x, y)
	double cosa = fabs(n.x) / sqrt(n.x*n.x + n.y*n.y);
	// углы в плоскости (y, z) (после первой ротации будет эта плоскость)
	double cosb = sqrt(n.x*n.x + n.y*n.y);
	// R1
	// cosa -sina 0
	// sina  cosa 0
	// 0     0    1
	// R2
	// cosb 0 -sinb
	// 0    1  0
	// sinb 0  cosb
	
	Point rn = n;
	double a = -sign(rn.x*rn.y)*acos(cosa);
	rn = rn.rotate_z(a); // -> move to (x,z)
	double s = rn.scalar(Point(0, 1, 0));
	assert(s < 1e-10);

	double b = -sign(rn.x*rn.z)*acos(cosb);
	rn = rn.rotate_y(b);

	assert(fabs(rn.x*rn.x - 1.0) < 1e-10);
	assert(fabs(rn.y) < 1e-10);
	assert(fabs(rn.z) < 1e-10);

	fprintf(stdout, ">\n");
	n.print();
	rn.print();
}

extern "C" int test_triag_rotate(int argc, char ** argv)
{
	//для простоты считаем, что треугольник лежит на сфере
	vector<Triangle> tr;
	vector<Point> ps;
	build_icosahedron(tr, ps);
	// foreach triangle test it

	for (size_t i = 0; i < tr.size(); i++) {
		rotate_tr(tr[i], ps);
	}

	return 0;
}
