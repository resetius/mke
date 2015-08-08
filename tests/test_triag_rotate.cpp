#undef NDEBUG

#include <assert.h>
#include "mesh_builder.h"
#include "point.h"

using namespace std;
using phelm::Matrix;

static double sign(double a) {
	if (a < 0) {
		return -1;
	}
	else {
		return 1;
	}
}

static void rotate_tr(const Triangle &tr, vector<Point> & ps)
{
	// единичная нормаль
	Point n = (ps[tr.v1] - ps[tr.v2]) * (ps[tr.v1] - ps[tr.v3]);
	n = n / n.len();

	// углы в плоскости (x, y)
	double cosa = fabs(n.x) / sqrt(n.x*n.x + n.y*n.y);
	// углы в плоскости (y, z) (после первой ротации будет эта плоскость)
	double cosb = fabs(n.z);
	// R1
	// cosa -sina 0
	// sina  cosa 0
	// 0     0    1
	// R2
	// cosb 0 -sinb
	// 0    1  0
	// sinb 0  cosb
	
	Matrix m;
	Point rn = n;
	double a = -sign(rn.x*rn.y)*acos(cosa);
	rn = rn.rotate_z(a); // -> move to (x,z)
	double s = scalar(rn, Point(0, 1, 0));
	assert(s < 1e-10);

	double b =  sign(rn.x*rn.z)*acos(cosb);
	rn = rn.rotate_y(b);

	assert(fabs(rn.x) < 1e-10);
	assert(fabs(rn.y) < 1e-10);
	assert(fabs(rn.z*rn.z - 1.0) < 1e-10);

	m.rotate_z(a);
	m.rotate_y(b);

	Point rn2 = n.apply(m);

	fprintf(stdout, ">\n");
	n.print();
	rn.print();
	rn2.print();

	assert((rn2 - rn).len() < 1e-15);

	vector<Point> trps(3);
	vector<Point> rtrps(3);
	trps[0] = ps[tr.v1];
	trps[1] = ps[tr.v2];
	trps[2] = ps[tr.v3];

	double aa = (trps[1] - trps[0]).len();
	double ab = (trps[1] - trps[2]).len();
	double ac = (trps[0] - trps[2]).len();
	double p1 = aa + ab + ac;
	double S1 = sqrt(p1 *(p1 - aa)*(p1 - ab)*(p1 - ac));

	rtrps[0] = ps[tr.v1].rotate_z(a).rotate_y(b);
	rtrps[1] = ps[tr.v2].rotate_z(a).rotate_y(b);
	rtrps[2] = ps[tr.v3].rotate_z(a).rotate_y(b);

	double ba = (rtrps[1] - rtrps[0]).len();
	double bb = (rtrps[1] - rtrps[2]).len();
	double bc = (rtrps[0] - rtrps[2]).len();
	double p2 = ba + bb + bc;
	double S2 = sqrt(p2 *(p2 - ba)*(p2 - bb)*(p2 - bc));

	assert(fabs(S2 - S1) < 1e-10);
	assert(fabs(rtrps[0].z - rtrps[1].z) < 1e-10);
	assert(fabs(rtrps[0].z - rtrps[2].z) < 1e-10);

	fprintf(stdout, "{\n");
	trps[0].print();
	trps[1].print();
	trps[2].print();
	fprintf(stdout, "} => \n");
	fprintf(stdout, "{\n");
	rtrps[0].print();
	rtrps[1].print();
	rtrps[2].print();
	fprintf(stdout, "} \n\n");
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
