#include "mesh_builder.h"

using namespace std;

void rotate_tr(const Triangle &tr, vector<Point> & ps)
{

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
