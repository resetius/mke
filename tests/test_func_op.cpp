#undef NDEBUG

#include <assert.h>
#include "func.h"
#include "mesh_builder.h"

using namespace std;
using phelm::FuncPtr;
using phelm::Symb;
using phelm::Cos;
using phelm::Sin;
using phelm::vars_t;

static void test_tr(const Triangle &tr, vector<Point> & ps)
{
	FuncPtr phi(new Symb("phi"));
	FuncPtr la(new Symb("la"));

	FuncPtr X(new Symb("x"));
	FuncPtr Y(new Symb("y"));
	FuncPtr p0 =
		(X - ps[tr.v2].x) * (ps[tr.v3].y - ps[tr.v2].y) -
		(Y - ps[tr.v2].y) * (ps[tr.v3].x - ps[tr.v2].x);

	FuncPtr p1 =
		(X - ps[tr.v1].x) * (ps[tr.v3].y - ps[tr.v1].y) -
		(Y - ps[tr.v1].y) * (ps[tr.v3].x - ps[tr.v1].x);

	FuncPtr p2 =
		(X - ps[tr.v1].x) * (ps[tr.v2].y - ps[tr.v1].y) -
		(Y - ps[tr.v1].y) * (ps[tr.v2].x - ps[tr.v1].x);

	p0->bind_args({ "x", "y" });
	p1->bind_args({ "x", "y" });
	p2->bind_args({ "x", "y" });
	p0 = p0 * p0->apply({ ps[tr.v1].x, ps[tr.v1].y });
	p1 = p1 * p1->apply({ ps[tr.v2].x, ps[tr.v2].y });
	p2 = p2 * p2->apply({ ps[tr.v3].x, ps[tr.v3].y });

	p0->print();
	printf("\n\n");
	FuncPtr x1 = FuncPtr(new Cos(phi))*FuncPtr(new Cos(la));
	FuncPtr y1 = FuncPtr(new Cos(phi))*FuncPtr(new Sin(la));

	vars_t compose;
	compose["x"] = x1;
	compose["y"] = y1;

	FuncPtr p01 = p0->apply(compose);
	p01->print();
	printf("\n\n");

	p01->diff("la")->print();
	printf("\n\n");
}

extern "C" int test_func_op(int argc, char ** argv)
{
	vector<Triangle> tr;
	vector<Point> ps;
	build_icosahedron(tr, ps);

	for (size_t i = 0; i < tr.size(); i++) {
		test_tr(tr[i], ps);
		break;
	}
	return 0;
}
