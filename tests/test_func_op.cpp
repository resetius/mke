#undef NDEBUG

#define _USE_MATH_DEFINES

#include <assert.h>
#include <math.h>
#include "func.h"
#include "util.h"
#include "mesh_builder.h"

using namespace std;
using phelm::FuncPtr;
using phelm::Symb;
using phelm::Cos;
using phelm::Sin;
using phelm::Pow;
using phelm::ATan2;
using phelm::ATan3;
using phelm::vars_t;

static void check(double a, double b)
{
	fprintf(stderr, "checking %le ?= %le\n", a, b);
	if (fabs(a - b) > 1e-14)
	{
		fprintf(stderr, "%.16lf != %.16lf\n", a, b);
		exit(-1);
	}
}

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
	// TODO: check ans
}

static void test_op() {
	vector<Triangle> tr;
	vector<Point> ps;
	build_icosahedron(tr, ps);

	for (size_t i = 0; i < tr.size(); i++) {
		test_tr(tr[i], ps);
		break;
	}
}

static void check_cos(FuncPtr cos_)
{
	cos_->bind_args({ "x" });
	check(cos_->apply({ 0.0 })->value(), 1);
	check(cos_->apply({ M_PI })->value(), -1);
	check(cos_->apply({ M_PI_2 })->value(), 0);
	check(cos_->apply({ -M_PI_2 })->value(), 0);
	check(cos_->apply({ M_PI / 4 })->value(), 0.7071067811865476);
	check(cos_->apply({ -M_PI / 4 })->value(), 0.7071067811865476);
	check(cos_->apply({ M_PI / 3 })->value(), 0.5);
	check(cos_->apply({ -M_PI / 3 })->value(), 0.5);
}

static void check_sin(FuncPtr sin_)
{
	sin_->bind_args({ "x" });
	check(sin_->apply({ 0.0 })->value(), 0);
	check(sin_->apply({ M_PI })->value(), 0);
	check(sin_->apply({ M_PI_2 })->value(), 1);
	check(sin_->apply({ -M_PI_2 })->value(), -1);
	check(sin_->apply({ M_PI / 4 })->value(), 0.7071067811865476);
	check(sin_->apply({ -M_PI / 4 })->value(), -0.7071067811865476);
	check(sin_->apply({ M_PI / 3 })->value(), 0.8660254037844386);
	check(sin_->apply({ -M_PI / 3 })->value(), -0.8660254037844386);
}

static void test_func() {
	FuncPtr X(new Symb("x"));
	FuncPtr Y(new Symb("y"));
	FuncPtr cos_(new Cos(X));
	check_cos(cos_);
	
	FuncPtr sin_(new Sin(X));
	check_sin(sin_);

	check_cos(sin_->diff("x"));
	check_sin(cos_->diff("x")*-1.0);

	X->bind_args({ "x" });
	check(X->apply({ 1.0 })->value(), 1.0);
	check(X->apply({ M_PI })->value(), M_PI);

	check(X->diff("y")->value(), 0.0);
	check(X->diff("x")->value(), 1.0);

	FuncPtr atan2_(new ATan2(Y, X));
	FuncPtr atan3_(new ATan3(Y, X));
	atan2_->bind_args({ "y", "x" });
	atan3_->bind_args({ "y", "x" });
	check(atan2_->apply({ 1, 0 })->value(), atan2(1, 0));
	check(atan2_->apply({ 1, 1 })->value(), atan2(1, 1));
	check(atan2_->apply({ 0, 1 })->value(), atan2(0, 1));
	check(atan2_->apply({ -1, 1 })->value(), atan2(-1, 1));
	check(atan2_->apply({ -1, 0 })->value(), atan2(-1, 0));
	check(atan2_->apply({ -1, -1 })->value(), atan2(-1, -1));
	check(atan2_->apply({ 0, -1 })->value(), atan2(0, -1));
	check(atan2_->apply({ 1, -1 })->value(), atan2(1, -1));

	check(atan3_->apply({ 1, 0 })->value(), atan3(1, 0));
	check(atan3_->apply({ 1, 1 })->value(), atan3(1, 1));
	check(atan3_->apply({ 0, 1 })->value(), atan3(0, 1));
	check(atan3_->apply({ -1, 1 })->value(), atan3(-1, 1));
	check(atan3_->apply({ -1, 0 })->value(), atan3(-1, 0));
	check(atan3_->apply({ -1, -1 })->value(), atan3(-1, -1));
	check(atan3_->apply({ 0, -1 })->value(), atan3(0, -1));
	check(atan3_->apply({ 1, -1 })->value(), atan3(1, -1));

	FuncPtr sqrt_(new Pow(X, 0.5));
	sqrt_->bind_args({ "x" });
	check(sqrt_->apply({ 0 })->value(), 0);
	check(sqrt_->apply({ 2 })->value(), sqrt(2));
	check(sqrt_->apply({ M_PI })->value(), sqrt(M_PI));

	FuncPtr dsqrt = sqrt_->diff("x");
	dsqrt->bind_args({ "x" });
	check(dsqrt->apply({ 2 })->value(), 0.5 / sqrt(2));
	check(dsqrt->apply({ M_PI })->value(), 0.5 / sqrt(M_PI));

	FuncPtr dsqrt1 = dsqrt->diff("x");
	dsqrt1->bind_args({ "x" });
	check(dsqrt1->apply({ 2 })->value(), 0.5*-0.5 * pow(2, -0.5-1));
}

extern "C" int test_func_op(int argc, char ** argv)
{
	test_op();
	test_func();
	return 0;
}
