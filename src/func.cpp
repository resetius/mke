/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2015 Alexey Ozeritsky (Алексей Озерицкий)
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
#include <assert.h>
#include "func.h"
#include "util.h"

using namespace phelm;

FuncPtr Func::apply(std::initializer_list<double> a) {
	std::vector<double> values;
	values.insert(values.end(), a.begin(), a.end());
	vars_t vars;
	assert(args.size() >= values.size());
	for (int i = 0; i < (int)args.size() && i < (int)values.size(); i++) {
		vars[args[i]] = FuncPtr(new Const(values[i]));
	}
	return apply(vars);
}

FuncPtr Func::apply(std::initializer_list<FuncPtr> a) {
	std::vector<FuncPtr> values;
	values.insert(values.end(), a.begin(), a.end());
	vars_t vars;
	assert(args.size() >= values.size());
	for (int i = 0; i < (int)args.size() && i < (int)values.size(); i++) {
		vars[args[i]] = values[i];
	}
	return apply(vars);
}

FuncPtr Symb::apply(const vars_t & vars) {
	vars_t::const_iterator it = vars.find(name);
	if (it == vars.end()) {
		return shared_from_this();
	}
	else {
		return it->second;
	}
}

FuncPtr BinOp::apply(const vars_t & vars) {
	FuncPtr newa = a->apply(vars);
	FuncPtr newb = b->apply(vars);

	if (newa->has_value() || newb->has_value()) {
		return op(newa, newb);
	}
	else {
		return FuncPtr(new BinOp(newa, newb, op, op_symb, lin_diff));
	}
}

FuncPtr BinOp::diff(const std::string & symb) const {
	if (has_symb(symb)) {
		if (lin_diff) {
			return op(a->diff(symb), b->diff(symb));
		}
		else {
			return op(a->diff(symb), b) + op(a, b->diff(symb));
		}
	}
	else {
		return FuncPtr(new Const(0.0));
	}
}

FuncPtr Cos::apply(const vars_t & vars) {
	FuncPtr newa = a->apply(vars);
	if (newa->has_value()) {
		return FuncPtr(new Const(cos(newa->value())));
	}
	else {
		return FuncPtr(new Cos(newa));
	}
}

FuncPtr Sin::apply(const vars_t & vars) {
	FuncPtr newa = a->apply(vars);
	if (newa->has_value()) {
		return FuncPtr(new Const(sin(newa->value())));
	}
	else {
		return FuncPtr(new Sin(newa));
	}
}

FuncPtr Cos::diff(const std::string & symb) const {
	if (a->has_symb(symb)) {
		FuncPtr t(new Sin(a));
		return t * -1.0;
	}
	else {
		return FuncPtr(new Const(0));
	}
}

FuncPtr Sin::diff(const std::string & symb) const {
	if (a->has_symb(symb)) {
		FuncPtr t(new Cos(a));
		return t;
	}
	else {
		return FuncPtr(new Const(0));
	}
}

FuncPtr phelm::operator * (const FuncPtr & a, const FuncPtr & b) {
	if (a->has_value() && b->has_value()) {
		return FuncPtr(new Const(a->value() * b->value()));
	}
	if (a->has_value()) {
		double aa = a->value();
		if (fabs(aa) < 1e-15) {
			return FuncPtr(new Const(0.0));
		}
		if (fabs(aa - 1.0) < 1e-15) {
			return b;
		}
	}
	if (b->has_value()) {
		double bb = b->value();
		if (fabs(bb) < 1e-15) {
			return FuncPtr(new Const(0.0));
		}
		if (fabs(bb - 1.0) < 1e-15) {
			return a;
		}
	}
	return FuncPtr(new Mul(a, b));
}

FuncPtr phelm::operator + (const FuncPtr & a, const FuncPtr & b) {
	if (a->has_value() && b->has_value()) {
		return FuncPtr(new Const(a->value() + b->value()));
	}
	if (a->has_value()) {
		double aa = fabs(a->value());
		if (aa < 1e-15) {
			return b;
		}
	}
	if (b->has_value()) {
		double bb = fabs(b->value());
		if (bb < 1e-15) {
			return a;
		}
	}
	return FuncPtr(new Add(a, b));
}

FuncPtr phelm::operator - (const FuncPtr & a, const FuncPtr & b) {
	if (a->has_value() && b->has_value()) {
		return FuncPtr(new Const(a->value() - b->value()));
	}
	if (b->has_value()) {
		double bb = fabs(b->value());
		if (bb < 1e-15) {
			return a;
		}
	}
	return FuncPtr(new Sub(a, b));
}

FuncPtr phelm::operator ^ (const FuncPtr & a, double n) {
	if (fabs(n) < 1e-15) {
		return FuncPtr(new Const(1));
	}

	if (fabs(n - 1) < 1e-15) {
		return a;
	}

	if (a->has_value()) {
		double v = pow(a->value(), n);
		return FuncPtr(new Const(v));
	}

	return FuncPtr(new Pow(a, n));
}

FuncPtr Pow::apply(const vars_t & vars) {
	FuncPtr b = a->apply(vars);
	if (b->has_value()) {
		return b ^ n;
	}
	else {
		return FuncPtr(new Pow(b, n));
	}
}

FuncPtr Pow::diff(const std::string & symb) const {	
	if (!has_symb(symb)) {
		return FuncPtr(new Const(0));
	}
	else {
		return (n * a->diff(symb)) * (a ^ (n - 1));
	}
}

FuncPtr ASin::apply(const vars_t & vars) {
	FuncPtr newa = a->apply(vars);
	if (newa->has_value()) {
		return FuncPtr(new Const(asin(newa->value())));
	}
	else {
		return FuncPtr(new ASin(newa));
	}
}

FuncPtr ASin::diff(const std::string & symb) const {
	throw std::runtime_error("not implemented");
}

FuncPtr ATan2::apply(const vars_t & vars) {
	FuncPtr newx = x->apply(vars);
	FuncPtr newy = y->apply(vars);
	if (newx->has_value() && newy->has_value()) {
		double v = atan2(newx->value(), newy->value());
		return FuncPtr(new Const(v));
	}
	else {
		return FuncPtr(new ATan2(newx, newy));
	}
}

FuncPtr ATan2::diff(const std::string & symb) const {
	throw std::runtime_error("not implemented");
}

FuncPtr ATan3::apply(const vars_t & vars) {
	FuncPtr newx = x->apply(vars);
	FuncPtr newy = y->apply(vars);
	if (newx->has_value() && newy->has_value()) {
		double v = atan3(newx->value(), newy->value());
		return FuncPtr(new Const(v));
	}
	else {
		return FuncPtr(new ATan3(newx, newy));
	}
}

FuncPtr ATan3::diff(const std::string & symb) const {
	throw std::runtime_error("not implemented");
}
