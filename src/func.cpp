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

#include "func.h"

using namespace phelm;

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

	if (newa->has_value() && b->has_value()) {
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
		return FuncPtr(new Const(cos(newa->value())));
	}
	else {
		return FuncPtr(new Sin(newa));
	}
}

FuncPtr Cos::diff(const std::string & symb) const {
	if (a->has_symb(symb)) {
		FuncPtr a(new Sin(a));
		return a * -1.0;
	}
	else {
		return FuncPtr(new Const(0));
	}
}

FuncPtr Sin::diff(const std::string & symb) const {
	if (a->has_symb(symb)) {
		FuncPtr a(new Cos(a));
		return a;
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
		double aa = fabs(a->value());
		if (aa < 1e-15) {
			return FuncPtr(new Const(0.0));
		}
		if (fabs(aa - 1.0) < 1e-15) {
			return b;
		}
	}
	if (b->has_value()) {
		double bb = fabs(b->value());
		if (bb < 1e-15) {
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
