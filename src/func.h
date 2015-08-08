#ifndef PP_FUNC_H
#define PP_FUNC_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2015 Alexey Ozeritsky
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

#include <map>
#include <string>
#include <stdexcept>
#include <memory>
#include <functional>
#include <stdio.h>
#include <vector>
#include <initializer_list>

namespace phelm {

class Func;
typedef std::shared_ptr<Func> FuncPtr;
typedef std::map<std::string, FuncPtr> vars_t;

class Func : public std::enable_shared_from_this<Func> {
	std::vector<std::string> args;
public:
	virtual ~Func() {}
	virtual FuncPtr apply(const vars_t & vars) = 0;
	virtual bool has_value() const {
		return false;
	}
	virtual operator double() const {
		throw std::logic_error("cannot cast to double");
	}
	double value() const {
		return (double)*this;
	}

	void bind_args(std::initializer_list<std::string> a) {
		args.clear();
		args.insert(args.end(), a.begin(), a.end());
	}

	FuncPtr apply(std::initializer_list<double> a);

	virtual void print(FILE * f = stdout) const = 0;
	virtual bool has_symb(const std::string & symb) const = 0;
	virtual FuncPtr diff(const std::string & symb) const = 0;
};

class Const : public Func {
	double value;

public:
	Const(double value) : value(value) {}
	FuncPtr apply(const vars_t & vars) {
		return shared_from_this();
	}
	operator double() const {
		return value;
	}
	void print(FILE * f) const {
		fprintf(f, "%lf", value);
	}
	bool has_symb(const std::string & symb) const {
		return false;
	}
	bool has_value() const {
		return true;
	}
	FuncPtr diff(const std::string & symb) const {
		return FuncPtr(new Const(0));
	}
};

class Symb : public Func {
	std::string name;

public:
	Symb(const std::string & name) : name(name) {}
	FuncPtr apply(const vars_t & vars);

	void print(FILE * f) const {
		fprintf(f, "%s", name.c_str());
	}

	bool has_symb(const std::string & symb) const {
		return symb == name;
	}
	FuncPtr diff(const std::string & symb) const {
		if (symb == name) {
			return FuncPtr(new Const(1.0));
		} else {
			return FuncPtr(new Const(0.0));
		}
	}
};

FuncPtr operator - (const FuncPtr & a, const FuncPtr & b);
FuncPtr operator + (const FuncPtr & a, const FuncPtr & b);
FuncPtr operator * (const FuncPtr & a, const FuncPtr & b);

typedef std::function<FuncPtr(FuncPtr, FuncPtr)> Op;

class BinOp : public Func {
protected:
	FuncPtr a;
	FuncPtr b;
	Op op;
	std::string op_symb;
	bool lin_diff;

public:

	BinOp(FuncPtr a, FuncPtr b, Op op, const std::string op_symb, bool lin_diff = true) : 
		a(a), b(b), op(op), op_symb(op_symb), lin_diff(lin_diff) {}

	FuncPtr apply(const vars_t & vars);

	void print(FILE * f) const {
		fprintf(f, "(");
		a->print(f);
		fprintf(f, "%s", op_symb.c_str());
		b->print(f);
		fprintf(f, ")");
	}

	bool has_symb(const std::string & symb) const {
		return a->has_symb(symb) || b->has_symb(symb);
	}

	FuncPtr diff(const std::string & symb) const;
};

class Mul : public BinOp {
public:
	Mul(FuncPtr a, FuncPtr b) : BinOp(a, b, [](FuncPtr x, FuncPtr y){
		return x * y;
	}, "*", false) {}
};

class Add : public BinOp {
public:
	Add(FuncPtr a, FuncPtr b) : BinOp(a, b, [](FuncPtr x, FuncPtr y){
		return x + y;
	}, "+") {}
};

class Sub : public BinOp {
public:
	Sub(FuncPtr a, FuncPtr b) : BinOp(a, b, [](FuncPtr x, FuncPtr y){
		return x - y;
	}, "-") {}
};

class Cos : public Func {
	FuncPtr a;
public:
	Cos(FuncPtr a) : a(a) {}

	FuncPtr apply(const vars_t & vars);

	void print(FILE * f) const {
		fprintf(f, "cos(");
		a->print(f);
		fprintf(f, ")");
	}

	bool has_symb(const std::string & symb) const {
		return a->has_symb(symb);
	}

	FuncPtr diff(const std::string & symb) const;
};

class Sin : public Func {
	FuncPtr a;
public:
	Sin(FuncPtr a) : a(a) {}

	FuncPtr apply(const vars_t & vars);

	void print(FILE * f) const {
		fprintf(f, "sin(");
		a->print(f);
		fprintf(f, ")");
	}

	bool has_symb(const std::string & symb) const {
		return a->has_symb(symb);
	}

	FuncPtr diff(const std::string & symb) const;
};

inline FuncPtr operator * (const FuncPtr & a, double b) {
	return FuncPtr(new Mul(a, FuncPtr(new Const(b))));
}

inline FuncPtr operator * (double b, const FuncPtr & a) {
	return FuncPtr(new Mul(a, FuncPtr(new Const(b))));
}

inline FuncPtr operator + (const FuncPtr & a, double b) {
	return FuncPtr(new Add(a, FuncPtr(new Const(b))));
}

inline FuncPtr operator - (const FuncPtr & a, double b) {
	return FuncPtr(new Sub(a, FuncPtr(new Const(b))));
}

}

#endif /*PP_FUNC_H*/
