#ifndef POLYNOM_H
#define POLYNOM_H
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
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
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

#include <vector>
#include <string.h>

typedef unsigned int uint;

struct Polynom {
	int deg_;                     /* !<degree    */
	int n_;                       /* !<dimension */

	int size_;                    /* !<число коэф-тов */
	std::vector < double > koef_; /* !<масив коэф-тов */

/**
 * создает новый полином
 * @param deg - степень
 * @param n - размерность (число независимых переменных, например полином от x, полином от x,y, полином от x,y,z)
 * @return полином
 */
	Polynom(int deg, int n)
		: deg_(deg), n_(n), size_((deg + 1) * (deg + 1)), koef_(size_)
	{
	}

	Polynom(int deg, int n, double * koef, int l)
		: deg_(deg), n_(n), size_((deg + 1) * (deg + 1)), koef_(size_)
	{
		memcpy(&koef_[0], koef, l * sizeof(double));
	}

	~Polynom() {}

	void print() const;

	double apply(const double * x) const;
	double apply(double x, double y) const;

	void operator /= (double k) {
		uint sz   = (uint)koef_.size();
		double k1 = 1.0 / k;
		for (uint i = 0; i < sz; ++i) {
			koef_[i] *= k1;
		}
	}
};

extern const Polynom P2X; //p(x, y) = x
extern const Polynom P2Y; //p(x, y) = y

struct Triangle;
struct Point;

//!возвращает производную полинома p по переменной i
Polynom diff(const Polynom & p, int i);
//!возвращает интеграл полинома p по треугольнику t
double integrate(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!int p cos x dx dy
double integrate_cos(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!int p sin x dx dy
double integrate_sin(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!int p / cos(x) dx dy
double integrate_1_cos(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!произведение полиномов
Polynom operator * (const Polynom &p1, const Polynom &p2);
//!разность полиномов
Polynom operator - (const Polynom &p1, const Polynom &p2);
//!сумма полиномов
Polynom operator + (const Polynom &p1, const Polynom &p2);

inline Polynom operator - (const Polynom &p1, double x)
{
	Polynom r(p1);
	r.koef_[0] -= x;
	return r;
}

inline Polynom operator * (const Polynom &p1, double x)
{
	Polynom r(p1);
	uint sz = (uint)r.koef_.size();
	for (uint i = 0; i < sz; ++i)
	{
		r.koef_[i] *= x;
	}
	return r;
}

#endif /* POLYNOM_H */
