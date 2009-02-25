#ifndef POLYNOM_H
#define POLYNOM_H
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (������� ���������)
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

struct Polynom {
	int deg_;                     /* !<degree    */
	int n_;                       /* !<dimension */

	int size_;                    /* !<����� ����-��� */
	std::vector < double > koef_; /* !<����� ����-��� */

/**
 * ������� ����� �������
 * @param deg - �������
 * @param n - ����������� (����� ����������� ����������, �������� ������� �� x, ������� �� x,y, ������� �� x,y,z)
 * @return �������
 */
	Polynom(int deg, int n);
	Polynom(int deg, int n, double * koef, int l);
	~Polynom();

	void print() const;

	double apply(const double * x) const;
	double apply(double x, double y) const;

	void operator /= (double k);
};

extern const Polynom P2X; //p(x, y) = x
extern const Polynom P2Y; //p(x, y) = y

struct Triangle;
struct Point;

//!���������� ����������� �������� p �� ���������� i
Polynom diff(const Polynom & p, int i);
//!���������� �������� �������� p �� ������������ t
double integrate(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!int p cos x dx dy
double integrate_cos(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!int p sin x dx dy
double integrate_sin(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!int p / cos(x) dx dy
double integrate_1_cos(const Polynom & p, const Triangle & t, const std::vector < Point > & ps);
//!������������ ���������
Polynom operator * (const Polynom &p1, const Polynom &p2);
//!�������� ���������
Polynom operator - (const Polynom &p1, const Polynom &p2);
//!����� ���������
Polynom operator + (const Polynom &p1, const Polynom &p2);

Polynom operator - (const Polynom &p1, double k);
Polynom operator * (const Polynom &p1, double k);

#endif /* POLYNOM_H */