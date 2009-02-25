#ifndef MKE_H
#define MKE_H
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

#include <stdio.h>

#include "polynom.h"

typedef unsigned int uint;

struct Point {
	double x;
	double y;

	Point(): x(0), y(0) {}
	Point(double x1, double y1): x(x1), y(y1) {}
	Point(double *x1): x(x1[0]), y(x1[1]) {}

	Point operator / (double k)
	{
		return Point(x / k, y / k);
	}

	Point operator * (double k)
	{
		return Point(x * k, y * k);
	}
};

inline Point operator + (const Point & p1, const Point & p2)
{
	return Point(p1.x + p2.x, p1.y + p2.y);
}

struct Triangle {
	int p[3];

	Triangle(int p1_, int p2_, int p3_)
	{
		p[0] = p1_;
		p[1] = p2_;
		p[2] = p3_;
	}

	double x(int i, const std::vector < Point > & ps) const {
		return ps[p[i]].x;
	}

	double y(int i, const std::vector < Point > & ps) const {
		return ps[p[i]].y;
	}

	void print(const std::vector < Point > & ps) const {
		fprintf(stderr, "(%.2lf,%.2lf)-(%.2lf,%.2lf)-(%.2lf,%.2lf)\n",
			x(0, ps), y(0, ps), x(1, ps), y(1, ps), x(2, ps), y(2, ps));
	}
};

struct Mesh {
	typedef std::vector < Triangle > triangles_t;
	typedef std::vector < Point > points_t; 
	typedef std::vector < int > points_flags_t;
	triangles_t tr;
	points_t ps;
	//�����
	//0 - ���������� �����
	//1 - ����� �� ������� �������
	points_flags_t ps_flags;

	// ����� -> ������������ � �����
	std::vector < std::vector < int > > adj;
	// ������ ���������� �����
	std::vector < int > inner;
	// ������ ������� �����
	std::vector < int > outer;
	// ������������ ������ � ������� ps ������ � ������� inner ��� outer
	std::vector < int > p2io;

	void load(FILE * f);

	/**
	 * ������ �������
	 * ���������� ��� �������� ������������
	 */
	std::vector < Polynom > elem1(const Triangle & t) const;
	/**
	 * ���� �����, �� ���������� ������ �������� ���������� �����
	 */
	std::vector < Polynom > elem1_inner(const Triangle & t) const;

	/**
	 * ������ �������
	 * ���������� ������� ������������ � �������� �����
	 */
	Polynom elem1(const Triangle & t, int p) const;

	// ���� ����� ������� ������������ ������ ��������, �� ����
	// ������� ���� �������, ������� ����� ��������� � ������������
	// ���������� �����
	// ���� ��� �� ��������

	//!������ �������
	Polynom elem2(const Triangle & t);
	//!������ �������
	Polynom elem3(const Triangle & t);
};

/**
 * to global coordinates
 */
typedef double (* x_t)(double u, double v);

void print_function(FILE * to, double * ans, const Mesh & m, 
					x_t x = 0, x_t y = 0, x_t z = 0);

class Matrix;

/**
 * Callback. ���������� ��� ���� ������� phi_i, phi_j, ������������
 * � ����� point �� ������������ tr
 */
typedef double (* right_part_cb_t)
(const Polynom & phi_i, 
 const Polynom & phi_j, 
 int point, /* ����� ����� */
 int tr,    /* ����� ������������ */
 const Mesh & mesh,
 void * user_data /* ���� ���� ������� ����� ������, 
				     �������� �������� F �� ���� ����*/
 );

/**
 * ������� ������� �������.
 * �������� right_part_cb_t ��� ���� ������� phi_i, phi_j, ������������
 * � ����� ����� point �� ������������ tr
 */
void generate_right_part(double * b, const Mesh & m, right_part_cb_t right_part_cb, void * user_data);


/**
 * Callback. ���������� ��� ���� ������� phi_i, phi_j, ������������
 * � ����� point �� ������������ tr
 */
typedef double (* integrate_cb_t)
(const Polynom & phi_i, 
 const Polynom & phi_j, 
 int point, /* ����� ����� */
 int tr,    /* ����� ������������ */
 const Mesh & mesh,
 void * user_data /* ���� ���� ������� ����� ������ */
 );

/**
 * ������� ������� �������.
 * �������� integrate_cb ��� ���� ������� phi_i, phi_j, ������������
 * � ����� ����� point �� ������������ tr
 */
void generate_matrix(Matrix & A, const Mesh & m, integrate_cb_t integrate_cb, void * user_data);

/**
 * ������ �������
 */
void mke_solve(double * answer, const double * bnd, double * rp, Matrix & A, const Mesh & m);

/* ��������� ������� ������� */
void mke_p2u(double * p, const double * u, const double * bnd, const Mesh & m);

/* ������� ������� ������� */
void mke_u2p(double * u, const double * p, const Mesh & m);

void sphere_chafe_solve(double * Ans, const double * X0,
						const Mesh & m, const double * bnd, 
						double tau, double sigma, double mu);

#endif /* MKE_H */
