#ifndef LAPL_H
#define LAPL_H
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

#include <vector>

#include "phelm.h"
#include "solver.h"

using phelm::Matrix;
using phelm::Mesh;
using phelm::Triangle;
using phelm::Polynom;

class SphereNorm {
public:
	const Mesh & m_;
	Matrix NORM_; // for fast norm calculation

	SphereNorm(const Mesh & m): m_(m), NORM_((int)m_.ps.size()) 
	{
		phelm::generate_full_matrix(NORM_, m_, phelm::sphere_scalar_cb, (void*)0);
	}

	double dist(const double * u, const double * v)
	{
		return phelm::fast_dist(u, v, m_, NORM_);
	}

	double norm(const double * u)
	{
		return phelm::fast_norm(u, m_, NORM_);
	}

	double scalar(const double * u, const double * v)
	{
		return phelm::fast_scalar(u, v, m_, NORM_);
	}
};

/**
 * ����������� ������
 */
class SphereLaplace {
public:
	Matrix idt_;
	Matrix laplace_;
	Matrix bnd1_; // L^-1
	Matrix bnd2_; // L^-1
	Matrix bnd3_; // L
	const Mesh & m_;

public:
	SphereLaplace(const Mesh & m);

	/**
	 * ������� �������� ������� ������� F �� ����������
	 * ������. � ������ ������� ������ ������ �������� �� bnd.
	 */
	void calc1(double * Ans, const double * F, const double * bnd);

	/**
	 * ������� �������� ������� ������� F �� ���������� ������.
	 * ���������� ������, ���������� ������ ���������� �����
	 */
	void calc2(double * Ans, const double * F);

/**
 * �������� �������� ������� �� ������� �������
 * @param Ans - �����
 * @param m - �����
 * @param F - ������ �����
 * @param bnd - ������� �������
 */
	void solve(double * Ans, const double * F, const double * bnd);
};

//@{ calculate right part or matrix elements
/* for sphere */
double slaplace(const Polynom & phi_i, const Polynom & phi_j, 
		const Triangle & trk, const Mesh::points_t & ps);

/* for plane */
double laplace(const Polynom & phi_i, const Polynom & phi_j,
		const Triangle & trk, const Mesh::points_t & ps);
//@}

class SphereChafe {
private:
	const Mesh & m_;
	SphereLaplace laplace_; /* ��������� */
	Matrix  A_;             /* ������� ����� ����� */

public:
	double tau_;
	double mu_;
	double sigma_;

public:
	SphereChafe(const Mesh & m, double tau, double sigma, double mu);
	~SphereChafe() {}

	void solve(double * Ans, const double * X0,
						const double * bnd, double t);
};

/**
 * ������ �� ������� �����
 */
class Laplace {
	Matrix idt_;      // inner
	Matrix laplace_;  // inner

//	Matrix fidt_;     // full
//	Matrix flaplace_; // full

	Matrix bnd1_; // L^-1
	Matrix bnd2_; // L^-1
	Matrix bnd3_; // L
	const Mesh & m_;

	void init_boundary_vectors();

public:
	Laplace(const Mesh & m);

	/**
	 * ������� �������� ������� ������� F �� ����������
	 * ������. � ������ ������� ������ ������ �������� �� bnd.
	 */
	void calc1(double * Ans, const double * F, const double * bnd);

	/**
	 * ������� �������� ������� ������� F �� ���������� ������.
	 * ���������� ������, ���������� ������ ���������� �����
	 */
	void calc2(double * Ans, const double * F);

	/**
	 * �������� �������� ������� �� ������� �������
	 * @param Ans - �����
	 * @param m - �����
	 * @param F - ������ �����
	 * @param bnd - ������� �������
	 */
	void solve(double * Ans, const double * F, const double * bnd);
};

class Chafe {
private:
	const Mesh & m_;
	Laplace laplace_; /* ��������� */
	Matrix A_;        /* ������� ����� ����� */

public:
	double tau_;
	double mu_;
	double sigma_;

public:
	Chafe(const Mesh & m, double tau, double sigma, double mu);
	~Chafe() {}

	void solve(double * Ans, const double * X0,
						const double * bnd, double t);
};

#endif /* LAPL_H */

