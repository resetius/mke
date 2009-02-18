#ifndef SOLVER_H
#define SOLVER_H
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

#ifdef SPARSE
#include <map>
#ifndef GMRES
#include <umfpack.h>
#endif
#endif

/**
 * �������� Ax = b
 */
class Matrix {
	int n_;

#ifdef SPARSE
	// ����������� �������
	// ������ �������� ��� � MatLab � UMFPACK
	std::vector < int > Ap_;    // ���������� ��������� � ��������
	std::vector < int > Ai_;    // ������� (����� ������) ��������� ���������
	std::vector < double > Ax_; // ��������� �������� �������

	//column_number -> row_number -> value
	typedef std::map < int , double > column_t;
	typedef std::vector < column_t > sparse_t;
	sparse_t A_;

#ifndef GMRES
	double Control_ [UMFPACK_CONTROL];
	double Info_ [UMFPACK_INFO];
	void *Symbolic_, *Numeric_ ;
#endif

#else
	std::vector < double > A_;  // �������
#endif

	//��������� Ap, Ai, Ax �� ������������� ������
	void make_sparse();
	void add_sparse(int i, int j, double a);
	void solve_sparse(double * b, double * x);

public:
	Matrix(int n);
	~Matrix();

	// A[i][j] += a;
	void add(int i, int j, double a);

	// Ax = b
	void solve(double * b, double * x);

	void print();
};

#endif /* SOLVER_H */
