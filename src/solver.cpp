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

#include <assert.h>
#include <stdio.h>

#include "solver.h"
#include "util.h"
#include "gmres.h"

typedef unsigned int uint;

#ifdef SPARSE
Matrix::Matrix(int n): n_(n), Ap_(n_ + 1), A_(n_)
#ifndef GMRES
	, Symbolic_(0), Numeric_(0)
#endif
{
#ifndef GMRES
	umfpack_di_defaults(Control_);
#endif
}
#else
Matrix::Matrix(int n): n_(n), A_(n_ * n_) {}
#endif

Matrix::~Matrix()
{
#if defined(SPARSE) && !defined(GMRES)
	umfpack_di_free_symbolic (&Symbolic_);
	umfpack_di_free_numeric (&Numeric_);
#endif
}

void Matrix::add(int i, int j, double a) {
#ifdef SPARSE
	add_sparse(i, j, a);
#else
	A_[i * n_ + j] += a;
#endif
}

void Matrix::mult_vector(double * out, const double * in)
{
#ifdef SPARSE
	if (Ax_.empty()) {
		make_sparse();
	}
	Sparse A;
	A.Ap = &Ap_[0];
	A.Ax = &Ax_[0];
	A.Ai = &Ai_[0];

	sparse_mult_vector_l(out, &A, in, n_);
#else
	matrix_mult_vector(out, &A_[0], in);
#endif
}

void Matrix::solve(double * x, const double * b)
{
	Timer t;

#ifdef SPARSE
	solve_sparse(x, b);
#else

#ifdef GMRES
	gmres(&x[0], &A_[0], &b[0], (Ax_t)matrix_mult_vector, n_, 100, 1000);
#else
	gauss(&A_[0], &b[0], &x[0], n_);
#endif

#endif

	fprintf(stderr, "solver time: %lf\n", t.elapsed());
}

void Matrix::print()
{
#ifdef SPARSE
	if (Ax_.empty()) {
		make_sparse();
	}

	Sparse A;
	A.Ap = &Ap_[0];
	A.Ax = &Ax_[0];
	A.Ai = &Ai_[0];
	sparse_print(&A, n_, stderr);

	// diag:
	{
		int i, i0, j;
		double * p = &Ax_[0];
		fprintf(stderr, "diag:\n");
		for (j = 0; j < n_; ++j) {
			for (i0 = Ap_[j]; i0 < Ap_[j + 1]; ++i0, ++p) {
				i = Ai_[i0];
				if (i == j) {
					fprintf(stderr, "%8.3le ", *p);
				}
			}
		}
		fprintf(stderr, "\n");
	}
#else
	print_matrix(&A_[0], n);
#endif
}

#ifdef SPARSE

void Matrix::add_sparse(int i, int j, double a)
{
	// transposed !
	A_[i][j] += a;
}

void Matrix::make_sparse()
{
	int nz = 0; // non-null elements
	int idx = 0;
	for (uint i = 0; i < A_.size(); ++i) {
		nz += A_[i].size();
	}
	Ax_.resize(nz);
	Ai_.resize(nz);

	Ap_[0] = 0;
	for (uint i = 0; i < A_.size(); ++i)
	{
		Ap_[i + 1] = Ap_[i] + A_[i].size();
		for (column_t::iterator it = A_[i].begin();
				it != A_[i].end(); ++it)
		{
			Ax_[idx] = it->second;
			Ai_[idx] = it->first;
			idx += 1;
		}

		{
			column_t t;
			A_[i].swap(t);
		}
	}

	{
		sparse_t t;
		A_.swap(t);
	}
}

void Matrix::solve_sparse(double * x, const double * b)
{
	if (Ax_.empty()) {
		make_sparse();
	}
#ifdef GMRES
	Sparse A;
	A.Ap = &Ap_[0];
	A.Ax = &Ax_[0];
	A.Ai = &Ai_[0];

	gmres(&x[0], &A, &b[0], (Ax_t)sparse_mult_vector_l, n_, 100, 1000);
#else // UMFPACK
	int status;

	if (Symbolic_ == 0) {
		status = umfpack_di_symbolic (n_, n_, &Ap_[0], &Ai_[0], &Ax_[0], &Symbolic_, Control_, Info_);
		assert(status == UMFPACK_OK);
	}

	if (Numeric_ == 0) {
		status = umfpack_di_numeric (&Ap_[0], &Ai_[0], &Ax_[0], Symbolic_, &Numeric_, Control_, Info_) ;
		assert(status == UMFPACK_OK);
	}

	status = umfpack_di_solve (UMFPACK_At, &Ap_[0], &Ai_[0], &Ax_[0], x, b, Numeric_, Control_, Info_);
	assert(status == UMFPACK_OK);
#endif
}

#endif

