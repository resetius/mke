/* -*- charset: utf-8 -*- */
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

#include <assert.h>

template < typename T >
void SparseMatrix < T > ::mult_vector(T * out, const T * in)
{
	if (Ax_.empty()) {
		make_sparse();
	}

	SparseCSR < T > A;
	A.Ap = &Ap_[0];
	A.Ax = &Ax_[0];
	A.Ai = &Ai_[0];
	A.n  = n_;
	A.nz = (int)Ax_.size();

	sparse_mult_vector_r(out, A, in);
}

template < typename T >
void SimpleMatrix < T > ::mult_vector(T * out, const T * in)
{
	matrix_mult_vector(out, &A_[0], in);
}

template < typename T >
void SimpleMatrix < T > ::solve(T * x, const T * b)
{
	Timer t;

	gauss(&A_[0], &b[0], &x[0], n_);

#ifdef _DEBUG
	fprintf(stderr, "solver time: %lf\n", t.elapsed());
#endif
}

template < typename T >
void SparseMatrix < T > ::print()
{
	if (Ax_.empty()) {
		make_sparse();
	}

	SparseCSR < T > A;
	A.Ap = &Ap_[0];
	A.Ax = &Ax_[0];
	A.Ai = &Ai_[0];
	A.n  = n_;
	A.nz = (int)Ax_.size();

	sparse_print(A, stderr);

	// diag:
	{
		int i, i0, j;
		T * p = &Ax_[0];
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
}

template < typename T >
void SimpleMatrix < T > ::print()
{
	print_matrix(&A_[0], n_);
}

template < typename T >
void SparseMatrix < T > ::add(int i, int j, T a)
{
	// transposed !
	A_[i][j] += a;
}

template < typename T >
void SparseMatrix < T > ::make_sparse()
{
	int nz = 0; // non-null elements
	int idx = 0;
	for (uint i = 0; i < A_.size(); ++i) {
		nz += (int)A_[i].size();
	}
	Ax_.resize(nz);
	Ai_.resize(nz);

	std::vector < T > Ax(nz);
	std::vector < int > Ai(nz);
	std::vector < int > Ap(Ap_.size());

	Ap[0] = 0;

	for (uint i = 0; i < A_.size(); ++i)
	{
		Ap[i + 1] = Ap[i] + (int)A_[i].size();
		for (typename column_t::iterator it = A_[i].begin();
				it != A_[i].end(); ++it)
		{
			Ax[idx] = it->second;
			Ai[idx] = it->first;
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

	vec_copy_from_host(&Ax_[0], &Ax[0], nz);
	vec_copy_from_host(&Ai_[0], &Ai[0], nz);
	vec_copy_from_host(&Ap_[0], &Ap[0], n_ + 1);
}

template < typename T >
void SparseMatrix < T > ::solve(T * x, const T * b)
{
	if (Ax_.empty()) {
		make_sparse();
	}

	SparseCSR < T > A;
	A.Ap = &Ap_[0];
	A.Ax = &Ax_[0];
	A.Ai = &Ai_[0];
	A.n  = n_;
	A.nz = (int)Ax_.size();

	gmres(&x[0], &A, &b[0], csr_mult_vector < T >, n_, 100, 1000);
}

#ifdef UMFPACK
template < typename T >
void UmfPackMatrix < T > ::solve(T * x, const T * b)
{
	if (base::Ax_.empty()) {
		base::make_sparse();
	}

	int status;

	if (Symbolic_ == 0) {
		status = umfpack_di_symbolic (base::n_, base::n_, 
				&base::Ap_[0], &base::Ai_[0], &base::Ax_[0], 
				&Symbolic_, Control_, Info_);
		assert(status == UMFPACK_OK);
	}

	if (Numeric_ == 0) {
		status = umfpack_di_numeric (&base::Ap_[0], &base::Ai_[0], 
				&base::Ax_[0], 
				Symbolic_, &Numeric_, Control_, Info_) ;
		assert(status == UMFPACK_OK);
	}

	status = umfpack_di_solve (UMFPACK_At, &base::Ap_[0], &base::Ai_[0], 
			&base::Ax_[0], x, b, Numeric_, Control_, Info_);
	assert(status == UMFPACK_OK);
}

#endif

