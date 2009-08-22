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

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::load(std::vector < row_t > & A)
{
	int idx = 0;
	n_  = A.size();
	nz_ = 0; // non-null elements
	for (uint i = 0; i < A.size(); ++i) {
		nz_ += (int)A[i].size();
	}

	Ax_.resize(nz_);
	Ai_.resize(nz_);
	Ap_.resize(n_);

	std::vector < T >   Ax(nz_);
	std::vector < int > Ai(nz_);
	std::vector < int > Ap(n_ + 1);

	Ap[0] = 0;

	for (uint i = 0; i < A.size(); ++i)
	{
		Ap[i + 1] = Ap[i] + (int)A[i].size();
		for (typename row_t::iterator it = A[i].begin();
				it != A[i].end(); ++it)
		{
			Ax[idx] = it->second;
			Ai[idx] = it->first;
			idx += 1;
		}

		{
			row_t t;
			A[i].swap(t);
		}
	}

	{
		sparse_t t;
		A.swap(t);
	}

	vec_copy_from_host(&Ax_[0], &Ax[0], nz_);
	vec_copy_from_host(&Ai_[0], &Ai[0], nz_);
	vec_copy_from_host(&Ap_[0], &Ap[0], n_ + 1);
}

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::mult(T * r, const T * x) const
{
	sparse_mult_vector_r(r, &Ap_[0], &Ai_[0], &Ax_[0], x, n_, nz_);
}

template < typename T, template < class > class Alloc >
void StoreELL < T, Alloc > ::load(std::vector < row_t > & A)
{
	nz_ = 0; // non-null elements
	n_  = A.size();
	for (uint i = 0; i < A.size(); ++i) {
		nz_ += (int)A[i].size();
	}
	cols_ = 0;
	for (uint i = 0; i < A.size(); ++i) {
		cols_ = std::max(cols_, (int)A[i].size());
	}
	stride_ = 32 * ((n_ + 32 - 1) / 32);

	Ax_.resize(cols_ * stride_);
	Ai_.resize(cols_ * stride_);

	std::vector < T > Ax(cols_ * stride_);
	std::vector < int > Ai(cols_ * stride_);

	for (uint i = 0; i < A.size(); ++i)
	{
		int idx = 0;
		for (typename row_t::iterator it = A[i].begin();
				it != A[i].end(); ++it)
		{
			Ax[idx * stride_ + i] = it->second;
			Ai[idx * stride_ + i] = it->first;
			idx++;
		}

		{
			row_t t;
			A[i].swap(t);
		}
	}

	{
		sparse_t t;
		A.swap(t);
	}

	vec_copy_from_host(&Ax_[0], &Ax[0], cols_ * stride_);
	vec_copy_from_host(&Ai_[0], &Ai[0], cols_ * stride_);
}

template < typename T, template < class > class Alloc >
void StoreELL < T, Alloc > ::mult(T * r, const T * x) const
{
	sparse_mult_vector_r(r, &Ai_[0], &Ax_[0], x, n_, cols_, stride_);
}

template < typename T >
void SimpleSolver < T > ::mult_vector(T * out, const T * in)
{
	matrix_mult_vector(out, &A_[0], in);
}

template < typename T >
void SimpleSolver < T > ::solve(T * x, const T * b)
{
	Timer t;

	gauss(&A_[0], &b[0], &x[0], n_);

#ifdef _DEBUG
	fprintf(stderr, "solver time: %lf\n", t.elapsed());
#endif
}

template < typename T >
void SimpleSolver < T > ::print()
{
	print_matrix(&A_[0], n_);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::mult_vector(T * out, const T * in)
{
	if (store_.mult.empty()) {
		store_.mult.load(A_);
	}

	store_.mult.mult(out, in);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::print()
{
	// implement !
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::add(int i, int j, T a)
{
	// transposed !
	A_[i][j] += a;
}

template < typename T, typename Invert >
void SparseSolver__Ax(T * r, const Invert * invert, const T * x, int n)
{
	invert->mult(r, x);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::solve(T * x, const T * b)
{
	if (store_.invert.empty()) {
		store_.invert.load(A_);
	}

	gmres(&x[0], &store_.invert, &b[0], SparseSolver__Ax < T, typename store_t::invert_t > , 
		store_.invert.n_, 100, 1000);
}

#ifdef SUPERLU

template < typename T >
struct SuperLu
{
	static void gstrf (superlu_options_t *options, SuperMatrix *A,
		    int relax, int panel_size, int *etree, void *work, int lwork,
			int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
			SuperLUStat_t *stat, int *info)
	{
		assert(0);
	}

	static void gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
        int *perm_c, int *perm_r, SuperMatrix *B,
        SuperLUStat_t *stat, int *info)
	{
		assert(0);
	}

	static void
	Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, 
		       T *nzval, int *rowind, int *colptr,
		       Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
		assert(0);
	}

	static void
	Create_Dense_Matrix(SuperMatrix *X, int m, int n, T *x, int ldx,
		    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
		assert(0);
	}
};

extern "C" void
sgstrf (superlu_options_t *options, SuperMatrix *A,
        int relax, int panel_size, int *etree, void *work, int lwork,
        int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
        SuperLUStat_t *stat, int *info);

extern "C" void
psgstrf (superlu_options_t *options, SuperMatrix *A,
        int relax, int panel_size, int *etree, void *work, int lwork,
        int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
        SuperLUStat_t *stat, int *info);

extern "C" void
pdgstrf (superlu_options_t *options, SuperMatrix *A,
        int relax, int panel_size, int *etree, void *work, int lwork,
        int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
        SuperLUStat_t *stat, int *info);

extern "C" void 
sgstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
        int *perm_c, int *perm_r, SuperMatrix *B,
        SuperLUStat_t *stat, int *info);

extern "C" void
sCreate_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, 
		       float *nzval, int *rowind, int *colptr,
		       Stype_t stype, Dtype_t dtype, Mtype_t mtype);

extern "C" void
sCreate_Dense_Matrix(SuperMatrix *X, int m, int n, float *x, int ldx,
		    Stype_t stype, Dtype_t dtype, Mtype_t mtype);

#ifdef SUPERLU_MT
#define MT_CALL(func) \
	p##func
#else
#define MT_CALL(func) \
	func
#endif

template < >
struct SuperLu < double >
{
	static void gstrf (superlu_options_t *options, SuperMatrix *A,
		    int relax, int panel_size, int *etree, void *work, int lwork,
			int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
			SuperLUStat_t *stat, int *info)
	{
		MT_CALL(dgstrf)(options, A, relax, panel_size, etree, work, lwork, perm_c, 
			perm_r, L, U, stat, info);
	}

	static void gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
        int *perm_c, int *perm_r, SuperMatrix *B,
        SuperLUStat_t *stat, int *info)
	{
		dgstrs(trans, L, U, perm_c, perm_r, B, stat, info);
	}

	static void
	Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, 
		       double *nzval, int *rowind, int *colptr,
		       Stype_t stype, Mtype_t mtype)
	{
		dCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr, stype, SLU_D, mtype);
	}

	static void
	Create_Dense_Matrix(SuperMatrix *X, int m, int n, double *x, int ldx,
		    Stype_t stype, Mtype_t mtype)
	{
		dCreate_Dense_Matrix(X, m, n, x, ldx, stype, SLU_D, mtype);
	}
};

template < >
struct SuperLu < float >
{
	static void gstrf (superlu_options_t *options, SuperMatrix *A,
		    int relax, int panel_size, int *etree, void *work, int lwork,
			int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
			SuperLUStat_t *stat, int *info)
	{
		MT_CALL(sgstrf)(options, A, relax, panel_size, etree, work, lwork, perm_c, 
			perm_r, L, U, stat, info);
	}

	static void gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
        int *perm_c, int *perm_r, SuperMatrix *B,
        SuperLUStat_t *stat, int *info)
	{
		sgstrs(trans, L, U, perm_c, perm_r, B, stat, info);
	}

	static void
	Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz, 
		       float *nzval, int *rowind, int *colptr,
		       Stype_t stype, Mtype_t mtype)
	{
		sCreate_CompCol_Matrix(A, m, n, nnz, nzval, rowind, colptr, stype, SLU_S, mtype);
	}

	static void
	Create_Dense_Matrix(SuperMatrix *X, int m, int n, float *x, int ldx,
		    Stype_t stype, Mtype_t mtype)
	{
		sCreate_Dense_Matrix(X, m, n, x, ldx, stype, SLU_S, mtype);
	}
};

template < typename T >
void SuperLUMatrix < T > ::solve(T * x, const T * b)
{
	if (base::Ax_.empty()) {
		fprintf(stderr, "make sparse\n");
		base::make_sparse();
	}

	SuperMatrix B;
	SuperLUStat_t stat;
	StatInit(&stat);
	int info;

	if (perm_c_.empty()) {
		int panel_size;
		int relax;

		SuperLu < T >::Create_CompCol_Matrix(&A_, (int)base::n_, (int)base::n_, (int)base::Ax_.size(),
							   &base::Ax_[0], &base::Ai_[0], &base::Ap_[0],
							   SLU_NC, SLU_GE);

		//this->print();
		//dPrint_CompCol_Matrix("A", &A_);

		superlu_options_t options;
		set_default_options(&options);

		perm_c_.resize(base::n_);
		perm_r_.resize(base::n_);
		etree_.resize(base::n_);

		get_perm_c(3, &A_, &perm_c_[0]);
		sp_preorder(&options, &A_, &perm_c_[0], &etree_[0], &AC_);

		panel_size = sp_ienv(1);
		relax      = sp_ienv(2);

		SuperLu < T > ::gstrf(&options, &AC_, relax, panel_size, &etree_[0], NULL, 0, 
			&perm_c_[0], &perm_r_[0], &L_, &U_, &stat, &info);

		assert (info == 0);

		NCPformat *store;
		store = (NCPformat *)L_.Store;
		fprintf(stderr, "  nz: %d, n: %d\n", (int)base::Ax_.size(), (int)base::n_);
		fprintf(stderr, "        : %lf%%\n", (double)base::Ax_.size() /(double)base::n_ /(double)base::n_ * 100.0);
		fprintf(stderr, "L nz: %d: %lf%%\n", store->nnz, (double)store->nnz /(double)base::n_ /(double)base::n_ * 100.0);
		store = (NCPformat *)U_.Store;
		fprintf(stderr, "U nz: %d: %lf%%\n", store->nnz, (double)store->nnz /(double)base::n_ /(double)base::n_ * 100.0);
	}

	memcpy(x, b, base::n_ * sizeof(T));
	SuperLu < T >::Create_Dense_Matrix(&B, base::n_, 1, x, base::n_, SLU_DN, SLU_GE);

	SuperLu < T > ::gstrs(TRANS, &L_, &U_, &perm_c_[0], &perm_r_[0], &B, &stat, &info);
	
	assert(info == 0);

	Destroy_SuperMatrix_Store(&B);
	StatFree(&stat);
}
#endif
