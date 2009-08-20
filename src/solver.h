#ifndef SOLVER_H
#define SOLVER_H
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

/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @section DESCRIPTION
 * The Matrix class and linear system solver.
 */

#include <vector>
#include <map>

#ifdef SPARSE
#ifndef GMRES
#define UMFPACK
#endif
#endif

#ifdef UMFPACK
#include <umfpack.h>
#endif

#ifdef SUPERLU
#include <slu_ddefs.h>
#endif

#include "util.h"
#include "gmres.h"

namespace phelm {

/**
 * @defgroup solver Linear equations solver.
 * Contains Linear equations solver class and
 * auxiliary routines.
 * @{
 */

/**
 * Matrix class.
 */
template < typename T >
class SparseMatrix {
protected:
	int n_;

	enum {
		CSR = 1,
		ELL = 2,
	};
	int format_;

	// разреженная матрица
	// формат хранения как в MatLab и UMFPACK
	// CSR or ELL format
	// for CSR only Ai and Ax are used !
	ArrayDevice < int > Ap_; // количества ненулевых в столбцах
	ArrayDevice < int > Ai_; // индексы (номер строки) ненулевых элементов
	ArrayDevice < T >   Ax_; // ненулевые элементы матрицы

	// for ELL format
	int cols_;
	int stride_;

	typedef std::map < int , T > row_t;
	// column number -> row_number -> value
	typedef std::vector < row_t > sparse_t;
	sparse_t A_;

	void make_sparse_csr();
	void make_sparse_ell();
	void make_sparse();

public:
	typedef T data_type;

	SparseMatrix(int n): n_(n), format_(ELL), Ap_(n + 1), A_(n)
	{
	}

	SparseMatrix(const int * Ap, const int * Ai, const T * Ax, int n, int nz): n_(n), format_(CSR)
	{
		Ap_.resize(n + 1); Ai_.resize(nz); Ax_.resize(nz);
		vec_copy(&Ap_[0], Ap, (n + 1));
		vec_copy(&Ai_[0], Ai, nz);
		vec_copy(&Ax_[0], Ax, nz);
	}

	/**
	 *  Add a number to element (i, j) (A[i][j] += a).
	 *  @param i - index
	 *  @param j - index
	 *  @param a - value
	 */
	void add(int i, int j, T a);

	/**
	 * Solve equation Ax = b.
	 * That function uses GMRES
	 * @param x - answer
	 * @param b - right part
	 */
	void solve(T * x, const T * b);

	/**
	 * Product of matrix by vector (out = A in).
	 * @param out - result
	 * @param in  - input vector
	 */
	void mult_vector(T * out, const T * in);

	/**
	 * print matrix to stdout.
	 */
	void print();
};

#ifdef UMFPACK

/**
 * Matrix class.
 */
template < typename T >
class UmfPackMatrix: public SparseMatrix < T >
{
	double Control_ [UMFPACK_CONTROL];
	double Info_ [UMFPACK_INFO];
	void *Symbolic_, *Numeric_ ;

	typedef SparseMatrix < T > base;
public:
	typedef T data_type;

	UmfPackMatrix(int n): SparseMatrix < T >(n), Symbolic_(0), Numeric_(0) 
	{
		base::format_ = base::CSR;
		umfpack_di_defaults(Control_);
	}

	UmfPackMatrix(const int * Ap, const int * Ai, const T * Ax, int n, int nz):
		SparseMatrix < T >(Ap, Ai, Ax, n, nz), Symbolic_(0), Numeric_(0)
	{
		umfpack_di_defaults(Control_);
	}

	~UmfPackMatrix()
	{
		umfpack_di_free_symbolic (&Symbolic_);
		umfpack_di_free_numeric (&Numeric_);
	}

	/**
	 * Solve equation Ax = b.
	 * That function uses UMFPACK
	 * @param x - answer
	 * @param b - right part
	 */
	void solve(T * x, const T * b);
};

template <>
class UmfPackMatrix < float >: public SparseMatrix < float >
{
	double Control_ [UMFPACK_CONTROL];
	double Info_ [UMFPACK_INFO];
	void *Symbolic_, *Numeric_ ;

	typedef SparseMatrix < float > base;

	std::vector < double > Ax_;

public:
	typedef float data_type;

	UmfPackMatrix(int n): SparseMatrix < float >(n), Symbolic_(0), Numeric_(0) 
	{
		base::format_ = base::CSR;
		umfpack_di_defaults(Control_);
	}

	UmfPackMatrix(const int * Ap, const int * Ai, const float * Ax, int n, int nz):
		SparseMatrix < float >(Ap, Ai, Ax, n, nz), Symbolic_(0), Numeric_(0)
	{
		umfpack_di_defaults(Control_);
	}

	~UmfPackMatrix()
	{
		umfpack_di_free_symbolic (&Symbolic_);
		umfpack_di_free_numeric (&Numeric_);
	}

	/**
	* Solve equation Ax = b.
	* That function uses UMFPACK
	* @param x - answer
	* @param b - right part
	*/
	void solve(float * x, const float * b)
	{
		if (base::Ax_.empty()) {
			base::make_sparse();
		}

		if (Ax_.empty()) {
			Ax_.resize(base::Ax_.size());

			for (int i = 0; i < (int)Ax_.size(); ++i)
			{
				Ax_[i] = (double)base::Ax_[i];
			}
		}

		int status = 0;

		if (Symbolic_ == 0) {
			status = umfpack_di_symbolic (base::n_, base::n_, 
				&base::Ap_[0], &base::Ai_[0], &Ax_[0], 
				&Symbolic_, Control_, Info_);
			assert(status == UMFPACK_OK);
		}

		if (Numeric_ == 0) {
			status = umfpack_di_numeric (&base::Ap_[0], &base::Ai_[0], 
				&Ax_[0], 
				Symbolic_, &Numeric_, Control_, Info_) ;
			assert(status == UMFPACK_OK);
		}

		std::vector < double > x1(base::n_);
		std::vector < double > b1(base::n_);

		for (int i = 0; i < base::n_; ++i)
		{
			b1[i] = (double)b[i];
		}

		status = umfpack_di_solve (UMFPACK_At, &base::Ap_[0], &base::Ai_[0], 
			&Ax_[0], &x1[0], &b1[0], Numeric_, Control_, Info_);
		assert(status == UMFPACK_OK);

		for (int i = 0; i < base::n_; ++i)
		{
			x[i] = (float)x1[i];
		}
	}
};
#endif

#ifdef SUPERLU
template < typename T >
class SuperLUMatrix: public SparseMatrix < T >
{
	typedef SparseMatrix < T > base;
	SuperMatrix A_, AC_, L_, U_, B_;

	std::vector < int > perm_c_;
	std::vector < int > perm_r_;
	std::vector < int > etree_;

public:
	typedef T data_type;

	SuperLUMatrix(int n): SparseMatrix < T >(n)
	{
		base::format_ = base::CSR;
	}

	SuperLUMatrix(const int * Ap, const int * Ai, const T * Ax, int n, int nz):
		SparseMatrix < T >(Ap, Ai, Ax, n, nz)
	{
	}

	~SuperLUMatrix()
	{
	}

	/**
	 * Solve equation Ax = b.
	 * That function uses SuperLU
	 * @param x - answer
	 * @param b - right part
	 */
	void solve(T * x, const T * b);
};
#endif

/**
 * Matrix class.
 */
template < typename T >
class SimpleMatrix
{
	int n_;
	Array < T, Allocator < T > > A_;  // матрица

public:
	typedef T data_type;

	SimpleMatrix(int n): n_(n), A_(n * n) {}

	/**
	 *  Add a number to element (i, j) (A[i][j] += a).
	 *  @param i - index
	 *  @param j - index
	 *  @param a - value
	 */
	void add(int i, int j, T a);

	/**
	 * Solve equation Ax = b.
	 * That function uses Gauss
	 * @param x - answer
	 * @param b - right part
	 */
	void solve(T * x, const T * b);

	/**
	 * Product of matrix by vector (out = A in).
	 * @param out - result
	 * @param in  - input vector
	 */
	void mult_vector(T * out, const T * in);

	/**
	 * print matrix to stdout.
	 */
	void print();
};

#if defined(UMFPACK) && !defined(GPGPU)
template < typename T >
class Matrix: public UmfPackMatrix < T >
{
public:
	Matrix(int n): UmfPackMatrix < T >(n) {}
};
#else
template < typename T >
class Matrix: public SparseMatrix < T > 
{
public:
	Matrix(int n): SparseMatrix < T > (n) {}
};
#endif

typedef Matrix < double > Matrixd;
typedef Matrix < float > Matrixf;

/**
 * Solve the system with A matrix (Ax=rp).
 * (Helper function) 
 * The function founds an answer on the inner part of the domain
 * and then sets boundary value of the answer to bnd    
 *
 * @param answer - the answer
 * @param bnd - boundary
 * @param rp - right part
 * @param A - the matrix of the system
 * @param m - mesh
 */
template < typename T, typename Matrix, typename Mesh >
void solve(T * answer, const T * bnd,
		   T * rp, Matrix & A, const Mesh & m)
{
	int sz  = (int)m.ps.size();
	int rs  = (int)m.inner.size();     // размерность
	Array < T, Allocator < T > > x(rs);      // ответ
	solve2(&x[0], rp, A, m);
	p2u(answer, &x[0], bnd, m);
}

/**
 * Solve the system with A matrix (Ax=rp).
 * (Helper function)
 * Found an answer on the inner part of the domain.    
 *
 * @param answer the answer
 * @param rp the right part
 * @param A the matrix of the system
 * @param m the mesh
 */
template < typename T, typename Matrix, typename Mesh >
void solve2(T * answer, T * rp, Matrix & A, const Mesh & m)
{
	int sz  = (int)m.ps.size();
	int rs  = (int)m.inner.size();     // размерность
	A.solve(answer, &rp[0]);
}

/** @} */ /* solver */

#include "impl/solver_impl.h"

}

#endif /* SOLVER_H */

