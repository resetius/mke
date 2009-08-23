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

#include "util.h"
#include "gmres.h"

namespace phelm {

/**
 * @defgroup solver Linear equations solver.
 * Contains Linear equations solver class and
 * auxiliary routines.
 * @{
 */

template < typename T, template < class > class Alloc = Allocator >
struct StoreCSR
{
	typedef T data_type;
	int n_;
	int nz_;
	Array < int, Alloc < int > > Ap_; // количества ненулевых в столбцах
	Array < int, Alloc < int > > Ai_; // индексы (номер строки) ненулевых элементов
	Array < T, Alloc < T > > Ax_;   // ненулевые элементы матрицы

	StoreCSR(): n_(0), nz_(0) {}
	StoreCSR(int n, int nz): n_(n), nz_(nz), Ap_(n_ + 1), Ai_(nz_), Ax_(nz_) {}

	void resize(int n, int nz)
	{
		n_  = n;
		nz_ = nz;
		Ap_.resize(n_ + 1);
		Ai_.resize(nz_);
		Ax_.resize(nz_);
	}

	bool empty() 
	{
		return n_ == 0;
	}

	typedef std::map < int , T > row_t;
	typedef std::vector < row_t > sparse_t;

	/**
	 * Fill CSR Matrix from unstructured data
	 */
	void load(const sparse_t & unstruct);

	void mult(T * r, const T * x) const;
};

template < typename T, template < class > class Alloc = Allocator >
struct StoreELL
{
	typedef T data_type;
	int n_;
	int nz_;
	int cols_;
	int stride_;

	Array < int, Alloc < int > > Ai_;
	Array < T, Alloc < T > > Ax_;

	StoreELL(): n_(0), nz_(0), cols_(0), stride_(0) {}

	StoreELL(int n, int nz, int cols): 
		n_(n), nz_(nz), cols_(cols), stride_(32 * ((cols_ + 32 - 1) / 32)),
		Ai_(cols_ * stride_), Ax_(cols_ * stride_)
	{
	}

	void resize(int n, int nz, int cols)
	{
		n_    = n;
		nz_   = nz;
		cols_ = cols;
		Ai_.resize(cols_ * stride_);
		Ax_.resize(cols_ * stride_);
	}

	bool empty() 
	{
		return n_ == 0;
	}

	typedef std::map < int , T > row_t;
	typedef std::vector < row_t > sparse_t;

	/**
	 * Fill ELL Matrix from unstructured data
	 */
	void load(const sparse_t & unstruct);

	void mult(T * r, const T * x) const;
};

template < typename Store1, typename Store2 >
struct DoubleStore
{
	typedef Store1 mult_t;
	typedef Store2 invert_t;

	Store1 mult;
	Store2 invert;
};

template < typename Store >
struct DoubleStore < Store, Store >
{
	typedef Store mult_t;
	typedef Store invert_t;

	Store both;
	Store & mult;
	Store & invert;
	DoubleStore (): mult(both), invert(both) {}
};

/**
 * Solver class.
 * В солвере может быть два контейнера с разными аллокаторами для обращения
 * и для умножения.
 */
template < typename T, typename MultStore, typename InvStore = MultStore >
class SparseSolver {
protected:
	typedef DoubleStore < MultStore, InvStore > store_t;
	store_t store_;
	typedef std::map < int , T > row_t;
	typedef std::vector < row_t > sparse_t;

	sparse_t A_;

public:
	typedef T data_type;

	SparseSolver(int n): A_(n)
	{
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
// implements UmfPackSolver
#include "impl/solver_umfpack.h"
#endif

#ifdef SUPERLU
// implements SuperLUSolver
#include "impl/solver_superlu.h"
#endif

/**
 * Matrix class.
 */
template < typename T >
class SimpleSolver
{
	int n_;
	Array < T, Allocator < T > > A_;  // матрица

public:
	typedef T data_type;

	SimpleSolver(int n): n_(n), A_(n * n) {}

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


#if defined(SUPERLU) && !defined(GPGPU)

template < typename T >
class Solver: public SuperLUSolver < T, StoreELL < T , Allocator >  >
{
	typedef SuperLUSolver < T, StoreELL < T , Allocator >  > base;
public:
	Solver(int n): base(n) {}
};

#elif defined(UMFPACK) && !defined(GPGPU)

template < typename T >
class Solver: public UmfPackSolver < T, StoreELL < T , Allocator >  >
{
	typedef UmfPackSolver < T, StoreELL < T , Allocator >  > base;
public:
	Solver(int n): base(n) {}
};

#else
template < typename T >
class Solver: public SparseSolver < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > > 
{
	typedef SparseSolver < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > >  base;
public:
	Solver(int n): base (n) {}
};
#endif

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

