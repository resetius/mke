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

#ifdef SPARSE
#include <map>
#ifndef GMRES
#include <umfpack.h>
#endif
#endif

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
class Matrix {
	int n_;

#ifdef SPARSE
	// разреженная матрица
	// формат хранения как в MatLab и UMFPACK
	std::vector < int > Ap_;    // количества ненулевых в столбцах
	std::vector < int > Ai_;    // индексы (номер строки) ненулевых элементов
	std::vector < double > Ax_; // ненулевые элементы матрицы

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
	std::vector < double > A_;  // матрица
#endif

	/* заполняет Ap, Ai, Ax по промежуточным данным */
	void make_sparse();
	void add_sparse(int i, int j, double a);
	void solve_sparse(double * x, const double * b);

public:
	/**
	 * create matrix NxN.
	 * @param n - dimension
	 */
	Matrix(int n);
	~Matrix();

	/**
	 *  Add a number to element (i, j) (A[i][j] += a).
	 *  @param i - index
	 *  @param j - index
	 *  @param a - value
	 */
	void add(int i, int j, double a);

	/**
	 * Solve equation Ax = b.
	 * В зависимости от параметров компиляции используется либо Gauss
	 * либо GMRES либо UMFPACK.
	 * Если определен SPARSE и GMRES, то используется GMRES
	 * Если определен SPARSE, то используется UMFPACK
	 * Если не определен SPARSE, то используется Gauss
	 * @param x - answer
	 * @param b - right part
	 */
	void solve(double * x, const double * b);

	/**
	 * Product of matrix by vector (out = A in).
	 * @param out - result
	 * @param in  - input vector
	 */
	void mult_vector(double * out, const double * in);

	/**
	 * print matrix to stdout.
	 */
	void print();
};

struct Mesh;

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
void solve(double * answer, const double * bnd,
		   double * rp, Matrix & A, const Mesh & m);

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
void solve2(double * answer, double * rp, Matrix & A, const Mesh & m);

/** @} */ /* solver */

}

#endif /* SOLVER_H */

