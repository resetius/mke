#ifndef PHELM_LA_H
#define PHELM_LA_H
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
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef GPGPU
#include <cuda_runtime_api.h>
#endif

namespace phelm {

/**
 * @defgroup la Linear Algebra functions and Classes.
 * @{
 */

/**
 * The Gauss method. Solves a linear system Ax=b.
 *
 * @param A - the matrix of the system.
 * @param b - the right part.
 * @param x - the answer.
 * @param n - dimension.
 * @return 0 on success.
 */
int gauss (double *A, double *b, double *x, int n);

/**
 * Print NxN matrix to stdout.
 * @param A - matrix
 * @param n - dimension
 */
void mat_print(const double * A, int n);

/**
 * Print NxN matrix to stdout.
 * @param A - matrix
 * @param n - dimension
 */
void mat_print(const float * A, int n);

/**
 * Print vector to stdout.
 * @param A - vector
 * @param n - dimension
 */
void vec_print(const double * A, int n);

/**
 * Print vector to stdout.
 * @param A - vector
 * @param n - dimension
 */
void vec_print(const float * A, int n);

/**
 * Product of NxN matrix and vector
 * @param r - output vector, r = Ax
 * @param A - intput matrix
 * @param x - input vector
 * @param n - dimension of matrix and vector
 */
void mat_mult_vector(double * r, const double * A, const double * x, int n);

/**
 * Product of NxN matrix and vector
 * @param r - output vector, r = Ax
 * @param A - intput matrix
 * @param x - input vector
 * @param n - dimension of matrix and vector
 */
void mat_mult_vector(float * r, const float * A, const float * x, int n);

/**
 * Sparse Matrix structure.
 * UMFPACK format.
 * The matrix can be stored by rows or by columns.
 *
 * If matrix is stored by rows then  Ap[i+1]-Ap[i] is the number of nonzero entries in row i, 
 * Ai contains indices of nonzero entries of row.
 *
 * If matrix is stored by columns then  Ap[i+1]-Ap[i] is the number of nonzero entries in column i, 
 * Ai contains indices of nonzero entries of column.
 */
template < typename T >
struct Sparse_t {
	typedef T data_type;
	int * Ap;    ///< the number of elements in columns or rows
	int * Ai;    ///< column or row indices
	data_type * Ax; ///< holds nonzero matrix entires
};

typedef Sparse_t < double > Sparse;
typedef Sparse_t < float > Sparsef;

/**
 * If matrix is stored by columns then left multiply by vector: r = x A
 * If matrix is stored by rows then right multiply by vector: r = A x
 * @param r - output vector
 * @param A - sparse matrix
 * @param x - intput vector
 * @param n - the size of vector and matrix
 */
void sparse_mult_vector_l(double * r, const Sparse * A,
						  const double * x, int n);

void sparse_mult_vector_l(float * r, const Sparsef * A,
						  const float * x, int n);

/**
 * If matrix is stored by columns then right multiply by vector: r = A x
 * If matrix is stored by rows then left multiply by vector: r = x A
 * @param r - output vector
 * @param A - sparse matrix
 * @param x - intput vector
 * @param n - the size of vector and matrix
 */
void sparse_mult_vector_r(double * r, const Sparse * A,
						  const double * x, int n);

void sparse_mult_vector_r(float * r, const Sparsef * A,
						  const float * x, int n);

template < typename Sparse >
struct MatMultiplier 
{
	typedef typename Sparse::data_type data_type;

	static void mult_vector_l(data_type * r, const Sparse * A,
		const data_type * x, int n)
	{
		sparse_mult_vector_l(r, A, x, n);
	}

	static void mult_vector_r(data_type * r, const Sparse * A,
		const data_type * x, int n)
	{
		sparse_mult_vector_r(r, A, x, n);
	}
};

/**
 * Print sparse matrix to file.
 * @param A - input sparse matrix
 * @param n - dimension of sparse matrix
 * @param f - output file
 */
void sparse_print(const Sparse * A, int n, FILE * f);

void sparse_print(const Sparsef * A, int n, FILE * f);

/**
 * Linear combination of two vectors.
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param k1 - coefficient
 * @param k2 - coefficient
 * @param n - dimension of vectors
 * @return r = k1 * a + k2 * b
 */
void vec_sum1(double * r, const double * a,
			  const double *b, double k1, double k2, int n);

void vec_sum1(float * r, const float * a,
			  const float *b, float k1, float k2, int n);

/**
 * Linear combination of two vectors.
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param k2 - coefficient
 * @param n - dimension of vectors
 * @return r = a + k2 * b
 */
void vec_sum2(double * r, const double * a,
			  const double *b, double k2, int n);

void vec_sum2(float * r, const float * a,
			  const float *b, float k2, int n);

/**
 * Product of vector by number.
 * @param a - output vector
 * @param b - input vector
 * @param k - input number
 * @param n - dimension of vector
 * @return a = b * k
 */
void vec_mult_scalar(double * a, const double * b, double k, int n);
void vec_mult_scalar(float * a, const float * b, float k, int n);

/**
 * Sum of two vectors.
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 * @return r = a + b
 */
void vec_sum(double * r, const double * a, const double *b, int n);
void vec_sum(float * r, const float * a, const float *b, int n);

/**
 * Vector norm.
  \f[
  \sqrt{\sum_{i=0}v_i^2}
  \f]
 * @param v - input vector
 * @param n - dimension of vector
 * @return vector norm
 */
double vec_norm2(const double *v, int n);
float vec_norm2(const float *v, int n);

/**
 * Inner product of two vectors.
   \f[
   \sum_{i=0}^{n}a_i b_i
   \f]
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 * @return inner product of a and b
 */
double vec_scalar2(const double * a, const double * b, int n);
float vec_scalar2(const float * a, const float * b, int n);

/**
 * Element by element vector multiplication.
 * r = a * b
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 */
void vec_mult(double * r, const double * a, const double * b, int n);
void vec_mult(float * r, const float * a, const float * b, int n);

/**
 * Difference of two vectors.
 * @param r - the output vector
 * @param a - the input vector
 * @param b - the input vector
 * @param n - the dimension of vectors
 * @return r = a - b
 */
void vec_diff(double * r, const double * a, const double * b, int n);
void vec_diff(float * r, const float * a, const float * b, int n);

/**
 * Copy vector a to vector b.
 * @param b - the output vector
 * @param a - the input vector
 * @param n - the dimension of vectors
 */
template < typename T >
void vec_copy(T * b, const T * a, int n)
{
#ifndef GPGPU
	memcpy(b, a, n * sizeof(T));
#else
	cudaMemcpy(b, a, n * sizeof(T), cudaMemcpyDeviceToDevice);
#endif
}

template < typename T >
void vec_copy_from_host(T * b, const T * a, int n)
{
#ifndef GPGPU
	memcpy(b, a, n * sizeof(T));
#else
	cudaMemcpy(b, a, n * sizeof(T), cudaMemcpyHostToDevice);
#endif
}

template < typename T >
void vec_copy_from_device(T * b, const T * a, int n)
{
#ifndef GPGPU
	memcpy(b, a, n * sizeof(T));
#else
	cudaMemcpy(b, a, n * sizeof(T), cudaMemcpyDeviceToHost);
#endif
}

template < typename T, typename Alloc >
class Array {
	size_t size_;
	T * data_;
	Alloc alloc_;

public:
	Array(): size_(0), data_(0) {}
	Array(size_t size): size_(size), data_(0) 
	{
		data_ = alloc_.allocate(size_);
	}

	~Array() {
		if (data_) {
			alloc_.deallocate(data_, size_);
		}
	}

	void resize(size_t size) {
		if (size > size_) {
			T * p = alloc_.allocate(size);
			if (data_) {
				vec_copy(p, data_, size);
				alloc_.deallocate(data_, size_);
			}
			data_ = p;
			size_ = size;
		}
	}

	size_t size() const { return size_; }
	bool empty() const { return size_ == 0; }

	T & operator [] (int i) {
		assert(i < (int)size_);
		return data_[i];
	}

	const T & operator [] (int i) const {
		assert(i < (int)size_);
		return data_[i];
	}
};

int check_device_supports_double();

/**
 * @}
 */

} /* namespace */

#endif /* PHELM_LA_H */
