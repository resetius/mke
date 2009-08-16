/* -*- charset: utf-8 -*- */
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

namespace phelm_private_ 
{
	template < typename Matrix >
	void mat_add(Matrix & A, int i, int j, double a, bool transpose)
	{
		if (transpose) {
			A.add(j, i, (typename Matrix::data_type)a);
		} else {
			A.add(i, j, (typename Matrix::data_type)a);
		}
	}

	template < typename Matrix, typename Container >
	void mat_add(Matrix & A, int i, int j, const Container & c, bool transpose)
	{
		typename Container::const_iterator b = c.begin();
		typename Container::const_iterator e = c.end();
		typename Container::const_iterator it;

		if (transpose) {
			for (it = b; it != e; ++it)
			{
				A.add(it->j, it->i, (typename Matrix::data_type)it->a);
			}
		} else {
			for (it = b; it != e; ++it)
			{
				A.add(it->i, it->j, (typename Matrix::data_type)it->a);
			}
		}
	}

	template < typename T >
	void vec_add(T * b, int i, double a)
	{
		b[i] += (T) a;
	}

	template < typename T, typename Container >
	void vec_add(T * b1, int i, const Container & c)
	{
		typename Container::const_iterator b = c.begin();
		typename Container::const_iterator e = c.end();
		typename Container::const_iterator it;

		for (it = b; it != e; ++it)
		{
			b1[it->i] += (T)it->a;
		}
	}
};
