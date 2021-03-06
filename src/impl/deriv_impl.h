/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2011 Alexey Ozeritsky
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

#include "phelm.h"
#include "array.h"

namespace deriv_impl 
{

	using namespace std;
	using namespace phelm;

	inline double id_cb (
		const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		int z,
		const Mesh & m,
		int point_i,
		int point_j,
		int i,
		int j,
		void *)
	{
		return integrate (phi_i * phi_j, trk, z);
	}

	template < typename T >
	struct diff_cb
	{
		int var_;

		diff_cb(int var): var_(var) {}

		double operator () ( 
			const Polynom & phi_i,
			const Polynom & phi_j,
			const Triangle & trk,
			int z,
			const Mesh & m,
			int point_i,
			int point_j,
			int i,
			int j,
			const T * u) const
		{
			Polynom poly = diff (phi_j, var_) * phi_i;
			double r = integrate(poly, trk, z);
			r *= u[point_j];
			/*if (m.ps_flags[point_j] == 1) {
				r += u[point_j] * 
					integrate(diff(phi_j * phi_i, var_),
					trk, m.ps);
			}*/
			return r;
		}
	};
}

template < typename T, typename Matrix >
Deriv < T, Matrix > ::Deriv (const Mesh & m): 
	m_(m), 
	idt_(m.size)
{
	generate_full_matrix (idt_, m, deriv_impl::id_cb, (void*) 0);
}

template < typename T, typename Matrix >
void Deriv < T, Matrix > ::calc_x(T * Ans, const T * u)
{
	//ArrayHost < T > rp(m_.size);
	std::vector < T > rp(m_.size);
	generate_full_right_part(&rp[0], m_, deriv_impl::diff_cb<T>(0), u);
	idt_.solve(&Ans[0], &rp[0]);
}

template < typename T, typename Matrix >
void Deriv < T, Matrix > ::calc_y(T * Ans, const T * u)
{
	std::vector < T > rp(m_.size);
	generate_full_right_part(&rp[0], m_, deriv_impl::diff_cb<T>(1), u);
	idt_.solve(&Ans[0], &rp[0]);
}

