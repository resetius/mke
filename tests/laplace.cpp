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
#include <math.h>

#include "mke.h"
#include "util.h"
#include "solver.h"
#include "laplace.h"

using namespace std;

struct laplace_right_part_cb_data
{
	const double * F;
	const double * bnd;
};

static double 
laplace_right_part_cb( const Polynom & phi_i,
                       const Polynom & phi_j,
                       int point, /* номер точки */
                       int trk_i, /* номер треугольника */
                       const Mesh & m,
                       laplace_right_part_cb_data * d)
{
	const double * F = d->F;
	const Triangle & trk    = m.tr[trk_i];
	double b;

	if (m.ps_flags[point] == 1) { // на границе
		int j0       = m.p2io[point]; //номер внешней точки
		const double * bnd = d->bnd;
		Polynom poly = diff(phi_j, 0) * diff(phi_i, 0) 
			+ diff(phi_j, 1) * diff(phi_i, 1);
		b = - bnd[j0] * integrate(poly, trk, m.ps);
	} else {
		b = - F[point] * integrate(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

static double 
laplace_integrate_cb( const Polynom & phi_i,
                      const Polynom & phi_j, 
                      int point, /* номер точки */
                      int trk_i, /* номер треугольника */
                      const Mesh & m,
                      void * user_data)
{
	const Triangle & trk  = m.tr[trk_i];
	Polynom poly = diff(phi_j, 0) * diff(phi_i, 0) + diff(phi_j, 1) * diff(phi_i, 1);
	double a = integrate(poly, trk, m.ps);
	return a;
}

void laplace_solve(double * Ans, const Mesh & m, 
				   const double * F, const double * bnd)
{
	//пока используем первый порядок
	int sz  = m.ps.size();
	int ntr = m.tr.size();
	int rs  = m.inner.size();     //размерность

	vector < double > b(rs);      // правая часть
	vector < double > x(rs);      // ответ
	Matrix A(rs);

	Timer full;

	laplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = bnd;

	generate_matrix(A, m, laplace_integrate_cb, 0);
	generate_right_part(&b[0], m, (right_part_cb_t)(laplace_right_part_cb), (void*)&d);
	//A.print();
	//vector_print(&b[0], rs);
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed());
	mke_solve(Ans, bnd, &b[0], A, m);
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed()); 
}
