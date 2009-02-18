/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (јлексей ќзерицкий)
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

using namespace std;

void laplace_solve(double * Ans, const Mesh & m, double * F, double * bnd)
{
	//пока используем первый пор€док
	int sz  = m.ps.size();
	int ntr = m.tr.size();
	int rs  = m.inner.size();     //размерность

	vector < double > b(rs);      // права€ часть
	vector < double > x(rs);      // ответ
	Matrix A(rs);

	Timer full;

	fprintf(stderr, "build 1: ");
	Timer build;

#pragma omp parallel for
	for (int i = 0; i < rs; ++i)
	{
		// по внутренним точкам
		int p = m.inner[i];
		for (size_t tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			Polynom phi_i           = m.elem1(trk, p);
			vector < Polynom > phik = m.elem1_inner(trk);

			for (size_t i0 = 0; i0 < phik.size(); ++i0) {
				int p2 = m.tr[trk_i].p[i0];
				//Polynom tmp = phi_i * phik[i0];
				b[i] += - F[p2] * integrate(phi_i * phik[i0], trk, m.ps);
				//b[i] += - F[p2] * integrate(phik[i0], trk, m.ps);
			}
		}
	}

	fprintf(stderr, "%lf \n", build.elapsed()); build.restart();
	//exit(0);
	fprintf(stderr, "build 2: ");

#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		// по внутренним точкам
		int p = m.inner[i];
		for (size_t tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			Polynom phi_i           = m.elem1(trk, p);
			vector < Polynom > phik = m.elem1(trk);

			for (size_t i0 = 0; i0 < phik.size(); ++i0) {
				int p2 = m.tr[trk_i].p[i0];
				if (m.ps_flags[p2] == 1) {
					//в случае ненулевых краевых условий это должно уйти в b[i]
					int j0 = m.p2io[p2]; //номер внешней точки
					Polynom poly = diff(phik[i0], 0) * diff(phi_i, 0)
						+ diff(phik[i0], 1) * diff(phi_i, 1);
					b[i] -= bnd[j0] * integrate(poly, trk, m.ps);
				} else {
					int j = m.p2io[p2]; //номер внутренней точки

					Polynom poly = diff(phik[i0], 0) * diff(phi_i, 0) 
						+ diff(phik[i0], 1) * diff(phi_i, 1);
					double a = integrate(poly, trk, m.ps);
					A.add(i, j, a);
				}
			}
		}
	}
	fprintf(stderr, "%lf \n", build.elapsed()); 
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed()); 

//	printf("matrix:\n");
//	matrix_print(&A[0], rs);

//	printf("vector:\n");
//	vector_print(&b[0], rs);

	fprintf(stderr, "solve %dx%d: \n", rs, rs);
//	A.print();

	Timer solve;
	A.solve(&b[0], &x[0]);
	fprintf(stderr, "done: %lf \n", solve.elapsed());
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed()); 

	for (int i = 0; i < sz; ++i) {
		if (m.ps_flags[i] == 1) {
			//внешн€€
			Ans[i] = bnd[m.p2io[i]];
		} else {
			//внутренн€€
			Ans[i] = x[m.p2io[i]];
		}
	}
}

