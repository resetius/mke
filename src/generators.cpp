/* -*- charset: utf-8 -*- */
/* $Id$ */

/* Copyright (c) 2015 Alexey Ozeritsky
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

#include "generators.h"
#include "solver.h"

using namespace phelm;

static void mat_add(
	SparseSolverAdder & A, const std::vector<Element> & Cont)
{
	for (auto it = Cont.begin(); it != Cont.end(); ++it) {
		A.add(it->i, it->j, it->a);
	}
}

void phelm::generate_matrixv(
	SparseSolverAdder & A, const Mesh & m,
	std::function<std::vector<Element>(
	const Polynom &,
	const Polynom &,
	const Triangle &,
	int,
	const Mesh &,
	int, int, int, int
	)> integrate_cb)
{
	int rs = (int)m.inner.size();    // размерность

	Timer t;
#pragma omp parallel for
	for (int i = 0; i < rs; ++i)   // номер строки
	{
		// по внутренним точкам
		int p = m.inner[i]; // get global point number
		int zone = m.tr[m.adj[p][0]].z;
		uint n_triangles = m.adj[p].size();

		assert(n_triangles > 0);
		assert(n_triangles > 1 || (m.order == 3 && m.tr[m.adj[p][0]].p[9] == p));

		for (uint tk = 0; tk < n_triangles; ++tk)
		{
			// по треугольникам в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk = m.tr[trk_i];
			const Triangle::basis_t & phik = trk.elem1(zone);
			const Polynom & phi_i = trk.elem1(p, zone);

			assert(phik.size() == trk.p.size());

			for (uint i0 = 0; i0 < phik.size(); ++i0)
			{
				int p2 = trk.p[i0];
				int j = m.p2io[p2]; // номер внутренней точки
				// номер столбца
				if (m.is_boundary(p2))
				{
					; // граница
				}
				else
				{
					mat_add(A, 
							integrate_cb(
								phi_i, phik[i0], trk, zone,
								m, p, p2, i, j));
				}
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_matrix: %lf \n", t.elapsed());
#endif
}

void phelm::generate_matrix(
	SparseSolverAdder & A, const Mesh & m,
	std::function<double(
	const Polynom &,
	const Polynom &,
	const Triangle &,
	int, const Mesh &,
	int, int, int, int
	)> integrate_cb)
{
	auto func = [integrate_cb](
		const Polynom & phi_i, const Polynom & phi_j,
		const Triangle & tr, int zone, const Mesh & m,
		int point_i, int point_j, int i, int j
		)
	{
		std::vector<Element> r;
		r.push_back(Element(i, j, 
			integrate_cb(phi_i, phi_j, tr, zone, m, point_i, point_j, i, j)));
		return r;
	};

	generate_matrixv(A, m, func);
}

void phelm::generate_boundary_matrixv(
	SparseSolverAdder & A, const Mesh & m,
	std::function<std::vector<Element>(
	const Polynom &,
	const Polynom &,
	const Triangle &,
	int,
	const Mesh &,
	int, int, int, int
	)> right_part_cb)
{
	int os = (int)m.outer.size(); // размер границы
	for (int j = 0; j < os; ++j)
	{
		// по внешним точкам
		int p2 = m.outer[j];
		int zone = m.tr[m.adj[p2][0]].z;

		for (uint tk = 0; tk < m.adj[p2].size(); ++tk)
		{
			// по треугольникам в точке
			int trk_j = m.adj[p2][tk];
			const Triangle & trk = m.tr[trk_j];
			const std::vector < Polynom > & phik = trk.elem1(zone);
			const Polynom & phi_j = trk.elem1(p2, zone);

			for (uint i0 = 0; i0 < phik.size(); ++i0)
			{
				int p = m.tr[trk_j].p[i0];

				if (m.is_boundary(p))
				{
					;
				}
				else
				{
					// p - внутренн€€ точка
					int i = m.p2io[p];
					mat_add(A,
							right_part_cb(phik[i0], phi_j,
							trk, zone, m, p, p2, i, j));
				}
			}
		}
	}
}

void phelm::generate_boundary_matrix(
	SparseSolverAdder & A, const Mesh & m,
	std::function<double(
	const Polynom &,
	const Polynom &,
	const Triangle &,
	int,
	const Mesh &,
	int, int, int, int
	)> integrate_cb)
{
	auto func = [integrate_cb](
		const Polynom & phi_i, const Polynom & phi_j,
		const Triangle & tr, int zone, const Mesh & m,
		int point_i, int point_j, int i, int j
		)
	{
		std::vector<Element> r;
		r.push_back(Element(i, j,
			integrate_cb(phi_i, phi_j, tr, zone, m, point_i, point_j, i, j)));
		return r;
	};

	generate_boundary_matrixv(A, m, func);
}
