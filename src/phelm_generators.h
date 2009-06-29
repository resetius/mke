/* -*- charset: utf-8 -*- */
/* $Id$ */
/**
 * @file phelm_generators.h
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @section LICENSE
 *
 * <pre>
 * Copyright (c) 2009 Alexey Ozeritsky (Алексей Озерицкий)
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
 * </pre>
 *
 * @section DESCRIPTION
 * finite-element matrix generators, right-part generator, mesh-vector norm, 
 * mesh-vectors inner-product, mesh-vectors distance.
 */

#include "solver.h"
#include "util.h"
#include "phelm_private.h"

namespace phelm {

/**
 * Создает матрицу системы.
 * Вызывает integrate_cb для всех функций phi_i, phi_j, определенных
 * в общей точке point на треугольнике tr
 * если установлен флаг transpose, то генерит 
 * транспонированную матрицу
 * 
 * параметры callback'а:
 *
 *	const Polynom & phi_i, 
 *	const Polynom & phi_j, 
 *	const Triangle & tr,   номер треугольника 
 *	const Mesh & mesh,     сетка 
 *	int point_i,           глобальный номер точки 
 *	int point_j,           глобальный номер точки 
 *	int i,                 номер строки матрицы
 *	int j,                 номер столбца матрицы
 *	void * user_data       сюда могу входить любые данные
 *
 */
template < typename Functor, typename Data >
void generate_matrix(Matrix & A, const Mesh & m, 
					 Functor integrate_cb, 
					 Data user_data, 
					 bool transpose = false)
{
	using namespace phelm_private_;
	int rs  = (int)m.inner.size();     // размерность

	Timer t;
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) { // номер строки
		// по внутренним точкам
		int p = m.inner[i];

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) 
		{
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				int j    = m.p2io[p2]; // номер внутренней точки
				                       // номер столбца
				if (m.ps_flags[p2] == 1) {
					; // граница
				} else {
					mat_add(A, i, j, 
						integrate_cb(phi_i, phik[i0], trk, 
							m, p, p2, i, j, user_data), transpose);
				}
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_matrix: %lf \n", t.elapsed()); 
#endif
}

/**
 * Generate finite-elements matrix.
 */
template < typename Functor, typename Data >
void generate_full_matrix(Matrix & A, const Mesh & m, 
						  Functor integrate_cb, 
						  Data user_data,
						  bool transpose = false)
{
	using namespace phelm_private_;
	int sz  = (int)m.ps.size();

	Timer t;
//#pragma omp parallel for
	for (int p = 0; p < sz; ++p) {
		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				mat_add(A, p, p2, 
					integrate_cb(phi_i, phik[i0], trk, m, 
						p, p2, p, p2, user_data), transpose);
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_full_matrix: %lf \n", t.elapsed()); 
#endif
}

/**
 * Generate right part. 
 */
template < typename Functor, typename Data >
void generate_right_part(double * b, const Mesh & m, 
						 Functor right_part_cb, 
						 Data user_data)
{
	using namespace phelm_private_;
	int rs  = (int)m.inner.size();     // размерность
	Timer t;

	// WARNING: если генерируем правую часть для системы уравнений,
	// то реальная размерность не известна, поэтому 
	// memset(b, 0) надо вызывать руками до вызова generate_right_part !

#pragma omp parallel for
	for (int i = 0; i < rs; ++i)
	{
		// по внутренним точкам
		int p = m.inner[i];
		b[i]  = 0.0;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			// если граница не меняется по времени, то надо делать отдельный callback
			// для вычисления постоянной поправки к правой части?
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				vec_add(b, i, 
					right_part_cb(phi_i, phik[i0], 
						trk, m, p, p2, i, m.p2io[p2], user_data));
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_right_part: %lf \n", t.elapsed());
#endif
}

/**
 * Generate right-part.
 */
template < typename Functor, typename Data >
void generate_full_right_part(double * b, const Mesh & m, 
							  Functor right_part_cb, 
							  Data user_data)
{
	using namespace phelm_private_;
	int sz  = (int)m.ps.size();     // размерность

	Timer t;

	// WARNING: если генерируем правую часть для системы уравнений,
	// то реальная размерность не известна, поэтому 
	// memset(b, 0) надо вызывать руками до вызова generate_right_part !

#pragma omp parallel for
	for (int p = 0; p < sz; ++p)
	{
		b[p]  = 0.0;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольника в точке
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			// если граница не меняется по времени, то надо делать отдельный callback
			// для вычисления постоянной поправки к правой части?
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				vec_add(b, p, 
					right_part_cb(phi_i, phik[i0], 
						trk, m, p, p2, p, p2, user_data));
			}
		}
	}
#ifdef _DEBUG
	fprintf(stderr, "generate_right_part: %lf \n", t.elapsed());
#endif
}

/**
 * генерирует матрицу для интеграции краевых условий в правую часть
 * inner.size() x outer.size()
 * в cb передается phi_j где j точка границы
 * phi_i, где i внутренняя точка
 */
template < typename Functor, typename Data >
void generate_boundary_matrix(Matrix & A, const Mesh & m, 
							  Functor right_part_cb, 
							  Data user_data,
							  bool transpose = false)
{
	using namespace phelm_private_;
	int os = (int)m.outer.size(); // размер границы
	for (int j = 0; j < os; ++j) {
		// по внешним точкам
		int p2 = m.outer[j];

		for (uint tk = 0; tk < m.adj[p2].size(); ++tk) {
			// по треугольникам в точке
			int trk_j = m.adj[p2][tk];
			const Triangle & trk    = m.tr[trk_j];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_j           = trk.elem1(p2);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p    = m.tr[trk_j].p[i0];

				if (m.ps_flags[p] == 1) {
					;
				} else {
					// p - внутренняя точка
					int i = m.p2io[p];
					mat_add(A, i, j, 
						right_part_cb(phik[i0], phi_j, 
							trk, m, p, p2, i, j, user_data), transpose);
				}
			}
		}
	}
}

template < typename Functor, typename Data >
void convolution(double * ans, const double * u, const double * v, 
				 const Mesh & m, Functor cb, Data user_data)
{
	using namespace phelm_private_;

	int sz  = (int)m.ps.size(); // размерность

//#pragma omp parallel for
	for (int i = 0; i < sz; ++i)
	{
		// по всем точкам
		int p = i;
		ans[i] = 0.0;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// по треугольникам в точке
			int trk_i               = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int j  = trk.p[i0];
				ans[i] += u[i] * v[j] * 
					cb(phi_i, phik[i0], trk, m, i, j, i, j, user_data);
					//cb(phik[i0], phi_i, trk, m, i, j, i, j, user_data);
			}
		}
	}
}

/* сеточное скалярное произведение двух функций */
template < typename Functor, typename Data >
double scalar(const double * u, const double * v, const Mesh & m, 
				  Functor cb, Data user_data)
{
	int sz  = (int)m.ps.size();
	double s = 0.0;
	std::vector < double > nr(sz);
	convolution(&nr[0], u, v, m, cb, user_data);
//#pragma omp parallel for reduction(+:s)
	for (int i = 0; i < sz; ++i) {
		s = s + nr[i];
	}
	return s;
}

template < typename Functor, typename Data >
void generate_scalar_matrix(Matrix & mat, const Mesh & m, 
							Functor cb, Data user_data)
{
	generate_full_matrix(mat, m, cb, user_data);
}

/* сеточная норма */
template < typename Functor, typename Data >
double norm(const double * u, const Mesh & m, 
				Functor cb, Data user_data)
{
	return sqrt(scalar(u, u, m, cb, user_data));
}

inline double norm(const double * u, const Mesh & m)
{
	return norm(u, m, generic_scalar_cb, (void*)0);
}

/* сеточное расстояние */
template < typename Functor, typename Data>
double dist(const double * u, const double * v, const Mesh & m, 
				Functor cb, Data user_data)
{
	int sz  = (int)m.ps.size(); // размерность
	std::vector < double > diff(sz);
	vec_diff(&diff[0], u, v, sz);
	return norm(&diff[0], m, cb, user_data);
}

inline double dist(const double * u, const double * v, const Mesh & m)
{
	return dist(u, v, m, generic_scalar_cb, (void*)0);
}

}

