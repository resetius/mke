#include "solver.h"
#include "util.h"
#include "mke_private.h"

namespace MKE {

/**
 * ������� ������� �������.
 * �������� integrate_cb ��� ���� ������� phi_i, phi_j, ������������
 * � ����� ����� point �� ������������ tr
 * ���� ���������� ���� transpose, �� ������� 
 * ����������������� �������
 * 
 * ��������� callback'�:
 *
 *	const Polynom & phi_i, 
 *	const Polynom & phi_j, 
 *	const Triangle & tr,   ����� ������������ 
 *	const Mesh & mesh,     ����� 
 *	int point_i,           ���������� ����� ����� 
 *	int point_j,           ���������� ����� ����� 
 *	int i,                 ����� ������ �������
 *	int j,                 ����� ������� �������
 *	void * user_data       ���� ���� ������� ����� ������
 *
 */

template < typename Functor, typename Data >
void generate_matrix(Matrix & A, const Mesh & m, 
					 Functor integrate_cb, 
					 Data user_data, 
					 bool transpose = false)
{
	using namespace MKE_Private_;
	int rs  = (int)m.inner.size();     // �����������

	Timer t;
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) { // ����� ������
		// �� ���������� ������
		int p = m.inner[i];

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) 
		{
			// �� ������������ � �����
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p2   = m.tr[trk_i].p[i0];
				int j    = m.p2io[p2]; // ����� ���������� �����
				                       // ����� �������
				if (m.ps_flags[p2] == 1) {
					; // �������
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

template < typename Functor, typename Data >
void generate_full_matrix(Matrix & A, const Mesh & m, 
						  Functor integrate_cb, 
						  Data user_data,
						  bool transpose = false)
{
	using namespace MKE_Private_;
	int sz  = (int)m.ps.size();

	Timer t;
//#pragma omp parallel for
	for (int p = 0; p < sz; ++p) {
		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// �� ������������ � �����
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

template < typename Functor, typename Data >
void generate_right_part(double * b, const Mesh & m, 
						 Functor right_part_cb, 
						 Data user_data)
{
	using namespace MKE_Private_;
	int rs  = (int)m.inner.size();     // �����������
	Timer t;

	// WARNING: ���� ���������� ������ ����� ��� ������� ���������,
	// �� �������� ����������� �� ��������, ������� 
	// memset(b, 0) ���� �������� ������ �� ������ generate_right_part !

#pragma omp parallel for
	for (int i = 0; i < rs; ++i)
	{
		// �� ���������� ������
		int p = m.inner[i];
		b[i]  = 0.0;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// �� ������������ � �����
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			// ���� ������� �� �������� �� �������, �� ���� ������ ��������� callback
			// ��� ���������� ���������� �������� � ������ �����?
			
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

template < typename Functor, typename Data >
void generate_full_right_part(double * b, const Mesh & m, 
							  Functor right_part_cb, 
							  Data user_data)
{
	using namespace MKE_Private_;
	int sz  = (int)m.ps.size();     // �����������

	Timer t;

	// WARNING: ���� ���������� ������ ����� ��� ������� ���������,
	// �� �������� ����������� �� ��������, ������� 
	// memset(b, 0) ���� �������� ������ �� ������ generate_right_part !

#pragma omp parallel for
	for (int p = 0; p < sz; ++p)
	{
		b[p]  = 0.0;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// �� ������������ � �����
			int trk_i = m.adj[p][tk];
			const Triangle & trk    = m.tr[trk_i];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_i           = trk.elem1(p);
			
			// ���� ������� �� �������� �� �������, �� ���� ������ ��������� callback
			// ��� ���������� ���������� �������� � ������ �����?
			
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
 * ���������� ������� ��� ���������� ������� ������� � ������ �����
 * inner.size() x outer.size()
 * � cb ���������� phi_j ��� j ����� �������
 * phi_i, ��� i ���������� �����
 */
template < typename Functor, typename Data >
void generate_boundary_matrix(Matrix & A, const Mesh & m, 
							  Functor right_part_cb, 
							  Data user_data,
							  bool transpose = false)
{
	using namespace MKE_Private_;
	int os = (int)m.outer.size(); // ������ �������
	for (int j = 0; j < os; ++j) {
		// �� ������� ������
		int p2 = m.outer[j];

		for (uint tk = 0; tk < m.adj[p2].size(); ++tk) {
			// �� ������������� � �����
			int trk_j = m.adj[p2][tk];
			const Triangle & trk    = m.tr[trk_j];
			const std::vector < Polynom > & phik = trk.elem1();
			const Polynom & phi_j           = trk.elem1(p2);

			for (uint i0 = 0; i0 < phik.size(); ++i0) {
				int p    = m.tr[trk_j].p[i0];

				if (m.ps_flags[p] == 1) {
					;
				} else {
					// p - ���������� �����
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
	using namespace MKE_Private_;

	int sz  = (int)m.ps.size(); // �����������

//#pragma omp parallel for
	for (int i = 0; i < sz; ++i)
	{
		// �� ���� ������
		int p = i;
		ans[i] = 0.0;

		for (uint tk = 0; tk < m.adj[p].size(); ++tk) {
			// �� ������������� � �����
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

/* �������� ��������� ������������ ���� ������� */
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

/* �������� ����� */
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

/* �������� ���������� */
template < typename Functor, typename Data>
double dist(const double * u, const double * v, const Mesh & m, 
				Functor cb, Data user_data)
{
	int sz  = (int)m.ps.size(); // �����������
	std::vector < double > diff(sz);
	vec_diff(&diff[0], u, v, sz);
	return norm(&diff[0], m, cb, user_data);
}

inline double dist(const double * u, const double * v, const Mesh & m)
{
	return dist(u, v, m, generic_scalar_cb, (void*)0);
}

}

