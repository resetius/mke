#include <vector>
#include <math.h>

#include "mke.h"
#include "jacobian.h"
#include "util.h"

using namespace std;

static double 
jacobian(const Polynom & phi_i, const Polynom & phi_j, const Triangle & trk, 
	 const Mesh & m, int i, int j, void * data)
{
	Polynom pt1 = diff(phi_i, 1) * diff(phi_j, 0);
	Polynom pt2 = diff(phi_i, 0) * diff(phi_j, 1);

	Point p = m.ps[i].p[0];
	return (pt1.apply(p.x, p.y) - pt2.apply(p.x, p.y)) / cos(p.x);
}

static double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		void *)
{
	return integrate_cos(phi_i * phi_j, trk, m.ps);
}

static double diff_1_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate(poly, trk, m.ps);
	return r;
}

static double diff_1_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 0) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate_cos(poly, trk, m.ps);
	return r;
}

static double diff_2_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate(poly, trk, m.ps);
	return r;
}

static double diff_2_cos_rp(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int i,
		int j,
		const double * u)
{
	Polynom poly = diff(phi_j, 1) * phi_i;
	double v = (u) ? u[j] : 1;
	double r = v * integrate_cos(poly, trk, m.ps);
	return r;
}

void Jacobian::calc2(double * Ans, const double * u, const double * v)
{
	int rs = (int)m_.inner.size();
	int sz = (int)m_.ps.size();
	int os = (int)m_.outer.size();
#if 0
	vector < double > rp(sz);
	convolution(&rp[0], u, v, m_, (scalar_cb_t)jacobian, 0);
	mke_u2p(Ans, &rp[0], m_);
#endif

#if 1
	vector < double > v_in(rs);
	vector < double > u_in(rs);
	vector < double > v_in_bnd(os);
	vector < double > u_in_bnd(os);

	vector < double > rp(rs);
	vector < double > pt1(rs);
	vector < double > pt2(rs);
	vector < double > tmp(rs);

	mke_u2p(&v_in[0], v, m_);
	mke_u2p(&u_in[0], u, m_);
	mke_proj_bnd(&v_in_bnd[0], &v_in[0], m_);
	mke_proj_bnd(&u_in_bnd[0], &u_in[0], m_);

	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_cos_rp, (void*)u);
	//diff2_cos_.mult_vector(&tmp[0], &u_in[0]);
	//diff2_cos_rp_.mult_vector(&rp[0], &u_in_bnd[0]);
	//vector_sum(&rp[0], &rp[0], &tmp[0], rp.size());
	//idt_.solve(&pt1[0], &rp[0]);

	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_rp, (void*)v);
	idt_.solve(&tmp[0], &rp[0]);
	vector_mult(&pt1[0], &pt1[0], &tmp[0], (int)pt1.size());

	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_1_cos_rp, (void*)u);
	idt_.solve(&pt2[0], &rp[0]);
	generate_right_part(&rp[0], m_, (right_part_cb_t)diff_2_rp, (void*)v);
	idt_.solve(&tmp[0], &rp[0]);
	vector_mult(&pt2[0], &pt2[0], &tmp[0], (int)pt1.size());

	vector_diff(Ans, &pt1[0], &pt2[0], (int)pt1.size());
#endif
}

/**
 * J(u,v)=1/cos(phi) (du/d\la dv/d\phi - du/d\phi dv/d\la)
 */
Jacobian::Jacobian(const Mesh & m): m_(m), 
	idt_((int)m.inner.size()),
	diff1_(m.inner.size()),
	diff2_(m.inner.size()),
	diff1_cos_(m.inner.size()),
	diff2_cos_(m.inner.size()),
	diff1_rp_(m.inner.size()),
	diff2_rp_(m.inner.size()),
	diff1_cos_rp_(m.inner.size()),
	diff2_cos_rp_(m.inner.size())
{
	generate_matrix(idt_, m, id_cb, 0);

	generate_matrix(diff1_, m_, (right_part_cb_t)diff_1_rp, 0);
	generate_matrix(diff2_, m_, (right_part_cb_t)diff_2_rp, 0);
	generate_matrix(diff1_cos_, m_, (right_part_cb_t)diff_1_cos_rp, 0);
	generate_matrix(diff2_cos_, m_, (right_part_cb_t)diff_2_cos_rp, 0);

	generate_boundary_matrix(diff1_rp_, m_, (right_part_cb_t)diff_1_rp, 0);
	generate_boundary_matrix(diff2_rp_, m_, (right_part_cb_t)diff_2_rp, 0);
	generate_boundary_matrix(diff1_cos_rp_, m_, (right_part_cb_t)diff_1_cos_rp, 0);
	generate_boundary_matrix(diff2_cos_rp_, m_, (right_part_cb_t)diff_2_cos_rp, 0);
}

void Jacobian::calc1(double * Ans, const double * u, const double * v, const double * bnd)
{
	vector < double > p1(m_.inner.size());
	calc2(&p1[0], u, v);
	mke_p2u(Ans, &p1[0], bnd, m_);
}

