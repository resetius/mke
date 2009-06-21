#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "solver.h"

using MKE::Matrix;
using MKE::Mesh;

class SphereJacobian {
	const Mesh & m_;
	Matrix idt_;
	Matrix diff1_;
	Matrix diff2_;
	Matrix diff1_cos_;
	Matrix diff2_cos_;

	Matrix diff1_rp_;
	Matrix diff2_rp_;
	Matrix diff1_cos_rp_;
	Matrix diff2_cos_rp_;

public:
	SphereJacobian(const Mesh & m);

	/**
	 * ������� J(u, v) �� ���������� ������. 
	 * � ������ ������� ������ ������ �������� �� bnd.
	 */
	void calc1(double * Ans, const double * u, const double * v, const double * bnd);

	/**
	 * ����������� � J(u, v) (=-J(u, v))
	 */
	void calc1t(double * Ans, const double * u, const double * v, const double * bnd);

	/**
	 * ������� J(u, v) �� ���������� ������. 
	 * ���������� ������, ���������� ������ ���������� �����
	 */
	void calc2(double * Ans, const double * u, const double * v);

	/**
	 * ����������� � J(u, v) (=-J(u, v))
	 */
	void calc2t(double * Ans, const double * u, const double * v);
};

#endif /* JACOBIAN_H */
