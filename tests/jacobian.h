#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "solver.h"

class Jacobian {
	const Mesh & m_;
	Matrix idt_;

public:
	Jacobian(const Mesh & m);

	/**
	 * ������� J(u, v) �� ���������� ������. 
	 * � ������ ������� ������ ������ �������� �� bnd.
	 */
	void calc1(double * Ans, const double * u, const double * v, const double * bnd);

	/**
	 * ������� J(u, v) �� ���������� ������. 
	 * ���������� ������, ���������� ������ ���������� �����
	 */
	void calc2(double * Ans, const double * u, const double * v);
};

#endif /* JACOBIAN_H */
