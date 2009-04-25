#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "solver.h"

class Jacobian {
	const Mesh & m_;
	Matrix idt_;

public:
	Jacobian(const Mesh & m);

	/**
	 * Находит J(u, v) во внутренних точках. 
	 * В точках границы просто кладет значение из bnd.
	 */
	void calc1(double * Ans, const double * u, const double * v, const double * bnd);

	/**
	 * Находит J(u, v) во внутренних точках. 
	 * Возвращает вектор, содержащий ТОЛЬКО внутренние точки
	 */
	void calc2(double * Ans, const double * u, const double * v);
};

#endif /* JACOBIAN_H */
