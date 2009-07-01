#ifndef JACOBIAN_H
#define JACOBIAN_H
/* -*- charset: utf-8 -*- */
/* $Id$ */

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
 * 
 * @section DESCRIPTION
 * The spherical jacobian.
 */

#include "solver.h"

using phelm::Matrix;
using phelm::Mesh;

/**
 * @ingroup aux
 * @{
 */

/**
 * Calculates the jacobian on a sphere.
 \f[
 J(u,v) = \frac{1}{cos(\varphi)} (\frac{\partial u}{\partial \lambda}\frac{\partial v}{\partial \varphi} 
 - \frac{\partial u}{\partial \varphi}\frac{\partial v}{\partial \lambda})
 \f] 
 */
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

	Matrix diff1_t_;
	Matrix diff2_t_;
	Matrix diff1_cos_t_;
	Matrix diff2_cos_t_;

	Matrix diff1_rp_t_;
	Matrix diff2_rp_t_;
	Matrix diff1_cos_rp_t_;
	Matrix diff2_cos_rp_t_;

public:
	/**
	 * constructor.
	 */
	SphereJacobian(const Mesh & m);

	/**
	 * Находит J(u, v) во внутренних точках. 
	 * В точках границы просто кладет значение из bnd.
	 */
	void calc1(double * Ans, const double * u, const double * v, const double * bnd);

	/**
	 * сопряженный к J(u, v) (=-J(u, v))
	 */
	void calc1t(double * Ans, const double * u, const double * v, const double * bnd);

	/**
	 * Находит J(u, v) во внутренних точках. 
	 * Возвращает вектор, содержащий ТОЛЬКО внутренние точки
	 */
	void calc2(double * Ans, const double * u, const double * v);

	/**
	 * сопряженный к J(u, v) (=-J(u, v))
	 */
	void calc2t(double * Ans, const double * u, const double * v);
};

/** @} */

#endif /* JACOBIAN_H */

