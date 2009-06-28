#ifndef GMRES_H
#define GMRES_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/**
 * @file gmres.h
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @section LICENSE
 *
 * Copyright (c) 2009 Alexey Ozeritsky \n
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without \n
 * modification, are permitted provided that the following conditions \n
 * are met: \n
 * 1. Redistributions of source code must retain the above copyright \n
 *    notice, this list of conditions and the following disclaimer. \n
 * 2. Redistributions in binary form must reproduce the above copyright \n
 *    notice, this list of conditions and the following disclaimer in the \n
 *    documentation and/or other materials provided with the distribution. \n
 * 3. Redistributions in any form must be accompanied by information on \n
 *    how to obtain complete source code for the Phelm software and any \n
 *    accompanying software that uses the Phelm software.  The source code \n
 *    must either be included in the distribution or be available for no \n
 *    more than the cost of distribution plus a nominal fee, and must be \n
 *    freely redistributable under reasonable conditions.  For an \n
 *    executable file, complete source code means the source code for all \n
 *    modules it contains.  It does not include source code for modules or \n
 *    files that typically accompany the major components of the operating \n
 *    system on which the executable file runs. \n
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR \n
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES \n
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. \n
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, \n
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT\n
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE \n,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY \n
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT \n
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF \n
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. \n
 *
 * @section DESCRIPTION
 *
 * Generalized minimal residual method (GMRES) \n
 * http://en.wikipedia.org/wiki/GMRES
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * callback
 * y = Ax
 */
typedef void (*Ax_t)(double * y, const void * A, const double * x, int n);

/**
 * Solve equation Ax=b
 *
 * @param x right part
 * @param A matrix
 * @param b right part
 * @param Ax callback that calculates y=Ax
 * @param n dimension of system
 * @param k_dim Krylov dimension
 * @param max_id maximum numner of iterations
 */
void gmres(double * x, const void * A, const double * b, 
			 Ax_t Ax, int n, int k_dim, int max_it);

#ifdef __cplusplus
}
#endif

#endif /* GMRES_H */

