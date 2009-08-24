#ifndef MKE_H
#define MKE_H
/* -*- charset: utf-8 -*- */
/*$Id: phelm.h,v 8a841ed5b1e5 2009/08/24 14:07:14 aozeritsky $*/

/**
 * @file 
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision: 8a841ed5b1e5 $
 *
 * @page License
 * @section LICENSE
 *
 * @verbatim
  Copyright (c) 2009 Alexey Ozeritsky
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  3. Redistributions in any form must be accompanied by information on
     how to obtain complete source code for the Phelm software and any
     accompanying software that uses the Phelm software.  The source code
     must either be included in the distribution or be available for no
     more than the cost of distribution plus a nominal fee, and must be
     freely redistributable under reasonable conditions.  For an
     executable file, complete source code means the source code for all
     modules it contains.  It does not include source code for modules or
     files that typically accompany the major components of the operating
     system on which the executable file runs.
 
  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  @endverbatim
 * If you want to use Phelm in closed-source software 
 * contact me by sending email to aozeritsky@gmail.com.
 *
 * @mainpage Phelm Documentation
 * @section into_sec Introduction
 * This library is meant to solve partial differential equations.
 * It implements the Finite Element Method.
 * The library is able to solve partial differential equations 
 * on two-dimensional smooth manifolds. 
 *
 * To solve a problem you need the following:
 * - triangulate the manifold;
 * - partition the manifold into subdomains;
 * - specify a local coordinate system in each subdomain;
 * - specify a function of surface integral in the local coordinates.
 *
 * The library already contains triangulation builders for spherical surfaces 
 * and flat rectangular domains.
 *
 * The sample result of mesh builder is shown below: 
 * @image html zones.png
 * In that case the sphere is built on 4 parts.
 * 
 * @section ex Usage examples
 * These examples are to give you some tips on Phelm features.
    -# @ref test_laplace.cpp "Laplace equation on a flat domain"
  \f{eqnarray*}
  \Delta u &=& f(x, y) \\
  u|_{\partial\Omega}&=&u_0
  \f}
    -# @ref test_laplace.cpp "Laplace equation on a sphere"
  \f{eqnarray*}
  \Delta \psi &=& f(\varphi, \lambda) \\
  \Delta \psi &=& \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi\\
  \psi|_{\partial\Omega}&=&\psi_0 \\
  \f}
    -# @ref test_system_laplace.cpp "Double Laplace equations on a flat domain"
  \f{eqnarray*}
  \Delta u + v &=& f(x, y)\\
  u + \Delta v &=& g(x, y)\\
  u|_{\partial\Omega}&=&u_0\\
  v|_{\partial\Omega}&=&v_0\\
  \f}
    -# @ref test_chafe.cpp "Chafe-Infante equation on a flat domain"
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{\partial\Omega}&=&a \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
    -# @ref test_chafe.cpp "Chafe-Infante equation on a sphere"
    -# @ref test_barvortex.cpp "The Barotropic vorticity equation on a sphere"
  \f{eqnarray*}
  \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi &=& f(\varphi, \lambda) \\
	\psi|_{t=0}=\psi_0
  \f}
    -# @ref test_baroclin.cpp  "The two-dimensional baroclinic atmosphere equations  on a sphere"
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	-\alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,\\
	u_1|_{t=0}&=&u_{10}\\
	u_2|_{t=0}&=&u_{20}\\
 \f}
 *
 * @page Build
 * @section build_sec Build
 * Phelm uses GMRES method to solve linear equations. If you do not like
 * that, install UMFPACK or place UMFPACK sources in 
 * /path-to-phelm/contrib/umfpack
 * 
 * @subsection Unix-like
 * @verbatim
  mkdir build-directory
  cd build-directory
  cmake -DCMAKE_BUILD_TYPE=Debug path-to-sources   # for Debug build
  cmake -DCMAKE_BUILD_TYPE=Release path-to-sources # for Release build
  make
  @endverbatim
 * @subsection Windows
 * @verbatim
  mkdir build-directory
  cd build-directory
  cmake -G "Visual Studio 2009" #place your version of Visual Studio here
  @endverbatim
 * @page Thanks
 * Thanks to
 *  - Andrey Kornev for consultations, <br/>
 *  - Andrey Ivanchikov for initial idea and consultations. 
 */

/**
 * 
  @example test_laplace.cpp
  Laplace equation on a flat domain
  \f{eqnarray*}
  \Delta u &=& f(x, y) \\
  u|_{\partial\Omega}&=&u_0
  \f}
  @example test_slaplace.cpp
  Laplace equation on a sphere
  \f{eqnarray*}
  \Delta \psi &=& f(\varphi, \lambda) \\
  \Delta \psi &=& \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi\\
  \psi|_{\partial\Omega}&=&\psi_0 \\
  \f}
  @example test_system_laplace.cpp
  Double Laplace equations on a flat domain
  \f{eqnarray*}
  \Delta u + v &=& f(x, y)\\
  u + \Delta v &=& g(x, y)\\
  u|_{\partial\Omega}&=&u_0\\
  v|_{\partial\Omega}&=&v_0\\
  \f}
  @example test_chafe.cpp
  Chafe-Infante equation on a flat domain
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{\partial\Omega}&=&a \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
  @example test_schafe.cpp
  Chafe-Infante equation on a sphere
  @example test_barvortex.cpp
  the Barotropic vorticity equation
  \f[
  \frac{\partial \Delta \varphi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
  \f]
  @example test_baroclin.cpp
  The two-dimensional baroclinic atmosphere equations
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	-\alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
  */

#include <stdio.h>
#include <vector>

#include "base.h"
#include "polynom.h"
#include "solver.h"
#include "mesh.h"
#include "generators.h"

#define PHELM_VERSION    2
#define PHELM_PATCHLEVEL 0
#define PHELM_SUBLEVEL   0

#endif /* MKE_H */

