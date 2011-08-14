/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky
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

namespace SphereLaplace_Private
{
double id_cb (const Polynom & phi_i,
              const Polynom & phi_j,
              const Triangle & trk,
              const Mesh & m,
              int point_i,
              int point_j,
              int, int,
              void *);


double
slaplace_integrate_cb ( const Polynom & phi_i,
                        const Polynom & phi_j,
                        const Triangle & trk,
                        const Mesh & m,
                        int point_i,
                        int point_j,
                        int, int,
                        void * user_data);

double
laplace_bnd1_cb ( const Polynom & phi_i,
                  const Polynom & phi_j,
                  const Triangle & trk,
                  const Mesh & m,
                  int point_i,
                  int point_j,
                  int, int,
                  void * d);

double
laplace_bnd2_cb ( const Polynom & phi_i,
                  const Polynom & phi_j,
                  const Triangle & trk,
                  const Mesh & m,
                  int point_i,
                  int point_j,
                  int, int,
                  void * );

}

template < typename T, typename Matrix >
void SphereLaplace < T, Matrix >::solve (T * Ans,
                                 const T * F, const T * bnd)
{
	using namespace phelm;
	//ïîêà èñïîëüçóåì ïåðâûé ïîðÿäîê
	int sz  = (int) m_.ps.size();
	int ntr = (int) m_.tr.size();
	int rs  = (int) m_.inner.size(); 

	ArrayDevice  b (rs); 
	ArrayDevice  x (rs); 

	Timer full;
#if 0
	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = bnd;
	generate_right_part (&b[0], m_, (right_part_cb_t) (slaplace_right_part_cb), (void*) &d);
#endif

#if 1
//	fprintf(stderr, "stage0\n"); vec_print(stderr, F, m_.size, "%8.3le ");
	u2p (&x[0], F, m_);
//	fprintf(stderr, "stage1\n"); vec_print(stderr, &x[0], rs, "%8.3le ");
	idt_.mult_vector (&b[0], &x[0]);
//	fprintf(stderr, "stage2\n"); vec_print(stderr, &b[0], rs, "%8.3le ");
	if (bnd)
	{
		bnd2_.mult_vector (&x[0], bnd);
//		fprintf(stderr, "stage3\n"); vec_print(stderr, &x[0], rs, "%8.3le ");
		vec_sum (&b[0], &b[0], &x[0], (int) x.size() );
//		fprintf(stderr, "stage4\n"); vec_print(stderr, &b[0], rs, "%8.3le ");
	}
//	vector < double > tmp(m_.outer.size()); // not necessary !
//	proj_bnd(&tmp[0], F, m_);           // not necessary !
//	bnd1_.mult_vector(&x[0], &tmp[0]);      // not necessary !
//	vector_sum(&b[0], &b[0], &x[0], x.size()); // not necessary !
#endif

	phelm::solve (Ans, bnd, &b[0], laplace_, m_);

//	vec_mult_scalar(Ans, Ans, (T)(6371.3 * 6371.3), rs); // mult R^2

#ifdef _DEBUG
	fprintf (stderr, "Total elapsed: %lf \n", full.elapsed() );
#endif
}

template < typename T, typename Matrix >
SphereLaplace < T, Matrix > ::SphereLaplace (const Mesh & m) : m_ (m),
		idt_ ( (int) m.inner.size() ),
		laplace_ ( (int) m.inner.size() ), bnd1_ ( (int) m.inner.size() ), bnd2_ ( (int) m.inner.size() ),
		bnd3_ ( (int) m.inner.size() )
{
  using namespace SphereLaplace_Private;
	generate_matrix (idt_, m, id_cb, (void*) 0);
	generate_matrix (laplace_, m, slaplace_integrate_cb, (void*) 0);
	generate_boundary_matrix (bnd1_, m_, laplace_bnd1_cb, (void*) 0);
	generate_boundary_matrix (bnd2_, m_, laplace_bnd2_cb, (void*) 0);
	generate_boundary_matrix (bnd3_, m_, slaplace_integrate_cb, (void*) 0);
}

template < typename T, typename Matrix >
void SphereLaplace < T, Matrix > ::calc2 (T * Ans, const T * F)
{
	using namespace phelm;

#if 0
	vector < double > rp (m_.inner.size() );

	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = 0;
	generate_right_part (&rp[0], m_, (right_part_cb_t) lp_rp, &d);
	idt_.solve (Ans, &rp[0]);
#endif
	int rs = (int) m_.inner.size();
	int os = (int) m_.outer.size();
	ArrayDevice in (rs);
	ArrayDevice out (rs);
	ArrayDevice tmp (os);
	u2p (&in[0], F, m_);
	proj_bnd (&tmp[0], F, m_);
	laplace_.mult_vector (&out[0], &in[0]);
	bnd3_.mult_vector (&in[0], &tmp[0]);
	vec_sum (&out[0], &out[0], &in[0], (int) in.size() );
	idt_.solve (Ans, &out[0]);

//	vec_mult_scalar(Ans, Ans, (T)(1.0 / 6371.3 / 6371.3), rs); // divide R^2
}

template < typename T, typename Matrix >
void SphereLaplace < T, Matrix > ::calc1 (T * Ans, const T * F, const T * bnd)
{
	ArrayDevice p1 (m_.inner.size() );

	calc2 (&p1[0], F);
#if 0
	vector < double > rp (m_.inner.size() );

	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = 0;
	generate_right_part (&rp[0], m_, (right_part_cb_t) lp_rp, &d);
	idt_.solve (&p1[0], &rp[0]);
#endif
	p2u (Ans, &p1[0], bnd, m_);
}
