
namespace SphereLaplace_Private
{
double id_cb(const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int, int,
		void *);


double 
slaplace_integrate_cb( const Polynom & phi_i,
                       const Polynom & phi_j, 
                       const Triangle & trk,
                       const Mesh & m,
                       int point_i,
                       int point_j,
                       int, int,
                       void * user_data);

double
laplace_bnd1_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int, int,
		void * d);

double
laplace_bnd2_cb( const Polynom & phi_i,
		const Polynom & phi_j,
		const Triangle & trk,
		const Mesh & m,
		int point_i,
		int point_j,
		int, int,
		void * );

}

template < typename T >
void SphereLaplace < T >::solve(T * Ans,
			  const T * F, const T * bnd)
{
	//пока используем первый порядок
	int sz  = (int)m_.ps.size();
	int ntr = (int)m_.tr.size();
	int rs  = (int)m_.inner.size();     //размерность

	ArrayDevice < T > b(rs);      // правая часть
	ArrayDevice < T > x(rs);      // ответ

	Timer full;
#if 0
	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = bnd;
	generate_right_part(&b[0], m_, (right_part_cb_t)(slaplace_right_part_cb), (void*)&d);
#endif

#if 1
	u2p(&x[0], F, m_);
	idt_.mult_vector(&b[0], &x[0]);
	if (bnd) {
		bnd2_.mult_vector(&x[0], bnd);
		vec_sum(&b[0], &b[0], &x[0], (int)x.size());
	}
//	vector < double > tmp(m_.outer.size()); // not necessary !
//	proj_bnd(&tmp[0], F, m_);           // not necessary !
//	bnd1_.mult_vector(&x[0], &tmp[0]);      // not necessary !
//	vector_sum(&b[0], &b[0], &x[0], x.size()); // not necessary !
#endif

	phelm::solve(Ans, bnd, &b[0], laplace_, m_);
#ifdef _DEBUG
	fprintf(stderr, "Total elapsed: %lf \n", full.elapsed()); 
#endif
}

template < typename T >
SphereLaplace < T > ::SphereLaplace(const Mesh & m): m_(m), 
	idt_((int)m.inner.size()),
	laplace_((int)m.inner.size()), bnd1_((int)m.inner.size()), bnd2_((int)m.inner.size()),
	bnd3_((int)m.inner.size())
{
	generate_matrix(idt_, m, SphereLaplace_Private::id_cb, (void*)0);
	generate_matrix(laplace_, m, SphereLaplace_Private::slaplace_integrate_cb, (void*)0);
	generate_boundary_matrix(bnd1_, m_, SphereLaplace_Private::laplace_bnd1_cb, (void*)0);
	generate_boundary_matrix(bnd2_, m_, SphereLaplace_Private::laplace_bnd2_cb, (void*)0);
	generate_boundary_matrix(bnd3_, m_, SphereLaplace_Private::slaplace_integrate_cb, (void*)0);
}

template < typename T >
void SphereLaplace < T > ::calc2(T * Ans, const T * F)
{
#if 0
	vector < double > rp(m_.inner.size());

	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = 0;
	generate_right_part(&rp[0], m_, (right_part_cb_t)lp_rp, &d);
	idt_.solve(Ans, &rp[0]);
#endif
	int rs = (int)m_.inner.size();
	int os = (int)m_.outer.size();
	ArrayDevice < T > in(rs);
	ArrayDevice < T > out(rs);
	ArrayDevice < T > tmp(os);
	u2p(&in[0], F, m_);
	proj_bnd(&tmp[0], F, m_);
	laplace_.mult_vector(&out[0], &in[0]);
	bnd3_.mult_vector(&in[0], &tmp[0]);
	vec_sum(&out[0], &out[0], &in[0], (int)in.size());
	idt_.solve(Ans, &out[0]);
}

template < typename T >
void SphereLaplace < T > ::calc1(T * Ans, const T * F, const T * bnd)
{
	ArrayDevice < T > p1(m_.inner.size());

	calc2(&p1[0], F);
#if 0
	vector < double > rp(m_.inner.size());

	slaplace_right_part_cb_data d;
	d.F   = F;
	d.bnd = 0;
	generate_right_part(&rp[0], m_, (right_part_cb_t)lp_rp, &d);
	idt_.solve(&p1[0], &rp[0]);
#endif
	p2u(Ans, &p1[0], bnd, m_);
}

namespace SphereChafe_Private
{

double 
schafe_integrate_cb( const Polynom & phi_i,
                     const Polynom & phi_j, 
                     const Triangle & trk,
                     const Mesh & m,
                     int point_i, int point_j,
					 int, int, 
                     SphereChafeConfig * d);

double f(/*double u,*/ double x, double y, double t, 
				double mu, double sigma);


template < typename T >
struct schafe_right_part_cb_data
{
	const T * F;
	const T * bnd;
	SphereChafeConfig  * d;
};

template < typename T >
double 
schafe_right_part_cb( const Polynom & phi_i,
                      const Polynom & phi_j,
                      const Triangle & trk,
                      const Mesh & m,
                      int point_i, int point_j,
					  int i, int j,
                      schafe_right_part_cb_data < T > * d)
{
	const T * F = d->F;
	double b;

	if (m.ps_flags[point_j] == 1) { // на границе
		int j0       = m.p2io[point_j]; //номер внешней точки
		const T * bnd = d->bnd;
		b = - (double)bnd[j0] * schafe_integrate_cb(phi_i, phi_j, 
			trk, m, point_i, point_j, i, j, d->d);
		//b = 0.0;
	} else {
		b = (double)F[point_j] * integrate_cos(phi_i * phi_j, trk, m.ps);
	}
	return b;
}

}

template < typename T >
SphereChafe < T > ::SphereChafe(const Mesh & m, double tau, double sigma, double mu)
	: SphereChafeConfig(tau, sigma, mu), m_(m), laplace_(m), A_((int)m.inner.size())
{
	/* Матрица левой части */
	/* оператор(u) = u/dt-mu \Delta u/2 + sigma u/2*/

	generate_matrix(A_, m_, SphereChafe_Private::schafe_integrate_cb, this);
}

/*
 * \f$\frac{du}{dt} = \mu \delta u - \sigma u + f (u)\f$
 */
template < typename T >
void SphereChafe < T > ::solve(T * Ans, const T * X0,
						const T * bnd, double t)
{
	int rs  = (int)m_.inner.size();
	int sz  = (int)m_.ps.size();
	ArrayDevice < T > u(rs);
	ArrayDevice < T > p(sz);
	ArrayHost   < T > hp(sz);
	ArrayDevice < T > delta_u(rs);

	ArrayDevice < T > rp(rs);
	ArrayDevice < T > crp(rs);

	// генерируем правую часть
	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)

	u2p(&u[0], X0, m_);
	laplace_.calc2(&delta_u[0], X0);

	// u/dt + mu \Delta u / 2
	vec_sum1(&delta_u[0], &u[0], &delta_u[0], (T)(1.0 / tau_), (T)(mu_ * 0.5), rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vec_sum1(&delta_u[0], &delta_u[0], &u[0], (T)(1.0), (T)(-sigma_ * 0.5), rs);

	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f(u)
#pragma omp parallel for
	for (int i = 0; i < rs; ++i) {
		int point = m_.inner[i];
		double x  = m_.ps[point].x();
		double y  = m_.ps[point].y();

		rp[i] = (T) SphereChafe_Private::f(/*u[i],*/ x, y, t, mu_, sigma_);
	}
	vec_copy_from_host(&crp[0], &rp[0], rs);
	vec_sum(&u[0], &delta_u[0], &crp[0], rs);

	// правую часть на границе не знаем !!!
	p2u(&p[0], &u[0], 0 /*bnd*/, m_);
	vec_copy_from_device(&hp[0], &p[0], sz);

	SphereChafe_Private::schafe_right_part_cb_data < T > data2;
	data2.F   = &p[0];
	data2.bnd = bnd;
	data2.d   = this;
	generate_right_part(&rp[0], m_, 
		SphereChafe_Private::schafe_right_part_cb < T >, &data2);

//	fprintf(stderr, "rp: \n");vector_print(&delta_u[0], rs);
//	fprintf(stderr, "matrix:\n");A_.print();
//	laplace_.print();

	vec_copy_from_host(&crp[0], &rp[0], rs);
	phelm::solve(Ans, bnd, &rp[0], A_, m_);
}
