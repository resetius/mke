
#include <umfpack.h>

/**
 * UMFPACK Solver class.
 */
template < typename T, typename MultStore >
class UmfPackSolver: public SparseSolver < T, MultStore, StoreCSR < T, std::allocator > >
{
	double Control_ [UMFPACK_CONTROL];
	double Info_ [UMFPACK_INFO];
	void *Symbolic_, *Numeric_ ;

	typedef StoreCSR < T, std::allocator > CSR;
	typedef SparseSolver < T, MultStore, CSR > base;

public:
	typedef T data_type;
	typedef UmfPackSolver < T, MultStore > my_type;

	UmfPackSolver (int n) : base (n), Symbolic_ (0), Numeric_ (0)
	{
		umfpack_di_defaults (Control_);
	}

	~UmfPackSolver()
	{
		umfpack_di_free_symbolic (&Symbolic_);
		umfpack_di_free_numeric (&Numeric_);
	}

	/**
	 * Solve equation Ax = b.
	 * That function uses UMFPACK
	 * @param x - answer
	 * @param b - right part
	 */
	void solve (T * x, const T * b)
	{
		CSR & invert = base::store_.invert;
		base::prepare();

		int status;

		if (Symbolic_ == 0)
		{
			status = umfpack_di_symbolic (invert.n_,
			                              invert.n_,
			                              &invert.Ap_[0], &invert.Ai_[0],
			                              &invert.Ax_[0],
			                              &Symbolic_, Control_, Info_);
			assert (status == UMFPACK_OK);
		}

		if (Numeric_ == 0)
		{
			status = umfpack_di_numeric (&invert.Ap_[0], &invert.Ai_[0],
			                             &invert.Ax_[0],
			                             Symbolic_, &Numeric_, Control_, Info_) ;
			assert (status == UMFPACK_OK);
		}

		status = umfpack_di_solve (UMFPACK_At, &invert.Ap_[0],
		                           &invert.Ai_[0],
		                           &invert.Ax_[0], x, b, Numeric_, Control_, Info_);
		assert (status == UMFPACK_OK);
	}
};

template < typename MultStore >
class UmfPackSolver < float, MultStore >: public SparseSolver < float, MultStore, StoreCSR < float, std::allocator > >
{
	double Control_ [UMFPACK_CONTROL];
	double Info_ [UMFPACK_INFO];
	void *Symbolic_, *Numeric_ ;

	typedef StoreCSR < float, std::allocator > CSR;
	typedef SparseSolver < float, MultStore, CSR > base;

	std::vector < double > Ax_;

public:
	typedef float data_type;
	typedef UmfPackSolver < float, MultStore > my_type;

	UmfPackSolver (int n) : base (n), Symbolic_ (0), Numeric_ (0)
	{
		umfpack_di_defaults (Control_);
	}

	~UmfPackSolver()
	{
		umfpack_di_free_symbolic (&Symbolic_);
		umfpack_di_free_numeric (&Numeric_);
	}

	/**
	* Solve equation Ax = b.
	* That function uses UMFPACK
	* @param x - answer
	* @param b - right part
	*/
	void solve (float * x, const float * b)
	{
		CSR & invert = base::store_.invert;
		base::prepare();

		if (Ax_.empty() )
		{
			Ax_.resize (invert.Ax_.size() );

			for (int i = 0; i < (int) Ax_.size(); ++i)
			{
				Ax_[i] = (double) invert.Ax_[i];
			}
		}

		int status = 0;

		if (Symbolic_ == 0)
		{
			status = umfpack_di_symbolic (invert.n_, invert.n_,
			                              &invert.Ap_[0], &invert.Ai_[0], &Ax_[0],
			                              &Symbolic_, Control_, Info_);
			assert (status == UMFPACK_OK);
		}

		if (Numeric_ == 0)
		{
			status = umfpack_di_numeric (&invert.Ap_[0], &invert.Ai_[0],
			                             &Ax_[0],
			                             Symbolic_, &Numeric_, Control_, Info_) ;
			assert (status == UMFPACK_OK);
		}

		std::vector < double > x1 (invert.n_);
		std::vector < double > b1 (invert.n_);

		for (int i = 0; i < invert.n_; ++i)
		{
			b1[i] = (double) b[i];
		}

		status = umfpack_di_solve (UMFPACK_At, &invert.Ap_[0], &invert.Ai_[0],
		                           &Ax_[0], &x1[0], &b1[0], Numeric_, Control_, Info_);
		assert (status == UMFPACK_OK);

		for (int i = 0; i < invert.n_; ++i)
		{
			x[i] = (float) x1[i];
		}
	}
};
