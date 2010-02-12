#include <slu_ddefs.h>

extern "C" void
sgstrf (superlu_options_t *options, SuperMatrix *A,
        int relax, int panel_size, int *etree, void *work, int lwork,
        int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
        SuperLUStat_t *stat, int *info);

extern "C" void
psgstrf (superlu_options_t *options, SuperMatrix *A,
         int relax, int panel_size, int *etree, void *work, int lwork,
         int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
         SuperLUStat_t *stat, int *info);

extern "C" void
pdgstrf (superlu_options_t *options, SuperMatrix *A,
         int relax, int panel_size, int *etree, void *work, int lwork,
         int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
         SuperLUStat_t *stat, int *info);

extern "C" void
sgstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
        int *perm_c, int *perm_r, SuperMatrix *B,
        SuperLUStat_t *stat, int *info);

extern "C" void
sCreate_CompCol_Matrix (SuperMatrix *A, int m, int n, int nnz,
                        float *nzval, int *rowind, int *colptr,
                        Stype_t stype, Dtype_t dtype, Mtype_t mtype);

extern "C" void
sCreate_CompRow_Matrix (SuperMatrix *A, int m, int n, int nnz,
                        float *nzval, int *rowind, int *colptr,
                        Stype_t stype, Dtype_t dtype, Mtype_t mtype);

extern "C" void
sCreate_Dense_Matrix (SuperMatrix *X, int m, int n, float *x, int ldx,
                      Stype_t stype, Dtype_t dtype, Mtype_t mtype);

#ifdef SUPERLU_MT
#define MT_CALL(func) \
	p##func
#else
#define MT_CALL(func) \
	func
#endif


template < typename T >
struct SuperLu
{
	static void gstrf (superlu_options_t *options, SuperMatrix *A,
	                   int relax, int panel_size, int *etree, void *work, int lwork,
	                   int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
	                   SuperLUStat_t *stat, int *info)
	{
		assert (0);
	}

	static void gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
	                   int *perm_c, int *perm_r, SuperMatrix *B,
	                   SuperLUStat_t *stat, int *info)
	{
		assert (0);
	}

	static void
	Create_CompCol_Matrix (SuperMatrix *A, int m, int n, int nnz,
	                       T *nzval, int *rowind, int *colptr,
	                       Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
		assert (0);
	}

	static void
	Create_CompRow_Matrix (SuperMatrix *A, int m, int n, int nnz,
	                       T *nzval, int *rowind, int *colptr,
	                       Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
		assert (0);
	}

	static void
	Create_Dense_Matrix (SuperMatrix *X, int m, int n, T *x, int ldx,
	                     Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
		assert (0);
	}
};

template < >
struct SuperLu < double >
{
	static void gstrf (superlu_options_t *options, SuperMatrix *A,
	                   int relax, int panel_size, int *etree, void *work, int lwork,
	                   int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
	                   SuperLUStat_t *stat, int *info)
	{
		MT_CALL (dgstrf) (options, A, relax, panel_size, etree, work, lwork, perm_c,
		                  perm_r, L, U, stat, info);
	}

	static void gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
	                   int *perm_c, int *perm_r, SuperMatrix *B,
	                   SuperLUStat_t *stat, int *info)
	{
		dgstrs (trans, L, U, perm_c, perm_r, B, stat, info);
	}

	static void
	Create_CompCol_Matrix (SuperMatrix *A, int m, int n, int nnz,
	                       double *nzval, int *rowind, int *colptr,
	                       Stype_t stype, Mtype_t mtype)
	{
		dCreate_CompCol_Matrix (A, m, n, nnz, nzval, rowind, colptr, stype, SLU_D, mtype);
	}

	static void
	Create_CompRow_Matrix (SuperMatrix *A, int m, int n, int nnz,
	                       double *nzval, int *rowind, int *colptr,
	                       Stype_t stype, Mtype_t mtype)
	{
		dCreate_CompRow_Matrix (A, m, n, nnz, nzval, rowind, colptr, stype, SLU_D, mtype);
	}

	static void
	Create_Dense_Matrix (SuperMatrix *X, int m, int n, double *x, int ldx,
	                     Stype_t stype, Mtype_t mtype)
	{
		dCreate_Dense_Matrix (X, m, n, x, ldx, stype, SLU_D, mtype);
	}
};

template < >
struct SuperLu < float >
{
	static void gstrf (superlu_options_t *options, SuperMatrix *A,
	                   int relax, int panel_size, int *etree, void *work, int lwork,
	                   int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
	                   SuperLUStat_t *stat, int *info)
	{
		MT_CALL (sgstrf) (options, A, relax, panel_size, etree, work, lwork, perm_c,
		                  perm_r, L, U, stat, info);
	}

	static void gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U,
	                   int *perm_c, int *perm_r, SuperMatrix *B,
	                   SuperLUStat_t *stat, int *info)
	{
		sgstrs (trans, L, U, perm_c, perm_r, B, stat, info);
	}

	static void
	Create_CompCol_Matrix (SuperMatrix *A, int m, int n, int nnz,
	                       float *nzval, int *rowind, int *colptr,
	                       Stype_t stype, Mtype_t mtype)
	{
		sCreate_CompCol_Matrix (A, m, n, nnz, nzval, rowind, colptr, stype, SLU_S, mtype);
	}

	static void
	Create_CompRow_Matrix (SuperMatrix *A, int m, int n, int nnz,
	                       float *nzval, int *rowind, int *colptr,
	                       Stype_t stype, Mtype_t mtype)
	{
		sCreate_CompRow_Matrix (A, m, n, nnz, nzval, rowind, colptr, stype, SLU_S, mtype);
	}

	static void
	Create_Dense_Matrix (SuperMatrix *X, int m, int n, float *x, int ldx,
	                     Stype_t stype, Mtype_t mtype)
	{
		sCreate_Dense_Matrix (X, m, n, x, ldx, stype, SLU_S, mtype);
	}
};

template < typename T, typename MultStore >
class SuperLUSolver: public SparseSolver < T, MultStore, StoreCSR < T, std::allocator > >
{
	typedef StoreCSR < T, std::allocator > CSR;
	typedef SparseSolver < T, MultStore, CSR > base;

	SuperMatrix A_, AC_, L_, U_, B_;

	std::vector < int > perm_c_;
	std::vector < int > perm_r_;
	std::vector < int > etree_;

public:
	typedef T data_type;

	SuperLUSolver (int n) : base (n)
	{
	}

	~SuperLUSolver()
	{
	}

	/**
	 * Solve equation Ax = b.
	 * That function uses SuperLU
	 * @param x - answer
	 * @param b - right part
	 */
	void solve (T * x, const T * b);
};

template < typename T, typename MultStore >
void SuperLUSolver < T, MultStore > ::solve (T * x, const T * b)
{
	CSR & invert = base::store_.invert;

	base::prepare();

	SuperMatrix B;
	SuperLUStat_t stat;
	StatInit (&stat);
	int info;

	if (perm_c_.empty() )
	{
		int panel_size;
		int relax;

		SuperLu < T >::Create_CompCol_Matrix (&A_,
		                                      (int) invert.n_, (int) invert.n_, (int) invert.Ax_.size(),
		                                      &invert.Ax_[0], &invert.Ai_[0], &invert.Ap_[0],
		                                      SLU_NC, SLU_GE);

		//this->print();
		//dPrint_CompCol_Matrix("A", &A_);

		superlu_options_t options;
		set_default_options (&options);

		perm_c_.resize (invert.n_);
		perm_r_.resize (invert.n_);
		etree_.resize (invert.n_);

		get_perm_c (3, &A_, &perm_c_[0]);
		sp_preorder (&options, &A_, &perm_c_[0], &etree_[0], &AC_);

		panel_size = sp_ienv (1);
		relax      = sp_ienv (2);

		SuperLu < T > ::gstrf (&options, &AC_, relax, panel_size, &etree_[0], NULL, 0,
		                       &perm_c_[0], &perm_r_[0], &L_, &U_, &stat, &info);

		assert (info == 0);

		NCPformat *store;
		store = (NCPformat *) L_.Store;
		fprintf (stderr, "  nz: %d, n: %d\n", (int) invert.Ax_.size(), (int) invert.n_);
		fprintf (stderr, "        : %lf%%\n", (double) invert.Ax_.size() / (double) invert.n_ / (double) invert.n_ * 100.0);
		fprintf (stderr, "L nz: %d: %lf%%\n", store->nnz, (double) store->nnz / (double) invert.n_ / (double) invert.n_ * 100.0);
		store = (NCPformat *) U_.Store;
		fprintf (stderr, "U nz: %d: %lf%%\n", store->nnz, (double) store->nnz / (double) invert.n_ / (double) invert.n_ * 100.0);
	}

	memcpy (x, b, invert.n_ * sizeof (T) );
	SuperLu < T >::Create_Dense_Matrix (&B, invert.n_, 1, x, invert.n_, SLU_DN, SLU_GE);

	SuperLu < T > ::gstrs (TRANS, &L_, &U_, &perm_c_[0], &perm_r_[0], &B, &stat, &info);

	assert (info == 0);

	Destroy_SuperMatrix_Store (&B);
	StatFree (&stat);
}
