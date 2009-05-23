
namespace MKE_Private_ 
{
	inline void mat_add(Matrix & A, int i, int j, double a)
	{
		A.add(i, j, a);
	}

	template < typename Container >
	void mat_add(Matrix & A, int i, int j, const Container & c)
	{
		typename c::const_iterator b = c.begin();
		typename c::const_iterator e = c.end();
		typename c::const_iterator it;

		for (it = b; it != e; ++it)
		{
			A.add(it->i, it->j, it->a);
		}
	}

	inline void vec_add(double * b, int i, double a)
	{
		b[i] += a;
	}

	template < typename Container >
	void vec_add(double * b, int i, const Container & c)
	{
		typename c::const_iterator b = c.begin();
		typename c::const_iterator e = c.end();
		typename c::const_iterator it;

		for (it = b; it != e; ++it)
		{
			b[it->i] += it->a;
		}
	}
};
