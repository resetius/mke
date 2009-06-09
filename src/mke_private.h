
namespace MKE
{

namespace MKE_Private_ 
{
	inline void mat_add(Matrix & A, int i, int j, double a)
	{
		A.add(i, j, a);
	}

	template < typename Container >
	void mat_add(Matrix & A, int i, int j, const Container & c)
	{
		typename Container::const_iterator b = c.begin();
		typename Container::const_iterator e = c.end();
		typename Container::const_iterator it;

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
	void vec_add(double * b1, int i, const Container & c)
	{
		typename Container::const_iterator b = c.begin();
		typename Container::const_iterator e = c.end();
		typename Container::const_iterator it;

		for (it = b; it != e; ++it)
		{
			b1[it->i] += it->a;
		}
	}
};
}

