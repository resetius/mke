#ifndef ALLOC_CU_H
#define ALLOC_CU_H

#include <stdexcept>
#include <cublas.h>
#include <cuda_runtime_api.h>

namespace phelm {

template <class T> class cuda_allocator
{
public:
	typedef size_t    size_type;
	typedef ptrdiff_t difference_type;
	typedef T*        pointer;
	typedef const T*  const_pointer;
	typedef T&        reference;
	typedef const T&  const_reference;
	typedef T         value_type;

	template <class U> struct rebind
	{
		typedef cuda_allocator<U> other;
	};

	cuda_allocator() throw() {};
	cuda_allocator (const cuda_allocator&) throw() {};

	template <class U> cuda_allocator (const cuda_allocator<U>&) throw() {};
	~cuda_allocator() throw() {};

	pointer address (reference x) const { return &x; }
	const_pointer address (const_reference x) const { return &x; }

	pointer allocate (size_type size,
	                  const void * hint = 0)
	{
		pointer ret = 0;
		if (cublasAlloc((int)size, sizeof(T), (void**)&ret) != CUBLAS_STATUS_SUCCESS)
		{
			throw std::bad_alloc();
			return 0;
		}

		//cudaMemset(ret, 0, size * sizeof(T));
		return ret;
	}

	void deallocate (pointer p, size_type n)
	{
		cublasFree(p);
	}

	size_type max_size() const throw()
	{
		return (size_type)1 << (sizeof(size_type) * 8 - 1);
	}

	void construct (pointer p, const T& val)
	{
	}

	void destroy (pointer p)
	{
	}
};

template <class T1, class T2>
bool operator== (const cuda_allocator<T1>&, const cuda_allocator<T2>&)
throw()
{
	return true;
}

template <class T1, class T2>
bool operator!= (const cuda_allocator<T1>&, const cuda_allocator<T2>&)
throw()
{
	return false;
}

}

#endif /* ALLOC_CU_H */

