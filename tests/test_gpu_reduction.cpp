#include <stdlib.h>
#include <cublas.h>
#include "phelm.h"

using namespace phelm;

bool cmp(double a, double b, int n)
{
	return (fabs(a - b) < 1e-15);
}

bool cmp(float a, float b, int n)
{
	return (fabs(a - b)/n < 1e-7f);
}

template < typename T >
void gen_random_vec(T * v, int n)
{
	for (int i = 0; i < n; ++i) {
		v[i] = (T)(2.0 * ((double)rand() / (double)RAND_MAX) - 1.0);
	}
}

template < typename T >
T vec_scalar_host(T * a, T * b, int n)
{
	T sum = 0;
	for (int i = 0; i < n; ++i) {
		sum += a[i] * b[i];
	}
	return sum;
}

template < typename T >
bool test_reduction(int n)
{
	ArrayHost < T > a(n);
	ArrayHost < T > b(n);

	ArrayDevice < T > ca(n);
	ArrayDevice < T > cb(n);

	gen_random_vec(&a[0], n);
	gen_random_vec(&b[0], n);

	vec_copy_from_host(&ca[0], &a[0], n);
	vec_copy_from_host(&cb[0], &b[0], n);

	T host_scalar    = vec_scalar_host(&a[0], &b[0], n);
	T device_scalar  = vec_scalar2(&ca[0], &cb[0], n);
	T device_scalar2 = cublasSdot(n, &ca[0], 1, &cb[0], 1);

	if (!cmp(host_scalar, device_scalar, n) && !cmp(device_scalar, device_scalar2, n)) 
	{
		fprintf(stderr, "%d: %lf != %lf != %lf <=> host != device != cublas \n", n, 
			(double)host_scalar, (double)device_scalar, (double)device_scalar2);
		return false;
	}
	return true;
}

int main()
{
	bool ret = true;
	phelm_init();
	fprintf(stderr, "first\n");
	for (int i = 1025; i < 10000 && ret; ++i) {
		ret &= test_reduction < float > (i);
	}
	fprintf(stderr, "second\n");
	for (int i = 10000; i < 1000000 && ret; i *= 10) {
		ret &= test_reduction < float > (i);
	}
	phelm_shutdown();

	if (!ret) {
		return 1;
	}
	return 0;
}
