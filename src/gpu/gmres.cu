#include "linal_cuda.h"

namespace phelm {

template < typename T >
__global__ void
gmres_vec_sum_(T * x, const T * q, const T * y, int j, int n)
{
	int threads = gridDim.x  * blockDim.x;
	int k = blockDim.x * blockIdx.x + threadIdx.x;
	int i;
	for (;k < n; k += threads) {
		for (i = 0; i <= j; ++i) {
			x[k] += y[i] * q[i * n + k];
		}
	}
}

__host__ void gmres_vec_sum2(float * x, const float * q, const float * y, int j, int hz, int n)
{
	SPLAY(n);
	float * y1;
	cudaMalloc((void**)&y1, hz * sizeof(float));
	cudaMemcpy(y1, y, hz * sizeof(float), cudaMemcpyHostToDevice);
	gmres_vec_sum_ <<< blocks, threads >>> (x, q, y1, j, n);
	cudaFree(y1);
}

__host__ void gmres_vec_sum2(double * x, const double * q, const double * y, int j, int hz, int n)
{
	SPLAY(n);
	double * y1;
	cudaMalloc((void**)&y1, hz * sizeof(double));
	cudaMemcpy(y1, y, hz * sizeof(double), cudaMemcpyHostToDevice);
	gmres_vec_sum_ <<< blocks, threads >>> (x, q, y1, j, n);
	cudaFree(y1);
}

}
