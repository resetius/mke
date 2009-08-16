#include "linal_cuda.h"
#include "reduction.h"
#include "texture.h"

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

/*
__global__ void gmres_h_calc_()
{
	
}
*/

/*
__host__ void gmres_h_calc(float * ht, float * h, int hz, int n)
{
	int maxThreads = 256;
	int maxBlocks  = 64;
	int threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
	int blocks  = (n + (threads * 2 - 1)) / (threads * 2);
	blocks  = min(maxBlocks, blocks);

	{
		texture_reader(texA) AR(a, n);
		texture_reader(texB) BR(b, n);
		Multiplier < T, texture_reader(texA), texture_reader(texB) > m(AR, BR);

		reduce6 (threads, blocks, m, v2, n);
	}
}
*/

}
