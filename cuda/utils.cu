/*$Id$*/

/**
 * a = b * k
 */
__global__ static void vector_mult_scalar(double * a, const double * b, double k, int n)
{
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        a[i] = b[i] * k;
}

__global__ static void vector_sum1(double * r, const double * a, const double *b, double k1, double k2, int n)
{
        int i = blockIdx.x * blockDim.x + threadIdx.x;;
        r[i] = k1 * a[i] + k2 * b[i];
}

