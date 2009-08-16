
#include "shmem.h"

inline unsigned int nextPow2( unsigned int x ) 
{
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

template <unsigned int blockSize, typename IN, typename T>
__global__ void
reduce5_(IN v, T * g_odata, unsigned int n)
{
    SharedMemory<T> smem;
    T *sdata = smem.getPointer();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;

    sdata[tid] = (i < n) ? v.get(i) : 0;
    if (i + blockSize < n) 
        sdata[tid] += v.get(i+blockSize);

    __syncthreads();

    // do reduction in shared mem
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    
#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; EMUSYNC; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; EMUSYNC; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; EMUSYNC; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; EMUSYNC; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; EMUSYNC; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; EMUSYNC; }
    }
    
    // write result for this block to global mem 
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template < unsigned int blockSize, bool nIsPow2, typename IN, typename T >
__global__ void
reduce6_(IN v, T * g_odata, unsigned int n)
{
	SharedMemory<T> smem;
	T *sdata = smem.getPointer();

	// perform first level of reduction,
	// reading from global memory, writing to shared memory
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
	unsigned int gridSize = blockSize*2*gridDim.x;
	sdata[tid] = 0;

	// we reduce multiple elements per thread.  The number is determined by the 
	// number of active thread blocks (via gridDim).  More blocks will result
	// in a larger gridSize and therefore fewer elements per thread
	while (i < n)
	{
		sdata[tid] += v.get(i);
		// ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
		if (nIsPow2 || i + blockSize < n) 
			sdata[tid] += v.get(i+blockSize);
		i += gridSize;
	}
	__syncthreads();

	// do reduction in shared mem
	if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }

#ifndef __DEVICE_EMULATION__
	if (tid < 32)
#endif
	{
		if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; EMUSYNC; }
		if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; EMUSYNC; }
		if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; EMUSYNC; }
		if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; EMUSYNC; }
		if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; EMUSYNC; }
		if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; EMUSYNC; }
	}

	// write result for this block to global mem 
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

inline bool isPow2(unsigned int x)
{
    return ((x&(x-1))==0);
}

template < typename IN, typename T >
__host__ void
reduce6(int threads, int blocks, IN in, T * out, unsigned int n)
{
	int smemSize = threads * sizeof(T);
	if (isPow2(n))
	{
		switch (threads)
		{
		case 512:
			reduce6_<512, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 256:
			reduce6_<256, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 128:
			reduce6_<128, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 64:
			reduce6_< 64, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 32:
			reduce6_< 32, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 16:
			reduce6_< 16, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  8:
			reduce6_<  8, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  4:
			reduce6_<  4, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  2:
			reduce6_<  2, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  1:
			reduce6_<  1, true><<< blocks, threads, smemSize >>>(in, out, n); break;
		}
	}
	else
	{
		switch (threads)
		{
		case 512:
			reduce6_<512, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 256:
			reduce6_<256, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 128:
			reduce6_<128, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 64:
			reduce6_< 64, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 32:
			reduce6_< 32, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case 16:
			reduce6_< 16, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  8:
			reduce6_<  8, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  4:
			reduce6_<  4, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  2:
			reduce6_<  2, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		case  1:
			reduce6_<  1, false><<< blocks, threads, smemSize >>>(in, out, n); break;
		}
	}
}

template < typename IN, typename T >
__host__ void
reduce5(int threads, int blocks, IN in, T * out, unsigned int n)
{
	int smemSize = threads * sizeof(T);
	switch (threads)
	{
	case 512:
		reduce5_<512><<< blocks, threads, smemSize >>>(in, out, n); break;
	case 256:
		reduce5_<256><<< blocks, threads, smemSize >>>(in, out, n); break;
	case 128:
		reduce5_<128><<< blocks, threads, smemSize >>>(in, out, n); break;
	case 64:
		reduce5_< 64><<< blocks, threads, smemSize >>>(in, out, n); break;
	case 32:
		reduce5_< 32><<< blocks, threads, smemSize >>>(in, out, n); break;
	case 16:
		reduce5_< 16><<< blocks, threads, smemSize >>>(in, out, n); break;
	case  8:
		reduce5_<  8><<< blocks, threads, smemSize >>>(in, out, n); break;
	case  4:
		reduce5_<  4><<< blocks, threads, smemSize >>>(in, out, n); break;
	case  2:
		reduce5_<  2><<< blocks, threads, smemSize >>>(in, out, n); break;
	case  1:
		reduce5_<  1><<< blocks, threads, smemSize >>>(in, out, n); break;
	}
}
