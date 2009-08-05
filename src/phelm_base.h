#ifndef PHELM_BASE_H
#define PHELM_BASE_H

#include <vector>
#ifdef GPGPU
#include "alloc_cu.h"
#endif

namespace phelm {

void phelm_init();
void phelm_shutdown();

/**
 * @ingroup misc
 * allocator.
 */
#ifndef GPGPU
#define phelm_allocator std::allocator
#else
#define phelm_allocator cuda_allocator
#endif

/**
 * @ingroup misc
 * vector.
 */
typedef std::vector < double, phelm_allocator < double > > vec;

}

#endif /* PHELM_BASE_H */
