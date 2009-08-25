#ifndef TEXTURE_H
#define TEXTURE_H

/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Redistributions in any form must be accompanied by information on
 *    how to obtain complete source code for the Phelm software and any
 *    accompanying software that uses the Phelm software.  The source code
 *    must either be included in the distribution or be available for no
 *    more than the cost of distribution plus a nominal fee, and must be
 *    freely redistributable under reasonable conditions.  For an
 *    executable file, complete source code means the source code for all
 *    modules it contains.  It does not include source code for modules or
 *    files that typically accompany the major components of the operating
 *    system on which the executable file runs.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdexcept>

#define register_texture(type, name) \
	texture < type > tex_##name; \
	struct reader_tex_##name \
	{ \
		size_t off; \
		type * ptr; \
		int size; \
		reader_tex_##name (const type * p, int n): ptr((type*)p), size(n)  \
		{ \
			bind(); \
		} \
		void bind() \
		{ \
			if (cudaBindTexture(&off, tex_##name, ptr, size * sizeof(ptr[0])) != cudaSuccess) { \
				throw std::runtime_error("cannot bind texture \n"); \
			} \
			off /= sizeof(ptr[0]); \
		} \
		void unbind() \
		{ \
			if (cudaUnbindTexture(tex_##name) !=  cudaSuccess) \
			{ \
				throw std::runtime_error("cannot unbind texture \n"); \
			} \
		} \
		__device__ type get (size_t i) \
			{ return tex1Dfetch(tex_##name, off + i); } \
	};

#define texture_reader(name) \
	reader_tex_##name

template < typename T >
struct simple_reader
{
	const T * ptr;
	size_t off;
	simple_reader(const T * p): ptr(p), off(0) {}
	__device__ T get(size_t i) { return ptr[i]; }
};

#define WORD_ALIGN     (64)   /* alignment for 32-bit word */
#define LONG_ALIGN     (128)  /* alignment for 64-bit long */
#define _1DBUF_ALIGN   (256)  /* alignment for 1D buffer */
#define MAX_1DBUF_SIZE ((1<<27)-_1DBUF_ALIGN)

#endif /* TEXTURE_H */
