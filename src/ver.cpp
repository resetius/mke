#include <string>
#include <set>
#include <stdlib.h>
#include "ver.h"

class StringsStorage
{
	std::set < std::string > storage_;

	StringsStorage() {}
	~StringsStorage() {};

	StringsStorage (const StringsStorage &);
	StringsStorage & operator = (const StringsStorage &);

public:
	static StringsStorage & instance()
	{
		static StringsStorage s;
		return s;
	}

	void add (const std::string & str)
	{
		storage_.insert (str);
	}
};

extern "C"
const char * add_version_string (const char * const str)
{
	StringsStorage::instance().add (str);
	return str;
}

static const char * const license = add_version_string (
                                        " \
/** \n \
 * Copyright (c) 2009\n \
 *      Alexey Ozeritsky.  All rights reserved.\n \
 *\n \
 * Redistribution and use in source and binary forms, with or without\n \
 * modification, are permitted provided that the following conditions\n \
 * are met:\n \
 * 1. Redistributions of source code must retain the above copyright\n \
 *    notice, this list of conditions and the following disclaimer.\n \
 * 2. Redistributions in binary form must reproduce the above copyright\n \
 *    notice, this list of conditions and the following disclaimer in the\n \
 *    documentation and/or other materials provided with the distribution.\n \
 * 3. Redistributions in any form must be accompanied by information on\n \
 *    how to obtain complete source code for the MKE software and any\n \
 *    accompanying software that uses the MKE software.  The source code\n \
 *    must either be included in the distribution or be available for no\n \
 *    more than the cost of distribution plus a nominal fee, and must be\n \
 *    freely redistributable under reasonable conditions.  For an\n \
 *    executable file, complete source code means the source code for all\n \
 *    modules it contains.  It does not include source code for modules or\n \
 *    files that typically accompany the major components of the operating\n \
 *    system on which the executable file runs.\n \
 *\n \
 * THIS SOFTWARE IS PROVIDED BY ORACLE CORPORATION ``AS IS'' AND ANY EXPRESS\n \
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n \
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR\n \
 * NON-INFRINGEMENT, ARE DISCLAIMED.  IN NO EVENT SHALL ORACLE CORPORATION\n \
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR\n \
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF\n \
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS\n \
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN\n \
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)\n \
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF\n \
 * THE POSSIBILITY OF SUCH DAMAGE.\n \
 */\n \
");

