#include "netcdf_cmd.hpp"
#include "netcdf_cmd.h"

int CMD_Parser::next()
{
	if (cur_ == argc_) return 0;
	yylval.str = argv_[cur_];
	if (!strcmp(argv_[cur_], "help")) {
		return HELP;
	} else if (!strcmp(argv_[cur_], "info")) {
		return INFO;
	} else {
		return yylval.str[0];
	}
	cur_ ++;
}

