#ifndef NETCDF_CMD_H
#define NETCDF_CMD_H

#include <string.h>

struct CMD_Parser {
	int cur_;
	int argc_;
	char ** argv_;

	CMD_Parser (int argc, char * argv[]): cur_(0), argc_(argc), argv_(argv) {}
	int next();
};

int yyparse(CMD_Parser * );

#endif /* NETCDF_CMD_H */

