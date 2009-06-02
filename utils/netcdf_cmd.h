#ifndef NETCDF_CMD_H
#define NETCDF_CMD_H

#include <string.h>

struct CMD_Parser {
	int cur_;
	int argc_;
	char ** argv_;

	CMD_Parser (int argc, char * argv[]): cur_(0), argc_(argc), argv_(argv) {}
	int next();

	void info_dim(int number) {}
	void info_dim(const char * name) {}

	void info_var(int number) {}
	void info_var(const char * name) {}

	void info_att(int number) {}
	void info_att(const char * name) {}
};

int yyparse(CMD_Parser * );

#endif /* NETCDF_CMD_H */

