#ifndef NETCDF_CMD_H
#define NETCDF_CMD_H

#include <netcdf.hh>
#include <string.h>

struct CMD_Parser {
	int cur_;
	int argc_;
	char ** argv_;
	NcFile * f_;

	CMD_Parser (int argc, char * argv[]): f_(0), cur_(0), argc_(argc), argv_(argv) {}
	~CMD_Parser() { delete f_; }

	void parse();
	void open(const char * name);
	int next();

	void info();

	void info_dim(int number) {}
	void info_dim(const char * name) {}

	void info_var(int number) {}
	void info_var(const char * name) {}

	void info_att(int number) {}
	void info_att(const char * name) {}
};

int yyparse(CMD_Parser * );

#endif /* NETCDF_CMD_H */

