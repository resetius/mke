#ifndef NETCDF_CMD_H
#define NETCDF_CMD_H

#include <netcdf.hh>
#include <string.h>

struct CMD_Parser {
	int cur_;
	int argc_;
	char ** argv_;
	NcFile * f_;

	CMD_Parser (int argc, char * argv[]): f_(0), cur_(1), argc_(argc), argv_(argv) {}
	~CMD_Parser() { delete f_; }

	void parse();
	void open(const char * name);
	int next();

	void help();
	void info();

	void info_var(int number);
	void info_var(const char * name);
	void info_var(NcVar *);

	void info_att(int number);
	void info_att(const char * name);
	void info_att(NcAtt *);

	void check_file();

	void dump(const char * to, const char * what);
	void add_slice(const char * dim);
	void add_slice(const char * dim, double from, double to);
	void reset_slice();
};

int yyparse(CMD_Parser * );

#endif /* NETCDF_CMD_H */

