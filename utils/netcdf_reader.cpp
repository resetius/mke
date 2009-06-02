#include "netcdf_cmd.h"
#include "netcdf_cmd.hpp"

void CMD_Parser::info()
{
	if (!f_)
	{
		fprintf (stderr, "give me filename first!\n");
		exit (-1);
	}

	int dims = f_->num_dims();
	int vars = f_->num_vars();
	int atts = f_->num_atts();

	fprintf (stdout, "dims/vars/atts: %d/%d/%d\n", dims, vars, atts);
	fprintf (stdout, "dims:\n");
	for (int i = 0; i < dims; ++i)
	{
		NcDim * dim = f_->get_dim (i);
		if (dim) fprintf (stdout, "%d:\t%s: %ld\n", i, dim->name(), dim->size() );
	}

	fprintf (stdout, "vars:\n");
	for (int i = 0; i < vars; ++i)
	{
		NcVar * var = f_->get_var (i);
		if (var) fprintf (stdout, "%d:\t%s\n", i, var->name() );
	}

	fprintf (stdout, "atts:\n");
	for (int i = 0; i < vars; ++i)
	{
		NcAtt * att = f_->get_att (i);
		if (att) fprintf (stdout, "%d:\t%s\n", i, att->name() );
	}
}

void CMD_Parser::open (const char * file)
{
	if (f_)
	{
		delete f_;
		f_ = 0;
	}
	f_ = new NcFile (file);
	if (f_->is_valid() )
	{
		fprintf (stderr, "bad file %s\n", file);
		exit (-1);
	}
}

int CMD_Parser::next()
{
	if (cur_ == argc_) return 0;
	int ans = 0;

	yylval.str = argv_[cur_];

	if (!strcmp (argv_[cur_], "help") )
	{
		ans = HELP;
	}
	else if (!strcmp (argv_[cur_], "info") )
	{
		ans = INFO;
	}
	else if (!strcmp (argv_[cur_], "file") )
	{
		ans = FLE;
	}
	else
	{
		int i = atoi (argv_[cur_]);
		if (i != 0)
		{
			yylval.num = i;
			ans = NUMBER;
		}
		else if (!strcmp (argv_[cur_], "0") )
		{
			yylval.num = 0;
			ans = NUMBER;
		}
		else
		{
			ans = STR;
		}
	}
	cur_ ++;
	return ans;
}

void CMD_Parser::parse()
{
	while (yyparse (this) );
}

int main (int argc, char * argv[])
{
	CMD_Parser p (argc, argv);
	p.parse();
	return 0;
}

