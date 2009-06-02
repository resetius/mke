#include <sstream>
#include <stdlib.h>

#include "netcdf_cmd.h"
#include "netcdf_cmd.hpp"

using namespace std;

void CMD_Parser::help()
{
	fprintf(stderr, "help message\n");
}

void CMD_Parser::check_file()
{
	if (!f_)
	{
		fprintf (stderr, "give me filename first!\n");
		exit (-1);
	}
}

void CMD_Parser::info()
{
	check_file();

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
	for (int i = 0; i < atts; ++i)
	{
		NcAtt * att = f_->get_att (i);
		if (att) {
			fprintf (stdout, "%d:\t%s\n", i, att->name() );
			info_att(att);
		}
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
	if (!f_->is_valid())
	{
		fprintf (stderr, "bad file %s\n", file);
		exit (-1);
	}
}

void nctype_print(NcType type)
{
	switch (type) {
	case ncByte:
		fprintf(stdout, "byte\n");
		break;
	case ncChar:
		fprintf(stdout, "char\n");
		break;
	case ncShort:
		fprintf(stdout, "short\n");
		break;
	case ncInt:
		fprintf(stdout, "int\n");
		break;
	case ncFloat:
		fprintf(stdout, "float\n");
		break;
	case ncDouble:
		fprintf(stdout, "double\n");
		break;
	default:
		fprintf(stdout, "unknown\n");
		break;
	}
}

void CMD_Parser::info_var(int n)
{
	check_file();
	info_var(f_->get_var(n));
}

void CMD_Parser::info_var(const char * n)
{
	check_file();
	info_var(f_->get_var(n));
}

void CMD_Parser::info_var(NcVar * var)
{
	if (!var) {
		fprintf(stderr, "variable not found!\n");
		exit(-1);
	}

	int dims = var->num_dims();
	int atts = var->num_atts();

	fprintf (stdout, "dims/atts: %d/%d/%d\n", dims, atts);
	fprintf (stdout, "dims:\n");
	for (int i = 0; i < dims; ++i)
	{
		NcDim * dim = var->get_dim (i);
		if (dim) fprintf (stdout, "%d:\t%s: %ld\n", i, dim->name(), dim->size() );
	}

	fprintf (stdout, "atts:\n");
	for (int i = 0; i < atts; ++i)
	{
		NcAtt * att = var->get_att (i);
		if (att) { 
			fprintf (stdout, "%d:\t%s\n", i, att->name() );
			info_att(att);
		}
	}
}

void CMD_Parser::info_att(int num)
{
	check_file();
	info_att(f_->get_att(num));
}

void CMD_Parser::info_att(const char * name)
{
	check_file();
	info_att(f_->get_att(name));
}

void CMD_Parser::info_att(NcAtt * att)
{
	if (!att) {
		fprintf(stderr, "attribute not found!\n");
		exit(-1);
	}

	NcType type = att->type();
	int vals = att->num_vals();

	fprintf(stdout, "\tNcType: ");
	nctype_print(type);
	fprintf(stdout, "\tvals: %d\n", vals);

	NcValues * val = att->values();

	ostringstream str;
	val->print(str);
	fprintf(stdout, "\t%s\n", str.str().c_str());
}

template < typename T >
void get_and_write(int count, const long * counts, FILE * f, NcVar * var)
{
	T * d = new T[count];
	var->get(d, counts);
	fwrite(d, sizeof(T), count, f);
	delete [] d;
}

void CMD_Parser::dump(const char * to, const char * what)
{
	check_file();

	NcVar * var = f_->get_var(what);
	long from[]  = {0, 0, 0, 0, 0};
	long total[] = {0, 0, 0, 0, 0}; // set to max!!!

	if (!var) {
		fprintf(stderr, "variable not found!\n");
		exit(-1);
	}

	fprintf(stderr, "dump %s to %s\n", what, to);

	int dims  = var->num_dims();
	int elems = 1;
	for (int i = 0; i < dims & i < 5; ++i) {
		NcDim * dim = var->get_dim(i);
		if (dim) {
			Slice fake(f_, dim->name(), -1, -1);
			slices_t::iterator it = slices_.find(fake);
			if (it != slices_.end()) {
				from[i]  = it->from;
				total[i] = it->total;
			} else {
				total[i] = dim->size();
			}

			elems *= total[i];
		}
	}

	var->set_cur(from);
	
	FILE * f = fopen(to, "wb");
	NcType type = var->type();
	switch (type) {
	case ncByte:
		get_and_write < ncbyte > (elems, total, f, var);
		break;
	case ncChar:
		get_and_write < char > (elems, total, f, var);
		break;
	case ncShort:
		get_and_write < short > (elems, total, f, var);
		break;
	case ncLong:
		get_and_write < long > (elems, total, f, var);
		break;
	case ncFloat:
		get_and_write < float > (elems, total, f, var);
		break;
	case ncDouble:
		get_and_write < double > (elems, total, f, var);
		break;
	default:
		fprintf(stderr, "unknown type!\n");
		exit(-1);
		break;
	}
	fclose(f);
}

void CMD_Parser::add_slice(const char * dim)
{
	check_file();
	slices_.insert(Slice(f_, dim, -1, -1));
}

void CMD_Parser::add_slice(const char * dim, double from, double to)
{
	check_file();
	slices_.insert(Slice(f_, dim, from, to));
}

void CMD_Parser::reset_slice()
{
	slices_.clear();
}

int CMD_Parser::next()
{
	int ans = 0;
	if (cur_ == argc_) {
		return ans;
	}

	yylval.str = argv_[cur_];

	if (!strcmp (argv_[cur_], "-help") || !strcmp (argv_[cur_], "-h") || !strcmp (argv_[cur_], "--help"))
	{
		ans = HELP;
	}
	else if (!strcmp (argv_[cur_], "-info") )
	{
		ans = INFO;
	}
	else if (!strcmp (argv_[cur_], "-file") )
	{
		ans = FLE;
	}
	else if (!strcmp (argv_[cur_], "-dump") )
	{
		reset_slice();
		ans = DUMP;
	}
	else if (!strcmp (argv_[cur_], "-var") )
	{
		ans = VAR;
	}
	else if (!strcmp (argv_[cur_], "-att") )
	{
		ans = ATT;
	}
	else if (!strcmp (argv_[cur_], "-dim") )
	{
		ans = DIM;
	}
	else
	{
		char * s = argv_[cur_];
		int flag = 1;
		while (*s) {
			if (!(('0' <= *s) && (*s <= '9') || *s == '.')) {
				flag = 0;
				break;
			}
			s++;
		}


		if (flag)
		{
			yylval.num = atof(argv_[cur_]);
			ans = NUM;
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
	yyparse (this);
}

int main (int argc, char * argv[])
{
	CMD_Parser p (argc, argv);
	p.parse();
	return 0;
}

