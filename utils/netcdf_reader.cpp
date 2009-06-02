#include <netcdf.hh>

void info_file(const char * path)
{
	NcFile f(path);
	if (!f.is_valid()) {
		fprintf(stderr, "bad file %s\n", path);
		return;
	}

	int dims = f.num_dims();
	int vars = f.num_vars();
	int atts = f.num_atts();

	fprintf(stdout, "dims/vars/atts: %d/%d/%d\n", dims, vars, atts);
	fprintf(stdout, "dims:\n");
	for (int i = 0; i < dims; ++i) {
		NcDim * dim = f.get_dim(i);
		if (dim) fprintf(stdout, "%d:\t%s: %ld\n", i, dim->name(), dim->size());
	}

	fprintf(stdout, "vars:\n");
	for (int i = 0; i < vars; ++i) {
		NcVar * var = f.get_var(i);
		if (var) fprintf(stdout, "%d:\t%s\n", i, var->name());
	}

	fprintf(stdout, "atts:\n");
	for (int i = 0; i < vars; ++i) {
		NcAtt * att = f.get_att(i);
		if (att) fprintf(stdout, "%d:\t%s\n", i, att->name());
	}
}

int main(int argc, char * argv[])
{
	info_file(argv[1]);
	return 0;
}

