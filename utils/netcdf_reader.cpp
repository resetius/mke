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
		fprintf(stdout, "%s: %ld\n", dim->name(), dim->size());
	}

	fprintf(stdout, "vars:\n");
	for (int i = 0; i < vars; ++i) {
		NcVar * var = f.get_var(i);
		fprintf(stdout, "%s\n", var->name());
	}

	fprintf(stdout, "atts:\n");
	for (int i = 0; i < vars; ++i) {
		NcAtt * att = f.get_att(i);
		fprintf(stdout, "%s\n", att->name());
	}
}

int main(int argc, char * argv[])
{
	info_file(argv[0]);
	return 0;
}

