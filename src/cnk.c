#include <stdio.h>

int fk(int n)
{
	int i;
	int f;
	if (n == 0) return 1;

	f = 1;
	for (i = 1; i <= n; ++i) {
		f = f * i;
	}

	return f;
}

int main(int agrc, char * argv[])
{
	int m = 10;
	int n, k;

	printf("const static double Cnk[][%d] = \n{ \n", m);

	for (n = 0; n < m; ++n) {
		printf("{ ");
		for (k = 0; k < m; ++k) {
			double cnk = (double)fk(n) / (double)(fk(k) * fk(n - k));
			printf("%.16le, ", cnk);
		}
		printf("},\n");
	}

	printf("};\n");
	return 0;
}

