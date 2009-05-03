#include <stdio.h>
#include <windows.h>

int main(int argc, char * argv[])
{
	FILE * f = fopen(argv[1], "r");
	char buf[32768];

	while (fgets(buf, 32768, f)) {
		int slp = 0;
		if (strstr(buf, "# end")) {
			slp = 1;
		}
		fputs(buf, stdout);
		if (slp) {
			Sleep(100);
		}
	}	
}
