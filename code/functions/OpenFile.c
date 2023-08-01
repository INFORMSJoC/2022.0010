#include "OpenFile.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
FILE *Open_File(const char *name, const char *mode)  // This function opens the file "name" under mode "mode" (reading, writing, etc)
{
	FILE *file;
	int OK;
	if ((file = fopen(name, mode)) == NULL) {
		OK = 1;
		printf("\nError: File cannot be opened %d\n", OK);

	}
	return file;
}

