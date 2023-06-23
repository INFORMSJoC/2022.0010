#include "CreatIntVector.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
int *create_int_vector(int dim)
{
	int *ptr;

	if ((ptr = (int *)calloc(dim, sizeof(int))) == NULL) {
		int ok = 4;
		printf("\nError: Insuficient memory %d\n", ok);
		exit(8);
	}
	return ptr;
}