#include "CreatDoubleVector.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
double *create_double_vector(int dim)
{
	double *ptr;

	if ((ptr = (double *)calloc(dim, sizeof(double))) == NULL) {
		int ok = 5;
		printf("\nError: Insuficient memory %d\n", ok);
		exit(8);
	}
	return ptr;
}