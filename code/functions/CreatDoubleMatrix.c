#include "CreatDoubleVector.h"
#include "CreatDoubleMatrix.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
double **create_double_matrix(int rows, int Columns)
{
	int i;
	double **ptr;


	if ((ptr = (double **)calloc(rows, sizeof(double *))) == NULL) {
		int ok = 3;
		printf("\nError: Insuficient memory %d\n", ok);
		exit(8);
	}
	for (i = 0; i<rows; i++) {
		ptr[i] = create_double_vector(Columns);
	}
	return ptr;
}
