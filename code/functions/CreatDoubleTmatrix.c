#include "CreatDoubleVector.h"
#include "CreatDoubleMatrix.h"
#include "CreatDoubleTmatrix.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
double ***create_double_tmatrix(int thr, int rows, int Columns)
{
	int i;
	double ***ptr;

	if ((ptr = (double ***)calloc(thr, sizeof(double **))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i<thr; i++)
		ptr[i] = create_double_matrix(rows, Columns);
	return ptr;
}