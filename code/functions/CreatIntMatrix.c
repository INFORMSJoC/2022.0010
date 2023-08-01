#include "CreatIntVector.h"
#include "CreatIntMatrix.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
int **create_int_matrix(int rows, int Columns)
{
	int i;
	int **ptr;

	if ((ptr = (int **)calloc(rows, sizeof(int *))) == NULL) {
		int ok = 1;
		printf("\nError: Insuficient memory %d\n", ok);
		exit(8);
	}
	for (i = 0; i<rows; i++)
		ptr[i] = create_int_vector(Columns);
	return ptr;
}
