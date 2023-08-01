#include "CreatIntVector.h"
#include "CreatIntMatrix.h"
#include "CreatIntTmatrix.h"
#include "CreatIntFmatrix.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
int ****create_int_fmatrix(int fou, int thr, int rows, int Columns)
{
	int i;
	int ****ptr;

	if ((ptr = (int ****)calloc(fou, sizeof(int ***))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i<fou; i++)
		ptr[i] = create_int_tmatrix(thr, rows, Columns);
	return ptr;
}

