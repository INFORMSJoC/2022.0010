#include <ilcplex/cplex.h>

/* Bring in the declarations for the string and character functions,
malloc, and fabs. */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include "solveNVSAA.h"
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define EPSILON	 0.000000000001    //Added to the fractional solution obtained at the node
void free_and_null(char **ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */
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

double solveNVSAA(int I, int S, double *p, double **z, double *c, double d, int Moos, double **xis_oos) {
	int *x, **y0;
	x = create_int_vector(I);
	y0 = create_int_matrix(I, S);
	int i, s;
	int status;
	CPXENVptr env = NULL;
	CPXLPptr  lp = NULL;

	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		fprintf(stderr, "Failure in CPXopenCPLEX, status = %d.\n", status);
		
	}

	/* Create the master ILP */

	lp = CPXcreateprob(env, &status, "master_ILP.lp");
	if (lp == NULL) {
		fprintf(stderr, "Failure in CPXcreateprob, status = %d.\n", status);
		
	}
	char sense;
	int nzcnt, rmatbeg, *rmatind = NULL;
	double rhs, *rmatval = NULL;
	int index1 = 0;

	for (i = 0; i < I; i++) {
		x[i] = index1;
		index1++;
		double lb = 0.;
		double ub = CPX_INFBOUND;
		double cost = -1.*p[i]+c[i];
		status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			
		}
	}
	int num_x_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		for (s = 0; s < S; s++) {
			y0[i][s] = num_x_var + index1;
			index1++;
			double lb = 0.;
			double ub = CPX_INFBOUND;
			double cost = p[i] / S;
			status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				
			}
		}
	}
	int num_y0_var = index1;
	index1 = 0;

	rmatbeg = 0;

	rmatind = (int *)malloc(I * sizeof(int));
	
	rmatval = (double *)malloc(I * sizeof(double));
	
	sense = 'L';
	rhs = d;
	nzcnt = 0;

	for (i = 0; i < I; i++) {
		rmatind[nzcnt] = x[i];
		rmatval[nzcnt++] = c[i];
	}

	status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
		&rmatbeg, rmatind, rmatval, NULL, NULL);
	
	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);


	rmatind = (int *)malloc(2 * sizeof(int));
	
	rmatval = (double *)malloc(2 * sizeof(double));
	

	sense = 'L';
	for (s = 0; s < S; s++) {
		for (i = 0; i < I; i++) {
			rhs = z[i][s];
			nzcnt = 0;

			rmatind[nzcnt] = y0[i][s];
			rmatval[nzcnt++] = -1.;
			rmatind[nzcnt] = x[i];
			rmatval[nzcnt++] = 1.;

			status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
				&rmatbeg, rmatind, rmatval, NULL, NULL);
			
		}
	}
	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);

	CPXsetdblparam(env, CPX_PARAM_TILIM, 7200);
	CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);

	status = CPXlpopt(env, lp);
	
	int worker_lp_sol_stat = CPXgetstat(env, lp);
	int num_sub = CPXgetnumcols(env, lp);
	//printf("saa status is %d\t", worker_lp_sol_stat);
	double obj_oos = 0.;
	if (worker_lp_sol_stat == CPX_STAT_OPTIMAL) {
		double *x2;
		x2 = create_double_vector(num_sub);
		
		status = CPXgetx(env, lp, x2, 0, num_sub - 1);
		
		for (i = 0; i < I; i++) {
			for (s = 0; s < Moos; s++) {
				obj_oos += p[i] * MIN(x2[x[i]], xis_oos[i][s]) - c[i] * x2[x[i]];
			}
		}
		obj_oos = obj_oos / Moos;
		//printf("saa obj is %.2lf\t", obj_oos);
		free_and_null((char **)&x2);
	}
	if (lp != NULL) {
		int local_status = CPXfreeprob(env, &lp);
		if (local_status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n",
				local_status);
			status = local_status;
		}
	}

	/* Free the CPLEX environment, if necessary */

	if (env != NULL) {
		int local_status = CPXcloseCPLEX(&env);
		if (local_status) {
			fprintf(stderr,
				"Could not close CPLEX environment, status = %d.\n",
				local_status);
			status = local_status;
		}
	}
	/**/
	free_and_null((char **)&x);
	for (i = 0; i < I; i++) {
		free_and_null((char **)&y0[i]);
	}
	free_and_null((char **)&y0);
	
	return obj_oos;
}
