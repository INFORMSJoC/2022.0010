#include <ilcplex/cplex.h>
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
#include "solveNVDROW.h"
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define EPSILON	 0.000000000001    //Added to the fractional solution obtained at the node
void release_list(linklist *pHead)
{
	linklist *p = pHead;
	linklist *q;
	while (p)
	{
		q = p->next;
		free(p);
		p = q;
	}

	return;
}
int ***create_int_tmatrix(int thr, int rows, int Columns)
{
	int i;
	int ***ptr;

	if ((ptr = (int ***)calloc(thr, sizeof(int **))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i<thr; i++)
		ptr[i] = create_int_matrix(rows, Columns);
	return ptr;
}
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

double solveNVDROW(int I, int S, double *uz, double theta, double *p, double **z, double *c, double d, int Moos, double **xis_oos, double *prob, double *uu, double *prob_valid)//(int I, int S, double *ubar, double Gamma, double *bs, double **xis, double *cs, double d, int Moos, double **xis_oos) {
{	
	linklist *Ly1, *pl, *ql, *pl1;
	linklist *Ly11, *pl11;
	linklist *Leta1, *ple1;
	linklist *Leta2, *ple2;
	linklist *Leta3, *ple3;
	linklist *Leta4, *ple4;
	linklist *Leta5, *ple5;
	linklist *Leta6, *ple6;
	double obj_oos = 0.;

	Ly1 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Ly1->next = NULL;

	pl = Ly1;
	Ly11 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Ly11->next = NULL;

	pl11 = Ly11;

	Leta1 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Leta1->next = NULL;
	ple1 = Leta1;

	Leta2 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Leta2->next = NULL;
	ple2 = Leta2;

	Leta3 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Leta3->next = NULL;
	ple3 = Leta3;

	Leta4 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Leta4->next = NULL;
	ple4 = Leta4;

	Leta5 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Leta5->next = NULL;
	ple5 = Leta5;

	Leta6 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Leta6->next = NULL;
	ple6 = Leta6;
	double min_y1 = -1000.;
	double min_y11 = -1000.;
	double min_eta1 = -1000.;
	double min_eta2 = -1000.;
	double min_eta3 = -1000.;
	double min_eta4 = -1000.;
	double min_eta5 = -1000.;
	double min_eta6 = -1000.;

	int *x, alpha, *r, t, **y0, ***y1, **y2;
	int **eta1, **eta2, ***eta3, ***eta4, ***eta5, ***eta6;
	int *w1, **w2, **w3;
	double *tau;

	int tau1, tau2, *tau3, *tau4, **tau5, **tau6, **tau7, **tau8, ***tau9, ***tau10, **tau11, **tau12, ***tau13, ***tau14;
	tau3 = create_int_vector(S);
	tau4 = create_int_vector(S);
	tau5 = create_int_matrix(S, I);
	tau6 = create_int_matrix(S, I);
	tau7 = create_int_matrix(S, I);
	tau8 = create_int_matrix(S, I);
	tau9 = create_int_tmatrix(S, I, I);
	tau10 = create_int_tmatrix(S, I, I);
	tau11 = create_int_matrix(S, I);
	tau12 = create_int_matrix(S, I);
	tau13 = create_int_tmatrix(S, I, I);
	tau14 = create_int_tmatrix(S, I, I);
	int tau21, *tau31, *tau41, **tau51, **tau61, **tau71, **tau81, ***tau91, ***tau101, **tau111, **tau121, ***tau131, ***tau141;
	tau31 = create_int_vector(S);
	tau41 = create_int_vector(S);
	tau51 = create_int_matrix(S, I);
	tau61 = create_int_matrix(S, I);
	tau71 = create_int_matrix(S, I);
	tau81 = create_int_matrix(S, I);
	tau91 = create_int_tmatrix(S, I, I);
	tau101 = create_int_tmatrix(S, I, I);
	tau111 = create_int_matrix(S, I);
	tau121 = create_int_matrix(S, I);
	tau131 = create_int_tmatrix(S, I, I);
	tau141 = create_int_tmatrix(S, I, I);
	int indexd;
	x = create_int_vector(I);
	r = create_int_vector(S);
	y0 = create_int_matrix(I, S);
	y1 = create_int_tmatrix(I, I, S);
	y2 = create_int_matrix(I, S);
	eta1 = create_int_matrix(I, S);
	eta2 = create_int_matrix(I, S);
	eta3 = create_int_tmatrix(I, I, S);
	eta4 = create_int_tmatrix(I, I, S);
	eta5 = create_int_tmatrix(I, I, S);
	eta6 = create_int_tmatrix(I, I, S);
	w1 = create_int_vector(S);
	w2 = create_int_matrix(I, S);
	w3 = create_int_matrix(I, S);
	int **y21 = NULL;
	y21 = create_int_matrix(I, S);
	int ***y11;
	y11 = create_int_tmatrix(I, I, S);

	CPXENVptr env = NULL;
	CPXLPptr  lp = NULL;

	int iter, Iter;
	iter = 0;
	Iter = 1;
	int flag = 0;
	double avg_time1_iter = 0.;
	int status = 0;
	int i, j, s;
	linklist *Ly2, *ply2;
	linklist *Ly21, *ply21;
	Ly2 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Ly2->next = NULL;
	ply2 = Ly2;

	Ly21 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
	Ly21->next = NULL;
	ply21 = Ly21;
	while (flag == 0 && iter < Iter &&avg_time1_iter<7200) {
		env = CPXopenCPLEX(&status);
		if (env == NULL) {
			fprintf(stderr, "Failure in CPXopenCPLEX, status = %d.\n", status);
			goto TERMINATE;
		}

		/* Set MIP log interval to 1

		status = CPXsetintparam(env, CPXPARAM_MIP_Interval, 10000);
		if (status) {
		fprintf(stderr,
		"Failed to set CPXPARAM_MIP_Interval, status = %d.\n", status);
		goto TERMINATE;
		}*/

		/* Create the master ILP */

		lp = CPXcreateprob(env, &status, "master_ILP.lp");
		if (lp == NULL) {
			fprintf(stderr, "Failure in CPXcreateprob, status = %d.\n", status);
			goto TERMINATE;
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
				goto TERMINATE;
			}
		}
		int num_x_var = index1;
		index1 = 0;

		alpha = num_x_var + index1;
		index1++;
		double lb = 0.;
		double ub = CPX_INFBOUND;
		double cost = 1.;
		status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		int num_alpha_var = index1;
		index1 = 0;

		for (s = 0; s < S; s++) {
			r[s] = num_x_var + num_alpha_var + index1;
			index1++;
			lb = 0.;
			ub = CPX_INFBOUND;
			status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
			}
		}
		int num_r_var = index1;
		index1 = 0;

		t = num_r_var + num_x_var + num_alpha_var + index1;
		index1++;
		lb = 0.;
		ub = CPX_INFBOUND;
		status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		int num_t_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				y0[i][s] = num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_y0_var = index1;
		index1 = 0;


		double *lb_y1, *ub_y1;
		lb_y1 = create_double_vector(I*I*S);
		ub_y1 = create_double_vector(I*I*S);
		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					y1[i][j][s] = num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;
					//lb = -CPX_INFBOUND;
					//ub = CPX_INFBOUND;

				}
			}
		}
		int num_y1_var = index1;
		index1 = 0;
		for (i = 0; i < I*I*S; i++) {
			lb_y1[i] = 0.;
			ub_y1[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*I*S, NULL, lb_y1, ub_y1, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}

		for (j = 0; j < I; j++) {
			for (s = 0; s < S; s++) {
				y2[j][s] = num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;
				/*	lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
				}*/
			}
		}
		int num_y2_var = index1;
		index1 = 0;

		double *lb_y2, *ub_y2;
		lb_y2 = create_double_vector(I*S);
		ub_y2 = create_double_vector(I*S);


		for (i = 0; i < I*S; i++) {
			lb_y2[i] = 0.;
			ub_y2[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*S, NULL, lb_y2, ub_y2, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				eta1[i][s] = num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;
				/*	lb = 0.;
				ub = 0.;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
				}*/
			}
		}
		int num_eta1_var = index1;
		index1 = 0;
		for (i = 0; i < I*S; i++) {
			lb_y2[i] = 0.;
			ub_y2[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*S, NULL, lb_y2, ub_y2, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				eta2[i][s] = num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;
				/*	lb = 0.;
				ub = 0.;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
				}*/
			}
		}
		int num_eta2_var = index1;
		index1 = 0;
		for (i = 0; i < I*S; i++) {
			lb_y2[i] = 0.;
			ub_y2[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*S, NULL, lb_y2, ub_y2, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta3[i][j][s] = num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;
					/*	lb = 0.;
					ub = CPX_INFBOUND;
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
					}*/
				}
			}
		}
		int num_eta3_var = index1;
		index1 = 0;
		for (i = 0; i < I*I*S; i++) {
			lb_y1[i] = 0.;
			ub_y1[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*I*S, NULL, lb_y1, ub_y1, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta4[i][j][s] = num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;
					/*	lb = 0.;
					ub = CPX_INFBOUND;
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
					}*/
				}
			}
		}
		int num_eta4_var = index1;
		index1 = 0;
		for (i = 0; i < I*I*S; i++) {
			lb_y1[i] = 0.;
			ub_y1[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*I*S, NULL, lb_y1, ub_y1, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta5[i][j][s] = num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;
					/*	lb = 0.;
					ub = CPX_INFBOUND;
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
					}*/
				}
			}
		}
		int num_eta5_var = index1;
		index1 = 0;
		for (i = 0; i < I*I*S; i++) {
			lb_y1[i] = 0.;
			ub_y1[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*I*S, NULL, lb_y1, ub_y1, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta6[i][j][s] = num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;
					/*	lb = 0.;
					ub = CPX_INFBOUND;
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
					}*/
				}
			}
		}
		int num_eta6_var = index1;
		index1 = 0;
		for (i = 0; i < I*I*S; i++) {
			lb_y1[i] = 0.;
			ub_y1[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*I*S, NULL, lb_y1, ub_y1, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}

		for (s = 0; s < S; s++) {
			w1[s] = num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
			index1++;
			lb = 0.;
			ub = CPX_INFBOUND;
			status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
			}
		}

		int num_w1_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				w2[i][s] = num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_w2_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				w3[i][s] = num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_w3_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					y11[i][j][s] = num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;
					//	lb = 0.;
					//	ub = CPX_INFBOUND;

					//	status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					//	if (status) {
					//		fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					//		goto TERMINATE;
					//	}
				}
			}
		}
		int num_y11_var = index1;
		index1 = 0;
		for (i = 0; i < I*I*S; i++) {
			lb_y1[i] = 0.;
			ub_y1[i] = CPX_INFBOUND;
		}
		
		status = CPXnewcols(env, lp, I*I*S, NULL, lb_y1, ub_y1, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		int alpha1,*r1,**y01;
		r1 = create_int_vector(S);
		y01 = create_int_matrix(I, S);
		
		alpha1 = num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
		index1++;
		lb = 0.;
		ub = CPX_INFBOUND;
		cost = -1;
		status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		int num_alpha1_var = index1;
		index1 = 0;

		for (s = 0; s < S; s++) {
			r1[s] = num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
			index1++;
			lb = 0.;
			ub = CPX_INFBOUND;
			status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
			}
		}
		int num_r1_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				y01[i][s] = num_r1_var + num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_y01_var = index1;
		index1 = 0;

		int num_tau1 = num_y01_var + num_r1_var + num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;

		indexd = 0;

		tau2 = num_tau1 + indexd;
		tau21 = num_tau1 + indexd;
		indexd++;
		lb = 0.;
		ub = CPX_INFBOUND;
		status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
		int num_tau2 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			tau3[s] = num_tau2 + num_tau1 + indexd;
			tau31[s] = num_tau2 + num_tau1 + indexd;
			indexd++;
			lb = 0.;
			ub = CPX_INFBOUND;
			status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
			}
		}
		int num_tau3 = indexd;
		indexd = 0;


		for (s = 0; s < S; s++) {
			tau4[s] = num_tau3 + num_tau2 + num_tau1 + indexd;
			tau41[s] = num_tau3 + num_tau2 + num_tau1 + indexd;
			indexd++;
			lb = 0.;
			ub = CPX_INFBOUND;
			status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
			}
		}
		int num_tau4 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (j = 0; j < I; j++) {
				tau5[s][j] = num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				tau51[s][j] = num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_tau5 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (j = 0; j < I; j++) {
				tau6[s][j] = num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				tau61[s][j] = num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_tau6 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau7[s][i] = num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				tau71[s][i] = num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_tau7 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau8[s][i] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				tau81[s][i] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_tau8 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau9[s][i][j] = num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					tau91[s][i][j] = num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
					lb = 0.;
					ub = CPX_INFBOUND;
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		int num_tau9 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau10[s][i][j] = num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					tau101[s][i][j] = num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
					lb = 0.;
					ub = CPX_INFBOUND;
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		int num_tau10 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau11[s][i] = num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				tau111[s][i] = num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_tau11 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau12[s][i] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				tau121[s][i] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;
				lb = 0.;
				ub = CPX_INFBOUND;
				status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		int num_tau12 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau13[s][i][j] = num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					tau131[s][i][j] = num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
					lb = 0.;
					ub = CPX_INFBOUND;
					cost = 0.;
					if (i == j) {
						cost = -1.;
					}
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		int num_tau13 = indexd;
		indexd = 0;


		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau14[s][i][j] = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					tau141[s][i][j] = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
					lb = 0.;
					ub = CPX_INFBOUND;
					cost = 0.;
					if (i == j) {
						cost = 1.;
					}

					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		int num_tau14 = indexd;
		indexd = 0;

		rmatbeg = 0;

		rmatind = (int *)malloc(I * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc(I * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'L';
		rhs = d;
		nzcnt = 0;

		for (i = 0; i < I; i++) {
			rmatind[nzcnt] = x[i];
			rmatval[nzcnt++] = c[i];
		}

		status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
			&rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
			goto TERMINATE;
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);


		rmatind = (int *)malloc((S * 2 + 4) * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc((S * 2 + 4) * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		nzcnt = 0;
		rmatind[nzcnt] = tau2;
		rmatval[nzcnt++] = -1.;

		rmatind[nzcnt] = alpha;
		rmatval[nzcnt++] = 1.;
		rmatind[nzcnt] = alpha1;
		rmatval[nzcnt++] = -1.;
		for (s = 0; s < S; s++) {
			rmatind[nzcnt] = r[s];
			rmatval[nzcnt++] = -1. * prob[s];
		}
		for (s = 0; s < S; s++) {
			rmatind[nzcnt] = r1[s];
			rmatval[nzcnt++] = prob[s];
		}
		rmatind[nzcnt] = t;
		rmatval[nzcnt++] = -1.*theta;
		status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
			&rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
			goto TERMINATE;
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc((I*(4 + I) + 3 + I * I) * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc((I*(4 + I) + 3 + I * I) * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			nzcnt = 0;
			rmatind[nzcnt] = tau3[s];
			rmatval[nzcnt++] = 1.;

			rmatind[nzcnt] = r[s];
			rmatval[nzcnt++] = -1.;
			rmatind[nzcnt] = r1[s];
			rmatval[nzcnt++] = 1.;
			for (i = 0; i < I; i++) {
				rmatind[nzcnt] = y0[i][s];
				rmatval[nzcnt++] = p[i];
			}
			for (i = 0; i < I; i++) {
				rmatind[nzcnt] = y01[i][s];
				rmatval[nzcnt++] = -p[i];
			}
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = p[i] * z[j][s];
				}
			}
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = -1.*p[i] * z[j][s];
				}
			}
			for (i = 0; i < I; i++) {
				rmatind[nzcnt] = eta1[i][s];
				rmatval[nzcnt++] = uz[i] - z[i][s];
			}
			for (i = 0; i < I; i++) {
				rmatind[nzcnt] = eta2[i][s];
				rmatval[nzcnt++] = z[i][s]-uu[i];
			}
			status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
				&rmatbeg, rmatind, rmatval, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
				goto TERMINATE;
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc((I + 3) * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc((I + 3) * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			nzcnt = 0;
			rmatind[nzcnt] = tau4[s];
			rmatval[nzcnt++] = 1.;

			rmatind[nzcnt] = w1[s];
			rmatval[nzcnt++] = 1.;
			for (i = 0; i < I; i++) {
				rmatind[nzcnt] = y2[i][s];
				rmatval[nzcnt++] = p[i];
			}

			rmatind[nzcnt] = t;
			rmatval[nzcnt++] = -1.;

			status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
				&rmatbeg, rmatind, rmatval, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
				goto TERMINATE;
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc((I * 2 + 4) * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc((I * 2 + 4) * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			for (j = 0; j < I; j++) {
				nzcnt = 0;
				rmatind[nzcnt] = tau5[s][j];
				rmatval[nzcnt++] = 1.;

				rmatind[nzcnt] = w1[s];
				rmatval[nzcnt++] = -1.;
				for (i = 0; i < I; i++) {
					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = p[i];
				}
				for (i = 0; i < I; i++) {
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = -1.*p[i];
				}
				rmatind[nzcnt] = eta1[j][s];
				rmatval[nzcnt++] = -1.;
				rmatind[nzcnt] = eta2[j][s];
				rmatval[nzcnt++] = 1.;

				status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
					&rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc((I * 2 + 4) * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc((I * 2 + 4) * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			for (j = 0; j < I; j++) {
				nzcnt = 0;
				rmatind[nzcnt] = tau6[s][j];
				rmatval[nzcnt++] = 1.;
				rmatind[nzcnt] = w1[s];
				rmatval[nzcnt++] = -1.;
				for (i = 0; i < I; i++) {
					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = -1.*p[i];
				}
				for (i = 0; i < I; i++) {
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = p[i];
				}
				rmatind[nzcnt] = eta1[j][s];
				rmatval[nzcnt++] = 1.;
				rmatind[nzcnt] = eta2[j][s];
				rmatval[nzcnt++] = -1.;

				status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
					&rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc((I * 4 + 3) * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc((I * 4 + 3) * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				nzcnt = 0;
				rmatind[nzcnt] = tau7[s][i];
				rmatval[nzcnt++] = 1.;

				rmatind[nzcnt] = y0[i][s];
				rmatval[nzcnt++] = -1.;
				rmatind[nzcnt] = y01[i][s];
				rmatval[nzcnt++] = 1.;
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = -1.*z[j][s];
				}
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = z[j][s];
				}
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = eta3[i][j][s];
					rmatval[nzcnt++] = uz[j] - z[j][s];
				}
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = eta4[i][j][s];
					rmatval[nzcnt++] = z[j][s] - uu[j];
				}
				status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
					&rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc(3 * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc(3 * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				nzcnt = 0;
				rmatind[nzcnt] = tau8[s][i];
				rmatval[nzcnt++] = 1.;
				rmatind[nzcnt] = y2[i][s];
				rmatval[nzcnt++] = -1.;

				rmatind[nzcnt] = w2[i][s];
				rmatval[nzcnt++] = 1.;

				status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
					&rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc(6 * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc(6 * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					nzcnt = 0;
					rmatind[nzcnt] = tau9[s][i][j];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = w2[i][s];
					rmatval[nzcnt++] = -1.;

					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = -1.;
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = eta3[i][j][s];
					rmatval[nzcnt++] = -1.;
					rmatind[nzcnt] = eta4[i][j][s];
					rmatval[nzcnt++] = 1.;
					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc(6 * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc(6 * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					nzcnt = 0;
					rmatind[nzcnt] = tau10[s][i][j];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = w2[i][s];
					rmatval[nzcnt++] = -1.;

					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = -1.;
					rmatind[nzcnt] = eta3[i][j][s];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = eta4[i][j][s];
					rmatval[nzcnt++] = -1.;
					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc((I * 4 + 4) * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc((I * 4 + 4) * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				rhs = z[i][s];
				nzcnt = 0;
				rmatind[nzcnt] = tau11[s][i];
				rmatval[nzcnt++] = 1.;
				rmatind[nzcnt] = y0[i][s];
				rmatval[nzcnt++] = -1.;
				rmatind[nzcnt] = y01[i][s];
				rmatval[nzcnt++] = 1.;
				rmatind[nzcnt] = x[i];
				rmatval[nzcnt++] = 1.;
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = -1.*z[j][s];
				}
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = z[j][s];
				}
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = eta5[i][j][s];
					rmatval[nzcnt++] = uz[j] - z[j][s];
				}
				for (j = 0; j < I; j++) {
					rmatind[nzcnt] = eta6[i][j][s];
					rmatval[nzcnt++] = z[j][s] - uu[j];
				}
				status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
					&rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc(3 * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc(3 * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		rhs = 0.;
		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				nzcnt = 0;
				rmatind[nzcnt] = tau12[s][i];
				rmatval[nzcnt++] = 1.;
				rmatind[nzcnt] = y2[i][s];
				rmatval[nzcnt++] = -1.;
				rmatind[nzcnt] = w3[i][s];
				rmatval[nzcnt++] = 1.;

				status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
					&rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
					goto TERMINATE;
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc(6 * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc(6 * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					nzcnt = 0;
					if (i == j) {
						rhs = 1.;
					}
					else {
						rhs = 0.;
					}
					rmatind[nzcnt] = tau13[s][i][j];
					rmatval[nzcnt++] = 1.;

					rmatind[nzcnt] = w3[i][s];
					rmatval[nzcnt++] = -1.;

					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = -1.;
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = eta5[i][j][s];
					rmatval[nzcnt++] = -1.;
					rmatind[nzcnt] = eta6[i][j][s];
					rmatval[nzcnt++] = 1.;
					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);

		rmatind = (int *)malloc(6 * sizeof(int));
		if (rmatind == NULL) {
			fprintf(stderr, "No memory for rmatind array.\n");
			status = -1;
			goto TERMINATE;
		}
		rmatval = (double *)malloc(6 * sizeof(double));
		if (rmatval == NULL) {
			fprintf(stderr, "No memory for rmatval array.\n");
			status = -1;
			goto TERMINATE;
		}

		sense = 'E';
		//rhs = -1.;
		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					nzcnt = 0;
					if (i == j) {
						rhs = -1.;
					}
					else {
						rhs = 0.;
					}
					rmatind[nzcnt] = tau14[s][i][j];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = w3[i][s];
					rmatval[nzcnt++] = -1.;

					rmatind[nzcnt] = y1[i][j][s];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = y11[i][j][s];
					rmatval[nzcnt++] = -1.;
					rmatind[nzcnt] = eta5[i][j][s];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = eta6[i][j][s];
					rmatval[nzcnt++] = -1.;
					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
		}
		free_and_null((char **)&rmatind);
		free_and_null((char **)&rmatval);
		CPXsetdblparam(env, CPX_PARAM_TILIM, 7200);
		CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		//CPXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);
		//	CPXsetintparam(env, CPXPARAM_Simplex_Display, 2);
		//	CPXsetdblparam(env, CPXPARAM_Barrier_ConvergeTol, 0.01); // e-optimal solution (%gap)

		clock_t start_time1 = clock();//CPXprimopt
		status = CPXlpopt(env, lp);
		if (status) {
			fprintf(stderr, "Error in CPXprimopt: status = %d\n", status);
			goto TERMINATE;
		}
		int worker_lp_sol_stat = CPXgetstat(env, lp);
		int num_sub = CPXgetnumcols(env, lp);
		//printf("status is %d\t", worker_lp_sol_stat);
		if (worker_lp_sol_stat == CPX_STAT_OPTIMAL) {
			clock_t end_time1 = clock();
			avg_time1_iter += (double)(end_time1 - start_time1) / (double)CLOCKS_PER_SEC;

			double z1;
			double *x2;
			x2 = create_double_vector(num_sub);
			status = CPXgetobjval(env, lp, &z1);
			if (status) {
				fprintf(stderr, "Failed to obtain solution, status = %d.\n", status);
				goto TERMINATE;
			}
			
			//get optimal solution
			status = CPXgetx(env, lp, x2, 0, num_sub - 1);
			if (status) {
				fprintf(stderr, "Failed to obtain solution, status = %d.\n", status);
				goto TERMINATE;
			}


			min_y1 = 1000.;
			min_y11 = 1000.;
			int num_row = CPXgetnumrows(env, lp);
			tau = create_double_vector(num_row);
			status = CPXgetpi(env, lp, tau, 0, num_row - 1);
			indexd = 0;

			tau1 = indexd;
			indexd++;
			num_tau1 = indexd;
			indexd = 0;

			tau2 = num_tau1 + indexd;
			indexd++;
			num_tau2 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {

				tau3[s] = num_tau2 + num_tau1 + indexd;
				indexd++;

			}
			num_tau3 = indexd;
			indexd = 0;


			for (s = 0; s < S; s++) {

				tau4[s] = num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;

			}
			num_tau4 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (j = 0; j < I; j++) {
					tau5[s][j] = num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
				}
			}
			num_tau5 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (j = 0; j < I; j++) {
					tau6[s][j] = num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
				}
			}
			num_tau6 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					tau7[s][i] = num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
				}
			}
			num_tau7 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					tau8[s][i] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
				}
			}
			num_tau8 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					for (j = 0; j < I; j++) {
						tau9[s][i][j] = num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
						indexd++;
					}
				}
			}
			num_tau9 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					for (j = 0; j < I; j++) {
						tau10[s][i][j] = num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
						indexd++;
					}
				}
			}
			num_tau10 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					tau11[s][i] = num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
				}
			}
			num_tau11 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					tau12[s][i] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;
				}
			}
			num_tau12 = indexd;
			indexd = 0;

			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					for (j = 0; j < I; j++) {
						tau13[s][i][j] = num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
						indexd++;
					}
				}
			}
			num_tau13 = indexd;
			indexd = 0;


			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					for (j = 0; j < I; j++) {
						tau14[s][i][j] = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
						indexd++;
					}
				}
			}
			num_tau14 = indexd;
			indexd = 0;

			
				for (i = 0; i < I; i++) {
					for (s = 0; s < Moos; s++) {
						obj_oos += prob_valid[s] * p[i] * MIN(x2[x[i]], xis_oos[i][s]);
					}
					obj_oos -= c[i] * x2[x[i]];
					//printf("Xw is %.2lf\t", x2[x[i]]);
				}
				//obj_oos = obj_oos / Moos;
				//printf("\n");
				if (obj_oos < 1.) {
					for (i = 0; i < I; i++) {
						printf("X p is %.2lf\t%.2lf\t", x2[x[i]], p[i]);
					}
					printf("\n");
					for (i = 0; i < I; i++) {
						for (s = 0; s < Moos; s++) {
							printf("out sample is %.2lf\t", xis_oos[i][s]);
						}
					}
					printf("\n");
					for (s = 0; s < Moos; s++) {
						printf("prob is %.2lf\t", prob_valid[s]);
					}
					printf("\n");
				}
			//}

			free_and_null((char **)&x2);
			free_and_null((char **)&tau);

		}
		else {

			avg_time1_iter = 7300.;
		}
		//	solution = Open_File(outfile, "a+");
		//	fprintf(solution, "%.6lf\t%.6lf\t%.6lf\n", cputimeMIP1, fmin_y1, fmin_y11);

		//			fclose(solution);
		iter++;
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
		free_and_null((char **)&lb_y1);
		free_and_null((char **)&ub_y1);
		free_and_null((char **)&lb_y2);
		free_and_null((char **)&ub_y2);
	}
	/**/
	release_list(Ly1);
	release_list(Ly11);
	release_list(Leta1);
	release_list(Leta2);
	release_list(Leta3);
	release_list(Leta4);
	release_list(Leta5);
	release_list(Leta6);

	free_and_null((char **)&x);
	free_and_null((char **)&r);
	free_and_null((char **)&w1);
	free_and_null((char **)&tau3);
	free_and_null((char **)&tau4);
	for (s = 0; s < S; s++) {
		free_and_null((char **)&tau5[s]);
		free_and_null((char **)&tau6[s]);
		free_and_null((char **)&tau7[s]);
		free_and_null((char **)&tau8[s]);
		free_and_null((char **)&tau11[s]);
		free_and_null((char **)&tau12[s]);
	}
	free_and_null((char **)&tau5);
	free_and_null((char **)&tau6);
	free_and_null((char **)&tau7);
	free_and_null((char **)&tau8);
	free_and_null((char **)&tau11);
	free_and_null((char **)&tau12);
	for (s = 0; s < S; s++) {
		for (i = 0; i < I; i++) {
			free_and_null((char **)&tau9[s][i]);
			free_and_null((char **)&tau10[s][i]);
			free_and_null((char **)&tau13[s][i]);
			free_and_null((char **)&tau14[s][i]);
		}
		free_and_null((char **)&tau9[s]);
		free_and_null((char **)&tau10[s]);
		free_and_null((char **)&tau13[s]);
		free_and_null((char **)&tau14[s]);
	}
	free_and_null((char **)&tau9);
	free_and_null((char **)&tau10);
	free_and_null((char **)&tau13);
	free_and_null((char **)&tau14);
	free_and_null((char **)&tau31);
	free_and_null((char **)&tau41);
	for (s = 0; s < S; s++) {
		free_and_null((char **)&tau51[s]);
		free_and_null((char **)&tau61[s]);
		free_and_null((char **)&tau71[s]);
		free_and_null((char **)&tau81[s]);
		free_and_null((char **)&tau111[s]);
		free_and_null((char **)&tau121[s]);
	}
	free_and_null((char **)&tau51);
	free_and_null((char **)&tau61);
	free_and_null((char **)&tau71);
	free_and_null((char **)&tau81);
	free_and_null((char **)&tau111);
	free_and_null((char **)&tau121);
	for (s = 0; s < S; s++) {
		for (i = 0; i < I; i++) {
			free_and_null((char **)&tau91[s][i]);
			free_and_null((char **)&tau101[s][i]);
			free_and_null((char **)&tau131[s][i]);
			free_and_null((char **)&tau141[s][i]);
		}
		free_and_null((char **)&tau91[s]);
		free_and_null((char **)&tau101[s]);
		free_and_null((char **)&tau131[s]);
		free_and_null((char **)&tau141[s]);
	}
	free_and_null((char **)&tau91);
	free_and_null((char **)&tau101);
	free_and_null((char **)&tau131);
	free_and_null((char **)&tau141);
	for (i = 0; i < I; i++) {
		free_and_null((char **)&y0[i]);
		free_and_null((char **)&y2[i]);
		free_and_null((char **)&eta1[i]);
		free_and_null((char **)&eta2[i]);
		free_and_null((char **)&w2[i]);
		free_and_null((char **)&w3[i]);
		
	}
	free_and_null((char **)&y0);
	free_and_null((char **)&y2);
	free_and_null((char **)&eta1);
	free_and_null((char **)&eta2);
	free_and_null((char **)&w2);
	free_and_null((char **)&w3);
	
	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			free_and_null((char **)&y1[i][j]);
			free_and_null((char **)&y11[i][j]);
			free_and_null((char **)&eta3[i][j]);
			free_and_null((char **)&eta4[i][j]);
			free_and_null((char **)&eta5[i][j]);
			free_and_null((char **)&eta6[i][j]);
		}
		free_and_null((char **)&y1[i]);
		free_and_null((char **)&y11[i]);
		free_and_null((char **)&eta3[i]);
		free_and_null((char **)&eta4[i]);
		free_and_null((char **)&eta5[i]);
		free_and_null((char **)&eta6[i]);
	}
	free_and_null((char **)&y1);
	free_and_null((char **)&y11);
	free_and_null((char **)&eta3);
	free_and_null((char **)&eta4);
	free_and_null((char **)&eta5);
	free_and_null((char **)&eta6);
	
	return obj_oos;
TERMINATE:
	return status;
}
