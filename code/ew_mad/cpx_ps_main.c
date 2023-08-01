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
#include "CreatIntVector.h"
#include "CreatIntMatrix.h"
#include "CreatIntTmatrix.h"
#include "CreatIntFmatrix.h"
#include "CreatDoubleVector.h"
#include "CreatDoubleMatrix.h"
#include "CreatDoubleTmatrix.h"
#include "OpenFile.h"
#include "Free.h"
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define EPSILON	 0.000000000001    //Added to the fractional solution obtained at the node
typedef struct node
{
	int data;
	struct node *next;
}linklist;
linklist *Ly1, *pl, *ql, *pl1;
linklist *Ly11, *pl11;
linklist *Leta1, *ple1;
linklist *Leta2, *ple2;
linklist *Leta3, *ple3;
linklist *Leta4, *ple4;
linklist *Leta5, *ple5;
linklist *Leta6, *ple6;
linklist *Ly2, *ply2;
linklist *Ly21, *ply21;

void release_list(linklist *pHead);
void xiangjian(linklist *head1, linklist *head2);
void hebing(linklist *head1, linklist *head2);
void print(linklist *head);
double sigma;

int main(int  argc, char *argv[]) {
	int		   num_inst, confirm;
	FILE       *ini;
	clock_t	   start_time1, end_time1;
	double     cputimeMIP1 = 0.;
	int        status = 0;
	char       inst_name[40];
	strcpy(inst_name, "input.txt");
	/*
	if (argc == 1) {
	printf("Error: Input file not specified \n");
	exit(8);
	}
	ini = Open_File(argv[1], "r");
	*/
	ini = Open_File(inst_name, "r");
	FILE		*out;
	time_t		tt;																														//Variables for taking the time and knowing when the code was exectued
	struct tm	*tm;
	FILE		*solution;
	tt = time(NULL);
	tm = localtime(&tt);
	char outfile[40];
	char xvalue[40];
	char Solution[40];
	char log1[40];

	fscanf(ini, "%d", &num_inst);
	fscanf(ini, "%s", &outfile);
	//	fscanf(ini, "%s", &Progress);
	fscanf(ini, "%s", &xvalue);
	fscanf(ini, "%s", &Solution);
	fscanf(ini, "%s", &log1);
	//	fscanf(ini, "%s", &Solution1);
	//fscanf(ini, "%s", &Para);
	printf("num_inst is %d", num_inst);
	/*Printing the time and date these are being executed*/
	/************************************************/
	out = Open_File(outfile, "a+");
	fprintf(out, "\n %s\n", asctime(tm));
	fclose(out);

	solution = Open_File(xvalue, "a+");
	fprintf(solution, "\n %s\n", asctime(tm));
	fclose(solution);

	solution = Open_File(Solution, "a+");
	fprintf(solution, "\n %s\n", asctime(tm));
	fclose(solution);

	solution = Open_File(log1, "a+");
	fprintf(solution, "\n %s\n", asctime(tm));
	fclose(solution);
	int r, I, S;
	double d;
	char	   instance[40];           // External file containing all the input data of the problem
	char	   instance1[40];           // External file containing all the input data of the problem
	int num_r = 0;
	double avg_obj = 0.;
	double avg_time = 0.;
	double avg_obj_iter = 0.;
	double avg_time_iter = 0.;
	double avg_iter = 0.;
	int num_fea = 0;
	double avg_col = 0.;
	int avg_col_iter = 0;
	double avg_sparse = 0.;
	double avg_sparse_iter = 0.;
	double avg_tcol = 0.;
	double avg_tcol_iter = 0.;
	for (r = 0; r < num_inst; r++) {
		avg_obj_iter = 0.;
		avg_time_iter = 0.;
		cputimeMIP1 = 0.;
		avg_sparse_iter = 0.;
		avg_tcol_iter = 0.;
		confirm = 0;
		avg_col_iter = 0;
		/***********************Now reviewing the information for the instance**************************/
		if (fscanf(ini, "%d", &I) == 1) confirm++;
		if (fscanf(ini, "%d", &S) == 1) confirm++;
		if (fscanf(ini, "%s", &instance) == 1) confirm++;//data name
		if (fscanf(ini, "%s", &instance1) == 1) confirm++;//data name
		if (fscanf(ini, "%lf", &sigma) == 1) confirm++;
		if (confirm < 5)  goto TERMINATE;
		out = Open_File(outfile, "a+");
		fprintf(out, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", r + 1, I, S, instance, inst_name, sigma);
		fclose(out);

		solution = Open_File(xvalue, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", r + 1, I, S, instance, inst_name, sigma);
		fclose(solution);

		solution = Open_File(Solution, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", r + 1, I, S, instance, inst_name, sigma);
		fclose(solution);

		solution = Open_File(log1, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", r + 1, I, S, instance, inst_name, sigma);
		fclose(solution);
		d = 50 * I;
		double *p, *c, *uz, **z;
		p = create_double_vector(I);
		c = create_double_vector(I);
		uz = create_double_vector(I);
		z = create_double_matrix(I, S);

		FILE *in;
		int index = 0;
		char       path[50];
		sprintf(path, "./nsd/");
		strcat(path, instance);
		in = Open_File(path, "r");
		int OK = 0;
		index = 0;
		int i, j, s;
		if (OK == 0) {
			for (i = 0; i < I; i++) {
				fscanf(in, "%lf", &p[i]);
				index++;
			}
			if (index < I) {
				fprintf(stderr, "ERROR: Cannot read fixed cost \n");
				OK = 2;
				printf("ok is %d", OK);
				goto TERMINATE;
			}
			index = 0;


			for (i = 0; i < I; i++) {
				fscanf(in, "%lf", &uz[i]);
				index++;
			}
			if (index < I) {
				fprintf(stderr, "ERROR: Cannot read fixed cost \n");
				OK = 2;
				printf("ok is %d", OK);
				goto TERMINATE;
			}
			index = 0;

			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					fscanf(in, "%lf", &z[i][s]);
					index++;
				}
			}
			if (index < I * S) {
				fprintf(stderr, "ERROR: Cannot read fixed cost \n");
				OK = 2;
				printf("ok is %d", OK);
				goto TERMINATE;
			}
			index = 0;

		}
		fclose(in);

		for (i = 0; i < I; i++) {
			c[i] = 1.;
		}
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
		Ly2 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
		Ly2->next = NULL;
		ply2 = Ly2;

		Ly21 = (linklist*)malloc(sizeof(linklist)); //具有头结点  
		Ly21->next = NULL;
		ply21 = Ly21;


		int *x, alpha, *r, **t, **y0, ***y1, ***y2, **rho, **rho1;
		int **eta1, **eta2, ***eta3, ***eta4, ***eta5, ***eta6;
		int **w1, ***w2, ***w3;
		double *tau;

		int tau1, tau2, *tau3, **tau4, **tau5, **tau6, **tau7, ***tau8, ***tau9, ***tau10, **tau11, ***tau12, ***tau13, ***tau14;
		tau3 = create_int_vector(S);
		tau4 = create_int_matrix(S, I);
		tau5 = create_int_matrix(S, I);
		tau6 = create_int_matrix(S, I);
		tau7 = create_int_matrix(S, I);
		tau8 = create_int_tmatrix(S, I, I);
		tau9 = create_int_tmatrix(S, I, I);
		tau10 = create_int_tmatrix(S, I, I);
		tau11 = create_int_matrix(S, I);
		tau12 = create_int_tmatrix(S, I, I);
		tau13 = create_int_tmatrix(S, I, I);
		tau14 = create_int_tmatrix(S, I, I);
		int tau21, *tau31, **tau41, **tau51, **tau61, **tau71, ***tau81, ***tau91, ***tau101, **tau111, ***tau121, ***tau131, ***tau141;
		tau31 = create_int_vector(S);
		tau41 = create_int_matrix(S, I);
		tau51 = create_int_matrix(S, I);
		tau61 = create_int_matrix(S, I);
		tau71 = create_int_matrix(S, I);
		tau81 = create_int_tmatrix(S, I, I);
		tau91 = create_int_tmatrix(S, I, I);
		tau101 = create_int_tmatrix(S, I, I);
		tau111 = create_int_matrix(S, I);
		tau121 = create_int_tmatrix(S, I, I);
		tau131 = create_int_tmatrix(S, I, I);
		tau141 = create_int_tmatrix(S, I, I);
		int indexd;
		x = create_int_vector(I);
		r = create_int_vector(S);
		y0 = create_int_matrix(I, S);
		t = create_int_matrix(I, S);
		rho = create_int_matrix(I, S);
		rho1 = create_int_matrix(I, S);
		t = create_int_matrix(I, S);
		y1 = create_int_tmatrix(I, I, S);
		y2 = create_int_tmatrix(I, I, S);
		eta1 = create_int_matrix(I, S);
		eta2 = create_int_matrix(I, S);
		eta3 = create_int_tmatrix(I, I, S);
		eta4 = create_int_tmatrix(I, I, S);
		eta5 = create_int_tmatrix(I, I, S);
		eta6 = create_int_tmatrix(I, I, S);
		w1 = create_int_matrix(I, S);
		w2 = create_int_tmatrix(I, I, S);
		w3 = create_int_tmatrix(I, I, S);
		int ***y11;
		y11 = create_int_tmatrix(I, I, S);

		CPXENVptr env = NULL;
		CPXLPptr  lp = NULL;
		int iter, Iter;
		iter = 0;
		Iter = 1;

		while (iter < Iter) {
			env = CPXopenCPLEX(&status);
			if (env == NULL) {
				fprintf(stderr, "Failure in CPXopenCPLEX, status = %d.\n", status);
				goto TERMINATE;
			}

			/* Turn on output to the screen */

			status = CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
			if (status) {
				fprintf(stderr, "Failed to turn on screen indicator, status = %d.\n",
					status);
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
				double cost = -1.*p[i];
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

			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					t[i][s] = num_r_var + num_x_var + num_alpha_var + index1;
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

			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					for (s = 0; s < S; s++) {
						y2[i][j][s] = num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
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
			int num_y2_var = index1;
			index1 = 0;

			double *lb_y2, *ub_y2;
			lb_y2 = create_double_vector(I*S);
			ub_y2 = create_double_vector(I*S);


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

			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					w1[i][s] = num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
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

			int num_w1_var = index1;
			index1 = 0;

			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					for (s = 0; s < S; s++) {
						w2[i][j][s] = num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
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
			}
			int num_w2_var = index1;
			index1 = 0;

			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					for (s = 0; s < S; s++) {
						w3[i][j][s] = num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
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
			int alpha1, *r1, **y01;
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

			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					rho[i][s] = num_y01_var + num_r1_var + num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
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
			int num_rho_var = index1;
			index1 = 0;

			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					rho1[i][s] = num_rho_var + num_y01_var + num_r1_var + num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
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
			int num_rho1_var = index1;
			index1 = 0;

			int num_tau1 = num_rho1_var + num_rho_var + num_y01_var + num_r1_var + num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;

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
				for (j = 0; j < I; j++) {
					tau4[s][j] = num_tau3 + num_tau2 + num_tau1 + indexd;
					tau41[s][j] = num_tau3 + num_tau2 + num_tau1 + indexd;
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
					for (j = 0; j < I; j++) {
						tau8[s][i][j] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
						tau81[s][i][j] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
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
					for (j = 0; j < I; j++) {
						tau12[s][i][j] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
						tau121[s][i][j] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
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
			//1
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

			sense = 'E';
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

			//2
			rmatind = (int *)malloc((S * 2 + 3 + I * S * 3) * sizeof(int));
			if (rmatind == NULL) {
				fprintf(stderr, "No memory for rmatind array.\n");
				status = -1;
				goto TERMINATE;
			}
			rmatval = (double *)malloc((S * 2 + 3 + I * S * 3) * sizeof(double));
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
				rmatval[nzcnt++] = -1. / S;
			}
			for (s = 0; s < S; s++) {
				rmatind[nzcnt] = r1[s];
				rmatval[nzcnt++] = 1. / S;
			}
			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					rmatind[nzcnt] = t[i][s];
					rmatval[nzcnt++] = -1. / S * sigma;
				}
			}
			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					rmatind[nzcnt] = rho[i][s];
					rmatval[nzcnt++] = -1. / S * z[i][s];
				}
			}
			for (i = 0; i < I; i++) {
				for (s = 0; s < S; s++) {
					rmatind[nzcnt] = rho1[i][s];
					rmatval[nzcnt++] = 1. / S * z[i][s];
				}
			}

			status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
				&rmatbeg, rmatind, rmatval, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
				goto TERMINATE;
			}
			free_and_null((char **)&rmatind);
			free_and_null((char **)&rmatval);
			//3
			rmatind = (int *)malloc((I*(4 + I) + 3 + I * I + I * 2) * sizeof(int));
			if (rmatind == NULL) {
				fprintf(stderr, "No memory for rmatind array.\n");
				status = -1;
				goto TERMINATE;
			}
			rmatval = (double *)malloc((I*(4 + I) + 3 + I * I + I * 2) * sizeof(double));
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
					rmatval[nzcnt++] = z[i][s];
				}
				for (i = 0; i < I; i++) {
					rmatind[nzcnt] = rho[i][s];
					rmatval[nzcnt++] = -1.*z[i][s];
				}
				for (i = 0; i < I; i++) {
					rmatind[nzcnt] = rho1[i][s];
					rmatval[nzcnt++] = z[i][s];
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
			//4
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
				for (j = 0; j < I; j++) {
					nzcnt = 0;
					rmatind[nzcnt] = tau4[s][j];
					rmatval[nzcnt++] = 1.;

					rmatind[nzcnt] = w1[j][s];
					rmatval[nzcnt++] = 1.;
					for (i = 0; i < I; i++) {
						rmatind[nzcnt] = y2[i][j][s];
						rmatval[nzcnt++] = p[i];
					}

					rmatind[nzcnt] = t[j][s];
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
			//5
			rmatind = (int *)malloc((I * 2 + 6) * sizeof(int));
			if (rmatind == NULL) {
				fprintf(stderr, "No memory for rmatind array.\n");
				status = -1;
				goto TERMINATE;
			}
			rmatval = (double *)malloc((I * 2 + 6) * sizeof(double));
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

					rmatind[nzcnt] = w1[j][s];
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
					rmatind[nzcnt] = rho[j][s];
					rmatval[nzcnt++] = -1.;
					rmatind[nzcnt] = rho1[j][s];
					rmatval[nzcnt++] = 1.;
					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}
				}
			}
			//6

			sense = 'E';
			rhs = 0.;
			for (s = 0; s < S; s++) {
				for (j = 0; j < I; j++) {
					nzcnt = 0;
					rmatind[nzcnt] = tau6[s][j];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = w1[j][s];
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
					rmatind[nzcnt] = rho[j][s];
					rmatval[nzcnt++] = 1.;
					rmatind[nzcnt] = rho1[j][s];
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
			//7
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
						rmatval[nzcnt++] = z[j][s];
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
			//8
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
					for (j = 0; j < I; j++) {
						nzcnt = 0;
						rmatind[nzcnt] = tau8[s][i][j];
						rmatval[nzcnt++] = 1.;
						rmatind[nzcnt] = y2[i][j][s];
						rmatval[nzcnt++] = -1.;

						rmatind[nzcnt] = w2[i][j][s];
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
			//9
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
						rmatind[nzcnt] = w2[i][j][s];
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
			//10			

			sense = 'E';
			rhs = 0.;
			for (s = 0; s < S; s++) {
				for (i = 0; i < I; i++) {
					for (j = 0; j < I; j++) {
						nzcnt = 0;
						rmatind[nzcnt] = tau10[s][i][j];
						rmatval[nzcnt++] = 1.;
						rmatind[nzcnt] = w2[i][j][s];
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
			//11
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
						rmatval[nzcnt++] = z[j][s];
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
			//12
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
					for (j = 0; j < I; j++) {
						nzcnt = 0;
						rmatind[nzcnt] = tau12[s][i][j];
						rmatval[nzcnt++] = 1.;
						rmatind[nzcnt] = y2[i][j][s];
						rmatval[nzcnt++] = -1.;
						rmatind[nzcnt] = w3[i][j][s];
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
			//13
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

						rmatind[nzcnt] = w3[i][j][s];
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
			//14
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
						rmatind[nzcnt] = w3[i][j][s];
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
			CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);
			//	CPXsetdblparam(env, CPXPARAM_Barrier_ConvergeTol, 0.01); // e-optimal solution (%gap)

			start_time1 = clock();//CPXprimopt
			status = CPXlpopt(env, lp);
			if (status) {
				fprintf(stderr, "Error in CPXprimopt: status = %d\n", status);
				goto TERMINATE;
			}
			int worker_lp_sol_stat = CPXgetstat(env, lp);
			int num_sub = CPXgetnumcols(env, lp);
			printf("status is %d\t", worker_lp_sol_stat);
			if (worker_lp_sol_stat == CPX_STAT_OPTIMAL) {
				end_time1 = clock();
				cputimeMIP1 = (double)(end_time1 - start_time1) / (double)CLOCKS_PER_SEC;

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
				solution = Open_File(xvalue, "a+");

				for (j = 0; j < I; j++) {
					fprintf(solution, "%.4lf\t", x2[x[j]]);
				}
				fprintf(solution, "\n");
				fclose(solution);
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
					for (j = 0; j < I; j++) {
						tau4[s][j] = num_tau3 + num_tau2 + num_tau1 + indexd;
						indexd++;
					}
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
						for (j = 0; j < I; j++) {
							tau8[s][i][j] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							indexd++;
						}
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
						for (j = 0; j < I; j++) {
							tau12[s][i][j] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							indexd++;
						}
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

				for (i = 0; i < num_sub; i++) {
					if (x2[i] > 0.00001) {
						avg_sparse_iter = avg_sparse_iter + 1.;
					}
				}
				avg_tcol_iter = num_sub;
				free_and_null((char **)&x2);
				free_and_null((char **)&tau);
				avg_obj_iter = z1;
				avg_time_iter += cputimeMIP1;

				solution = Open_File(outfile, "a+");
				fprintf(solution, "%.4lf\t%.4lf\t%.0lf\t%d\n", avg_obj_iter, avg_time_iter, avg_sparse_iter, avg_col_iter);
				fclose(solution);

				//	print(Ly1);
			}
			else {
				double best_upper_bound, best_lower_bound;
				CPXgetobjval(env, lp, &best_upper_bound);
				CPXgetbestobjval(env, lp, &best_lower_bound);
				solution = Open_File(outfile, "a+");
				fprintf(solution, "%.4lf\t%.4lf\t%.2lf\t", best_upper_bound, best_lower_bound, ABS(best_upper_bound - best_lower_bound) / best_upper_bound * 100);
				fclose(solution);
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


		if (num_r < 50) {
			if (avg_time_iter < 7200) {
				avg_time += avg_time_iter;
				num_fea++;
				avg_obj += avg_obj_iter;
				avg_sparse += avg_sparse_iter;
				avg_tcol += avg_tcol_iter;
			}
			avg_iter += iter;
			avg_col += avg_col_iter;
			if (num_r == 49) {
				solution = Open_File(outfile, "a+");
				fprintf(solution, "\n");
				fprintf(solution, "%.4lf\t %.4lf\t %.4lf\t %.1lf\t %d\t %.0lf\t %.0lf\n", avg_obj / num_fea, avg_time / num_fea, avg_iter / 50., avg_col / 50., num_fea, avg_sparse*1. / num_fea, avg_tcol*1. / num_fea);
				fclose(solution);

			}
			num_r++;

		}
		if (num_r == 10) {
			num_r = 0;
			avg_obj = 0.;
			avg_time = 0.;
			avg_iter = 0.;
			num_fea = 0;
			avg_col = 0.;
			avg_sparse = 0.;
			avg_tcol = 0.;
		}
		/* Init the CPLEX environment */
		release_list(Ly1);
		release_list(Ly11);
		release_list(Ly2);
		release_list(Ly21);
		release_list(Leta1);
		release_list(Leta2);
		release_list(Leta3);
		release_list(Leta4);
		release_list(Leta5);
		release_list(Leta6);



		free_and_null((char **)&x);
		free_and_null((char **)&r);
		for (i = 0; i < I; i++) {
			free_and_null((char **)&w1[i]);
		}
		free_and_null((char **)&w1);
		free_and_null((char **)&p);
		free_and_null((char **)&c);
		free_and_null((char **)&uz);
		free_and_null((char **)&tau3);
		for (s = 0; s < S; s++) {
			free_and_null((char **)&tau4[s]);
			free_and_null((char **)&tau5[s]);
			free_and_null((char **)&tau6[s]);
			free_and_null((char **)&tau7[s]);
			for (i = 0; i < I; i++) {
				free_and_null((char **)&tau8[s][i]);
				free_and_null((char **)&tau12[s][i]);

			}
			free_and_null((char **)&tau8[s]);
			free_and_null((char **)&tau12[s]);
			free_and_null((char **)&tau11[s]);

		}
		free_and_null((char **)&tau4);
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
		for (s = 0; s < S; s++) {
			free_and_null((char **)&tau41[s]);
			free_and_null((char **)&tau51[s]);
			free_and_null((char **)&tau61[s]);
			free_and_null((char **)&tau71[s]);
			for (i = 0; i < I; i++) {
				free_and_null((char **)&tau81[s][i]);
				free_and_null((char **)&tau121[s][i]);

			}
			free_and_null((char **)&tau81[s]);
			free_and_null((char **)&tau121[s]);
			free_and_null((char **)&tau111[s]);

		}
		free_and_null((char **)&tau41);
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
		free_and_null((char **)&tau9);
		free_and_null((char **)&tau10);
		free_and_null((char **)&tau13);
		free_and_null((char **)&tau14);

		for (i = 0; i < I; i++) {
			free_and_null((char **)&y0[i]);
			for (j = 0; j < I; j++) {
				free_and_null((char **)&y2[i][j]);
				free_and_null((char **)&w2[i][j]);
				free_and_null((char **)&w3[i][j]);

			}
			free_and_null((char **)&y2[i]);
			free_and_null((char **)&t[i]);
			free_and_null((char **)&rho[i]);
			free_and_null((char **)&rho1[i]);
			free_and_null((char **)&eta1[i]);
			free_and_null((char **)&eta2[i]);
			free_and_null((char **)&w2[i]);
			free_and_null((char **)&w3[i]);
			free_and_null((char **)&z[i]);
		}
		free_and_null((char **)&y0);
		free_and_null((char **)&y2);
		free_and_null((char **)&t);
		free_and_null((char **)&rho);
		free_and_null((char **)&rho1);
		free_and_null((char **)&eta1);
		free_and_null((char **)&eta2);
		free_and_null((char **)&w2);
		free_and_null((char **)&w3);
		free_and_null((char **)&z);
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

	}
TERMINATE:
	return status;

}
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
void xiangjian(linklist *head1, linklist *head2)
{
	int flag = 0;
	linklist *p1, *p2, *q;
	p1 = head1;
	p2 = head1->next;
	while (p2)
	{
		for (q = head2->next; q; q = q->next)
			if (p2->data == q->data)
			{
				p1->next = p2->next;
				free(p2);
				p2 = p1->next;
				flag = 1;
				break;
			}
		if (flag == 0)
		{
			p1 = p1->next;
			p2 = p2->next;
		}
		else
			flag = 0;
	}

}
void print(linklist *head)
{
	linklist *p;
	p = head->next;
	while (p)
	{

		printf("%d  ", p->data);
		p = p->next;
	}
}
void hebing(linklist *head1, linklist *head2)
{
	int flag = 0;
	linklist *p, *q, *m;
	p = head1->next;
	while (p)
	{
		for (q = head2->next; q; q = q->next)
		{
			if (p->data == q->data)
			{
				flag = 1;
				break;
			}

		}
		if (flag == 0)
		{
			m = (linklist*)malloc(sizeof(linklist));
			m->data = p->data;
			m->next = head2->next;
			head2->next = m;
		}
		flag = 0;
		p = p->next;
	}

}

