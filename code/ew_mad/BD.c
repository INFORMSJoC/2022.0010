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
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define EPSILON	 0.000000000001    //Added to the fractional solution obtained at the node
double *create_double_vector(int);
double **create_double_matrix(int, int);
double ***create_double_tmatrix(int, int, int);
int    *create_int_vector(int);
int    **create_int_matrix(int, int);
int    ***create_int_tmatrix(int, int, int);
int    ****create_int_fmatrix(int, int, int, int);
static void free_and_null(char **ptr);
FILE   *Open_File(const char *name, const char *mode);
int *x, *r, **t, **y01, ***y11, ***y21, *alpha, **rho;
int **eta1, **eta2, ***eta3, ***eta4, ***eta5, ***eta6;
int **w1, ***w2, ***w3;
double ***tau;
double *p, *c, *uz, **z, sigma;
int r1, I, S;
double theta, d;
int tau1, tau2;
int tau3, *tau4, *tau5, *tau6, *tau7, **tau8, **tau9, **tau10, *tau11, **tau12, **tau13, **tau14;
int iter, Iter;
double UB, LB;
double lambda;
int compute_subproblem(double *obj, double *x1, int s);
int solve_projection(double obj, double *x1, double *xk);
int main(int  argc, char *argv[]) {
	int		   num_inst, confirm;
	FILE       *ini;
	clock_t	   start_time1, end_time1;
	double     cputimeMIP1 = 0.;
	clock_t	   start_time2, end_time2;
	double     cputimeMIP2 = 0.;
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
	char	   instance[40];           // External file containing all the input data of the problem
	char	   instance1[40];           // External file containing all the input data of the problem
	int num_r = 0;
	double avg_obj = 0.;
	double avg_time = 0.;
	double avg_time1 = 0.;
	double avg_obj_iter = 0.;
	double avg_time_iter = 0.;
	double avg_time1_iter = 0.;
	double prop_out = 0.;
	double avg_iter = 0.;
	lambda = 0.2;
	int num_fea = 0;
	for (r1 = 0; r1 < num_inst; r1++) {
		avg_obj_iter = 0.;
		avg_time_iter = 0.;
		avg_time1_iter = 0.;

		cputimeMIP1 = 0.;
		confirm = 0;
		LB = -10000.;
		UB = 10000.;
		/***********************Now reviewing the information for the instance**************************/
		if (fscanf(ini, "%d", &I) == 1) confirm++;
		if (fscanf(ini, "%d", &S) == 1) confirm++;
		if (fscanf(ini, "%s", &instance) == 1) confirm++;//data name
		if (fscanf(ini, "%s", &instance1) == 1) confirm++;//data name
		if (fscanf(ini, "%lf", &sigma) == 1) confirm++;

		//if (fscanf(ini, "%lf", &lambda) == 1) confirm++;
														  //	if (fscanf(ini, "%d", &S1) == 1) confirm++;
														 
		if (confirm < 5)  goto TERMINATE;
		out = Open_File(outfile, "a+");
		fprintf(out, "%d;\t %d\t %d\t %s\t%s\n", r1 + 1, I, S, instance, inst_name);
		fclose(out);

		solution = Open_File(xvalue, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\n", r1 + 1, I, S, instance, inst_name);
		fclose(solution);

		solution = Open_File(Solution, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\n", r1 + 1, I, S, instance, inst_name);
		fclose(solution);

		solution = Open_File(log1, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\n", r1 + 1, I, S, instance, inst_name);
		fclose(solution);
		d = 50 * I;
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
		for (i = 0; i < I; i++) {
			c[i] = 1.;
		}
		

		x = create_int_vector(I);
		r = create_int_vector(S);
		t = create_int_matrix(I, S);
		y01 = create_int_matrix(I, S);
		y11 = create_int_tmatrix(I, I, S);
		y21 = create_int_tmatrix(I, I, S);
		eta1 = create_int_matrix(I, S);
		eta2 = create_int_matrix(I, S);
		eta3 = create_int_tmatrix(I, I, S);
		eta4 = create_int_tmatrix(I, I, S);
		eta5 = create_int_tmatrix(I, I, S);
		eta6 = create_int_tmatrix(I, I, S);
		w1 = create_int_matrix(I, S);
		w2 = create_int_tmatrix(I, I, S);
		w3 = create_int_tmatrix(I, I, S);
		alpha = create_int_vector(S);
		rho = create_int_matrix(I, S);
		//need to modify
		int indexd = 0;
		tau4 = create_int_vector(I);
		tau5 = create_int_vector(I);
		tau6 = create_int_vector(I);
		tau7 = create_int_vector(I);
		tau8 = create_int_matrix(I, I);
		tau9 = create_int_matrix(I, I);
		tau10 = create_int_matrix(I, I);
		tau11 = create_int_vector(I);
		tau12 = create_int_matrix(I, I);
		tau13 = create_int_matrix(I, I);
		tau14 = create_int_matrix(I, I);

		tau3 = indexd;
		indexd++;
		int num_tau3 = indexd;
		indexd = 0;

		for (j = 0; j < I; j++) {
			tau4[j] = num_tau3 + indexd;
			indexd++;
		}
		int num_tau4 = indexd;
		indexd = 0;


		for (j = 0; j < I; j++) {
			tau5[j] = num_tau4 + num_tau3 + indexd;
			indexd++;
		}
		int num_tau5 = indexd;
		indexd = 0;

		for (j = 0; j < I; j++) {
			tau6[j] = num_tau5 + num_tau4 + num_tau3 + indexd;
			indexd++;
		}
		int num_tau6 = indexd;
		indexd = 0;

		for (i = 0; i < I; i++) {
			tau7[i] = num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
			indexd++;
		}
		int num_tau7 = indexd;
		indexd = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				tau8[i][j] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
				indexd++;
			}
		}
		int num_tau8 = indexd;
		indexd = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				tau9[i][j] = num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
				indexd++;
			}
		}
		int num_tau9 = indexd;
		indexd = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				tau10[i][j] = num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
				indexd++;
			}
		}
		int num_tau10 = indexd;
		indexd = 0;

		for (i = 0; i < I; i++) {
			tau11[i] = num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
			indexd++;
		}
		int num_tau11 = indexd;
		indexd = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				tau12[i][j] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
				indexd++;
			}
		}
		int num_tau12 = indexd;
		indexd = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				tau13[i][j] = num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
				indexd++;
			}
		}
		int num_tau13 = indexd;
		indexd = 0;


		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				tau14[i][j] = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
				indexd++;
			}
		}
		int num_tau14 = indexd;
		indexd = 0;


		int num_tau = num_tau14 + num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + indexd;
		Iter = 10000;
		iter = 0;
		tau = create_double_tmatrix(Iter, S, num_tau);
		int flag_solve = 0;

		while (avg_time1_iter < 7200) {
			CPXENVptr env = NULL;
			CPXLPptr  lp = NULL;
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

			for (s = 0; s < S; s++) {
				alpha[s] = num_x_var + index1;
				index1++;
				double lb = -100000.;
				double ub = CPX_INFBOUND;
				double cost = 1.;
				status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
				if (status) {
					fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
					goto TERMINATE;
				}
			}
			int num_alpha_var = index1;
			index1 = 0;

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

			int k;
			rmatind = (int *)malloc((I + 1) * sizeof(int));
			if (rmatind == NULL) {
				fprintf(stderr, "No memory for rmatind array.\n");
				status = -1;
				goto TERMINATE;
			}
			rmatval = (double *)malloc((I + 1) * sizeof(double));
			if (rmatval == NULL) {
				fprintf(stderr, "No memory for rmatval array.\n");
				status = -1;
				goto TERMINATE;
			}
			sense = 'L';
			for (s = 0; s < S; s++) {
				for (k = 0; k < iter; k++) {
					rhs = 0.;
					for (i = 0; i < I; i++) {
						rhs += tau[k][s][tau13[i][i]];
						rhs -= tau[k][s][tau14[i][i]];
						rhs += tau[k][s][tau11[i]] * z[i][s];
					}

					nzcnt = 0;
					rmatind[nzcnt] = alpha[s];
					rmatval[nzcnt++] = -1.;
					for (i = 0; i < I; i++) {
						rmatind[nzcnt] = x[i];
						rmatval[nzcnt++] = tau[k][s][tau11[i]];
					}
					for (i = 0; i < I; i++) {
						if (ABS(rhs) > 0.000001 || ABS(tau[k][s][tau11[i]]) > 0.000001) {
							status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
								&rmatbeg, rmatind, rmatval, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
								goto TERMINATE;
							}
							break;
						}
					}

					
				}
			}
			free_and_null((char **)&rmatind);
			free_and_null((char **)&rmatval);

			CPXsetdblparam(env, CPX_PARAM_TILIM, 7200);
			CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
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
				status = CPXgetx(env, lp, x2, 0, num_sub - 1);
				if (status) {
					fprintf(stderr, "Failed to obtain solution, status = %d.\n", status);
					goto TERMINATE;
				}
				double *xk;
				xk = create_double_vector(I);
				for (i = 0; i < I; i++) {
					xk[i] = x2[x[i]];
				}
				start_time2 = clock();//CPXprimopt

				LB = z1;

				if (ABS((UB - LB) / UB) <= 0.0001 || UB < LB) {
					flag_solve = 1;
				}
				if (flag_solve == 0) {
					
					double sum_sub = 0.;
					double obj_s;

					for (s = 0; s < S; s++) {
						compute_subproblem(&obj_s, xk, s);
						sum_sub += obj_s;
						sum_sub -= x2[alpha[s]];
					}
					sum_sub += z1;
					if (UB > sum_sub) {
						UB = sum_sub;
					}
				}
				end_time2 = clock();
				cputimeMIP2 = (double)(end_time2 - start_time2) / (double)CLOCKS_PER_SEC;

				
				avg_time1_iter += cputimeMIP1 + cputimeMIP2;

				solution = Open_File(log1, "a+");
				fprintf(solution, "%.4lf\t%.4lf\t%.4lf\t%.4lf\n", LB, UB, cputimeMIP1, cputimeMIP2);
				fclose(solution);
				solution = Open_File(xvalue, "a+");
				for (i = 0; i < I; i++) {
					fprintf(solution, "%.2lf\t", x2[x[i]]);
				}
				fprintf(solution, "\n");
				fclose(solution);
				free_and_null((char **)&x2);
			}
			else {
				avg_time1_iter += 7200;
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
			if (flag_solve == 1) {
				break;
			}
			iter++;
		}
		solution = Open_File(outfile, "a+");
		fprintf(solution, "%.4lf\t%.4lf\t%.4lf\t%d\n", LB, UB, avg_time1_iter, iter);
		fclose(solution);
		if (num_r < 10) {
			avg_obj += LB;
			if (avg_time1_iter < 7200) {
				avg_time += avg_time1_iter;
				num_fea++;
			}
				
			avg_iter += iter;
			if (num_r == 9) {
				solution = Open_File(outfile, "a+");
				fprintf(solution, "\n");
				fprintf(solution, "%.4lf\t%.4lf\t%.1lf\n", avg_obj / 10., avg_time / num_fea, avg_iter / 10.);
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
		}
		free_and_null((char **)&x);
		free_and_null((char **)&r);
		free_and_null((char **)&p);
		free_and_null((char **)&c);
		free_and_null((char **)&uz);
		free_and_null((char **)&tau4);
		free_and_null((char **)&tau5);
		free_and_null((char **)&tau6);
		free_and_null((char **)&tau7);
		free_and_null((char **)&tau11);


		for (s = 0; s < I; s++) {
			free_and_null((char **)&tau8[s]);
			free_and_null((char **)&tau9[s]);
			free_and_null((char **)&tau10[s]);
			free_and_null((char **)&tau12[s]);
			free_and_null((char **)&tau13[s]);
			free_and_null((char **)&tau14[s]);
		}
		free_and_null((char **)&tau8);
		free_and_null((char **)&tau12);
		free_and_null((char **)&tau9);
		free_and_null((char **)&tau10);
		free_and_null((char **)&tau13);
		free_and_null((char **)&tau14);

		for (i = 0; i < I; i++) {
			free_and_null((char **)&y01[i]);
			free_and_null((char **)&eta1[i]);
			free_and_null((char **)&eta2[i]);
			free_and_null((char **)&z[i]);
			free_and_null((char **)&w1[i]);
			free_and_null((char **)&t[i]);
			free_and_null((char **)&rho[i]);
		}
		free_and_null((char **)&y01);
		free_and_null((char **)&eta1);
		free_and_null((char **)&eta2);
		free_and_null((char **)&z);
		free_and_null((char **)&w1);
		free_and_null((char **)&t);
		free_and_null((char **)&rho);

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				free_and_null((char **)&y11[i][j]);
				free_and_null((char **)&y21[i][j]);
				free_and_null((char **)&eta3[i][j]);
				free_and_null((char **)&eta4[i][j]);
				free_and_null((char **)&eta5[i][j]);
				free_and_null((char **)&eta6[i][j]);
				free_and_null((char **)&w2[i][j]);
				free_and_null((char **)&w3[i][j]);

			}
			free_and_null((char **)&y11[i]);
			free_and_null((char **)&y21[i]);
			free_and_null((char **)&eta3[i]);
			free_and_null((char **)&eta4[i]);
			free_and_null((char **)&eta5[i]);
			free_and_null((char **)&eta6[i]);
			free_and_null((char **)&w2[i]);
			free_and_null((char **)&w3[i]);

		}
		free_and_null((char **)&y11);
		free_and_null((char **)&y21);
		free_and_null((char **)&eta3);
		free_and_null((char **)&eta4);
		free_and_null((char **)&eta5);
		free_and_null((char **)&eta6);
		free_and_null((char **)&w2);
		free_and_null((char **)&w3);


		for (i = 0; i < Iter; i++) {
			for (j = 0; j < S; j++) {
				free_and_null((char **)&tau[i][j]);
			}
			free_and_null((char **)&tau[i]);
		}
		free_and_null((char **)&tau);
	}
TERMINATE:
	return status;
}
int compute_subproblem(double *obj, double *x1, int s) {
	int i, j;
	int status = 0;
	int worker_lp_sol_stat = 0;
	CPXENVptr env = NULL;
	CPXLPptr  lp = NULL;
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		fprintf(stderr, "Failure in CPXopenCPLEX, status = %d.\n", status);
		goto TERMINATE;
	}

	/* Turn on output to the screen

	status = CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
	if (status) {
	fprintf(stderr, "Failed to turn on screen indicator, status = %d.\n",
	status);
	goto TERMINATE;
	}*/
	//chage to max problem
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
	r[s] = index1;
	index1++;
	double lb = -CPX_INFBOUND;
	double ub = CPX_INFBOUND;
	double cost = 1. / S;
	status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
	if (status) {
		fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
		goto TERMINATE;
	}

	int num_r_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		y01[i][s] = num_r_var + index1;
		index1++;
		lb = -CPX_INFBOUND;
		ub = CPX_INFBOUND;
		status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	int num_y0_var = index1;
	index1 = 0;


	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			y11[i][j][s] = num_y0_var + num_r_var + index1;
			index1++;
			lb = -CPX_INFBOUND;
			ub = CPX_INFBOUND;
			status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
			}

		}
	}
	int num_y1_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			y21[i][j][s] = num_y1_var + num_y0_var + num_r_var + index1;
			index1++;
			lb = -CPX_INFBOUND;
			ub = CPX_INFBOUND;
			status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
			if (status) {
				fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
				goto TERMINATE;
			}
		}
	}
	int num_y2_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		eta1[i][s] = num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
		index1++;
		lb = 0.;
		ub = CPX_INFBOUND;
		status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	int num_eta1_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		eta2[i][s] = num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
		index1++;
		lb = 0.;
		ub = CPX_INFBOUND;
		status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	int num_eta2_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			eta3[i][j][s] = num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
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
	int num_eta3_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			eta4[i][j][s] = num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
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
	int num_eta4_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			eta5[i][j][s] = num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
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
	int num_eta5_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			eta6[i][j][s] = num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
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
	int num_eta6_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		w1[i][s] = num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
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
		for (j = 0; j < I; j++) {
			w2[i][j][s] = num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
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
		for (j = 0; j < I; j++) {
			w3[i][j][s] = num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
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
		t[i][s] = num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
		index1++;
		lb = 0.;
		ub = CPX_INFBOUND;
		cost = 1. / S * sigma;
		status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	int num_t_var = index1;
	index1 = 0;

	for (i = 0; i < I; i++) {
		rho[i][s] = num_t_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_r_var + index1;
		index1++;
		lb = -CPX_INFBOUND;
		ub = CPX_INFBOUND;
		cost = 1. / S * z[i][s];
		status = CPXnewcols(env, lp, 1, &cost, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	int num_rho_var = index1;
	index1 = 0;

	rmatbeg = 0;


	//3
	rmatind = (int *)malloc((1 + I * I + I * 4) * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc((1 + I * I + I * 4) * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	sense = 'L';
	rhs = 0.;

	nzcnt = 0;

	rmatind[nzcnt] = r[s];
	rmatval[nzcnt++] = -1.;
	for (i = 0; i < I; i++) {
		rmatind[nzcnt] = y01[i][s];
		rmatval[nzcnt++] = p[i];
	}

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = p[i] * z[j][s];
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

	status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
		&rmatbeg, rmatind, rmatval, NULL, NULL);
	if (status) {
		fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
		goto TERMINATE;
	}

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//4
	rmatind = (int *)malloc((I + 2) * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc((I + 2) * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	rhs = 0.;

	for (j = 0; j < I; j++) {
		nzcnt = 0;

		rmatind[nzcnt] = w1[j][s];
		rmatval[nzcnt++] = 1.;
		for (i = 0; i < I; i++) {
			rmatind[nzcnt] = y21[i][j][s];
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


	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//5
	rmatind = (int *)malloc((I + 4) * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc((I + 4) * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	rhs = 0.;

	for (j = 0; j < I; j++) {
		nzcnt = 0;

		rmatind[nzcnt] = w1[j][s];
		rmatval[nzcnt++] = -1.;
		for (i = 0; i < I; i++) {
			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = p[i];
		}

		rmatind[nzcnt] = eta1[j][s];
		rmatval[nzcnt++] = -1.;
		rmatind[nzcnt] = eta2[j][s];
		rmatval[nzcnt++] = 1.;
		rmatind[nzcnt] = rho[j][s];
		rmatval[nzcnt++] = -1.;

		status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
			&rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
			goto TERMINATE;
		}
	}

	//6

	rhs = 0.;

	for (j = 0; j < I; j++) {
		nzcnt = 0;
		rmatind[nzcnt] = w1[j][s];
		rmatval[nzcnt++] = -1.;
		for (i = 0; i < I; i++) {
			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = -1.*p[i];
		}

		rmatind[nzcnt] = eta1[j][s];
		rmatval[nzcnt++] = 1.;
		rmatind[nzcnt] = eta2[j][s];
		rmatval[nzcnt++] = -1.;
		rmatind[nzcnt] = rho[j][s];
		rmatval[nzcnt++] = 1.;

		status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
			&rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
			goto TERMINATE;
		}
	}

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//7
	rmatind = (int *)malloc((I * 3 + 1) * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc((I * 3 + 1) * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	rhs = 0.;

	for (i = 0; i < I; i++) {
		nzcnt = 0;

		rmatind[nzcnt] = y01[i][s];
		rmatval[nzcnt++] = -1.;
		for (j = 0; j < I; j++) {
			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = -1.*z[j][s];
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

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//8
	rmatind = (int *)malloc(2 * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc(2 * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	rhs = 0.;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			nzcnt = 0;
			rmatind[nzcnt] = y21[i][j][s];
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

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//9
	rmatind = (int *)malloc(4 * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc(4 * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	rhs = 0.;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			nzcnt = 0;
			rmatind[nzcnt] = w2[i][j][s];
			rmatval[nzcnt++] = -1.;

			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = -1.;
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

	//10			

	rhs = 0.;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			nzcnt = 0;
			rmatind[nzcnt] = w2[i][j][s];
			rmatval[nzcnt++] = -1.;

			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = 1.;
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

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//11
	rmatind = (int *)malloc((I * 3 + 1) * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc((I * 3 + 1) * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	for (i = 0; i < I; i++) {
		rhs = z[i][s] - x1[x[i]];
		nzcnt = 0;
		rmatind[nzcnt] = y01[i][s];
		rmatval[nzcnt++] = -1.;

		for (j = 0; j < I; j++) {
			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = -1.*z[j][s];
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

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//12
	rmatind = (int *)malloc(2 * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc(2 * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}

	rhs = 0.;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			nzcnt = 0;
			rmatind[nzcnt] = y21[i][j][s];
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

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	//13
	rmatind = (int *)malloc(4 * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc(4 * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}


	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			nzcnt = 0;
			if (i == j) {
				rhs = 1.;
			}
			else {
				rhs = 0.;
			}

			rmatind[nzcnt] = w3[i][j][s];
			rmatval[nzcnt++] = -1.;

			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = -1.;
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

	//rhs = -1.;

	for (i = 0; i < I; i++) {
		for (j = 0; j < I; j++) {
			nzcnt = 0;
			if (i == j) {
				rhs = -1.;
			}
			else {
				rhs = 0.;
			}

			rmatind[nzcnt] = w3[i][j][s];
			rmatval[nzcnt++] = -1.;

			rmatind[nzcnt] = y11[i][j][s];
			rmatval[nzcnt++] = 1.;
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

	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);

	CPXsetdblparam(env, CPX_PARAM_TILIM, 7200);
	CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
	status = CPXlpopt(env, lp);
	if (status) {
		fprintf(stderr, "Error in CPXprimopt: status = %d\n", status);
		goto TERMINATE;
	}
	worker_lp_sol_stat = CPXgetstat(env, lp);
	int num_sub = CPXgetnumcols(env, lp);
	//printf("status is %d\t", worker_lp_sol_stat);
	if (worker_lp_sol_stat == CPX_STAT_OPTIMAL) {
		status = CPXgetobjval(env, lp, obj);
		//printf("obj is %.2lf\t", *obj);
		if (status) {
			fprintf(stderr, "Failed to obtain solution, status = %d.\n", status);
			goto TERMINATE;
		}
		int num_row = CPXgetnumrows(env, lp);
		double *tauk;
		tauk = create_double_vector(num_row);
		status = CPXgetpi(env, lp, tauk, 0, num_row - 1);
		for (i = 0; i < num_row; i++) {
			tau[iter][s][i] = -1.*tauk[i];
		}
		free_and_null((char **)&tauk);

	}
	else {
		printf("The subproblem is infeasible or unbouded");
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
TERMINATE:
	return worker_lp_sol_stat;
}

int solve_projection(double obj, double *x1, double *xk) {
	int status = 0;
	//while (avg_time1_iter < 7200 && ABS((UB - LB) / UB)>0.0001) {
	CPXENVptr env1 = NULL;
	CPXLPptr  lp1 = NULL;
	env1 = CPXopenCPLEX(&status);
	if (env1 == NULL) {
		fprintf(stderr, "Failure in CPXopenCPLEX, status = %d.\n", status);
		goto TERMINATE;
	}

	/* Turn on output to the screen */

	status = CPXsetintparam(env1, CPXPARAM_ScreenOutput, CPX_ON);
	if (status) {
		fprintf(stderr, "Failed to turn on screen indicator, status = %d.\n",
			status);
		goto TERMINATE;
	}
	lp1 = CPXcreateprob(env1, &status, "master_ILP.lp");
	if (lp1 == NULL) {
		fprintf(stderr, "Failure in CPXcreateprob, status = %d.\n", status);
		goto TERMINATE;
	}
	char sense;
	int nzcnt, rmatbeg, *rmatind = NULL;
	double rhs, *rmatval = NULL;
	int index1 = 0;
	int i, s;

	for (i = 0; i < I; i++) {
		x[i] = index1;
		index1++;
		double lb = 0.;
		double ub = CPX_INFBOUND;
		double cost = -1.*p[i];
		status = CPXnewcols(env1, lp1, 1, NULL, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	int num_x_var = index1;
	index1 = 0;

	for (s = 0; s < S; s++) {
		alpha[s] = num_x_var + index1;
		index1++;
		double lb = -100000.;
		double ub = CPX_INFBOUND;
		double cost = 1.;
		status = CPXnewcols(env1, lp1, 1, NULL, &lb, &ub, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	int num_alpha_var = index1;
	index1 = 0;

	int rho;
	rho = num_x_var + num_alpha_var + index1;
	index1++;
	double lb = 0.;
	double ub = CPX_INFBOUND;
	double cost = 1.;
	status = CPXnewcols(env1, lp1, 1, &cost, &lb, &ub, NULL, NULL);
	if (status) {
		fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
		goto TERMINATE;
	}
	int num_rho_var = index1;
	index1 = 0;
	rmatbeg = 0;

	int *quadrow;
	int *quadcol;
	double *quadval;
	quadrow = create_int_vector(I);
	quadcol = create_int_vector(I);
	quadval = create_double_vector(I);
	rmatind = create_int_vector(I + 1);
	rmatval = create_double_vector(I + 1);
	rhs = 0;
	for (i = 0; i < I + 1; i++) {
		rhs -= xk[i] * xk[i];
	}
	sense = 'L';
	nzcnt = 0;
	int nzcnt1 = 0;
	rmatind[nzcnt1] = rho;
	rmatval[nzcnt1++] = -1.;
	for (i = 0; i < I; i++) {
		rmatind[nzcnt1] = x[i];
		rmatval[nzcnt1++] = -2.*xk[i];
	}
	for (i = 0; i < I; i++) {
		quadrow[nzcnt] = x[i];
		quadcol[nzcnt] = x[i];
		quadval[nzcnt++] = 1.;
	}
	
	status = CPXaddqconstr(env1, lp1, nzcnt1, nzcnt, rhs, sense, rmatind, rmatval, quadrow, quadcol, quadval, NULL);
	if (status) {
		fprintf(stderr, "Error in CPXaddrows: status = %d\n", status);
		goto TERMINATE;
	}
	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);
	free_and_null((char **)&quadrow);
	free_and_null((char **)&quadcol);
	free_and_null((char **)&quadval);

	rmatind = create_int_vector(I + S);
	rmatval = create_double_vector(I + S);

	printf("obj and ub is %.2lf\t%.2lf\t", obj, UB);
	rhs = (1 - lambda)*obj + lambda * UB;
	sense = 'L';
	nzcnt = 0;
	for (i = 0; i < I; i++) {
		rmatind[nzcnt] = x[i];
		rmatval[nzcnt++] = -1.*p[i];
	}
	for (i = 0; i < S; i++) {
		rmatind[nzcnt] = alpha[i];
		rmatval[nzcnt++] = 1.;
	}
	status = CPXaddrows(env1, lp1, 0, 1, nzcnt, &rhs, &sense,
		&rmatbeg, rmatind, rmatval, NULL, NULL);
	if (status) {
		fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
		goto TERMINATE;
	}
	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);

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

	status = CPXaddrows(env1, lp1, 0, 1, nzcnt, &rhs, &sense,
		&rmatbeg, rmatind, rmatval, NULL, NULL);
	if (status) {
		fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
		goto TERMINATE;
	}
	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);

	int k;
	rmatind = (int *)malloc((I + 1) * sizeof(int));
	if (rmatind == NULL) {
		fprintf(stderr, "No memory for rmatind array.\n");
		status = -1;
		goto TERMINATE;
	}
	rmatval = (double *)malloc((I + 1) * sizeof(double));
	if (rmatval == NULL) {
		fprintf(stderr, "No memory for rmatval array.\n");
		status = -1;
		goto TERMINATE;
	}
	sense = 'L';
	for (s = 0; s < S; s++) {
		for (k = 0; k < iter; k++) {
			rhs = 0.;
			for (i = 0; i < I; i++) {
				rhs += tau[k][s][tau13[i][i]];
				rhs -= tau[k][s][tau14[i][i]];
				rhs += tau[k][s][tau11[i]] * z[i][s];
			}

			nzcnt = 0;
			rmatind[nzcnt] = alpha[s];
			rmatval[nzcnt++] = -1.;
			for (i = 0; i < I; i++) {
				rmatind[nzcnt] = x[i];
				rmatval[nzcnt++] = tau[k][s][tau11[i]];
			}
			for (i = 0; i < I; i++) {
				if (ABS(rhs) > 0.000001 || ABS(tau[k][s][tau11[i]]) > 0.000001) {
					status = CPXaddrows(env1, lp1, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}
					break;
				}
			}


		}
	}
	free_and_null((char **)&rmatind);
	free_and_null((char **)&rmatval);

	CPXsetdblparam(env1, CPX_PARAM_TILIM, 7200);
	//CPXsetintparam(env1, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
	status = CPXbaropt(env1, lp1);
	if (status) {
		fprintf(stderr, "Error in CPXprimopt: status = %d\n", status);
		goto TERMINATE;
	}
	int worker_lp_sol_stat = CPXgetstat(env1, lp1);
	int num_sub = CPXgetnumcols(env1, lp1);
	printf("status is %d\t", worker_lp_sol_stat);
	if (worker_lp_sol_stat == CPX_STAT_OPTIMAL) {

		status = CPXgetx(env1, lp1, x1, 0, num_sub - 1);
		if (status) {
			fprintf(stderr, "Failed to obtain solution, status = %d.\n", status);
			goto TERMINATE;
		}

	}
	else if (worker_lp_sol_stat == CPX_STAT_NUM_BEST) {
		status = CPXgetx(env1, lp1, x1, 0, num_sub - 1);
		if (status) {
			fprintf(stderr, "Failed to obtain solution, status = %d.\n", status);
			goto TERMINATE;
		}
	}
	else {
		status = 1;
		goto TERMINATE;
	}
	if (lp1 != NULL) {
		int local_status = CPXfreeprob(env1, &lp1);
		if (local_status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n",
				local_status);
			status = local_status;
		}
	}

	/* Free the CPLEX environment, if necessary */

	if (env1 != NULL) {
		int local_status = CPXcloseCPLEX(&env1);
		if (local_status) {
			fprintf(stderr,
				"Could not close CPLEX environment, status = %d.\n",
				local_status);
			status = local_status;
		}
	}
TERMINATE:
	return status;
	//	}

}

static void free_and_null(char **ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */
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
