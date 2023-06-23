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
#include "CGsolver.h"
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define EPSILON	 0.000000000001    //Added to the fractional solution obtained at the node

void release_list(linklist *pHead);
void xiangjian(linklist *head1, linklist *head2);
void hebing(linklist *head1, linklist *head2);
void print(linklist *head);
double theta;


//int CGsovler(double *p, double *uz, double **z, linklist *Ly1, linklist *Ly11, linklist *Leta1, linklist *Leta2, linklist *Leta3, linklist *Leta4, linklist *Leta5, linklist *Leta6, CPXENVptr env, CPXLPptr lp, int *worker_lp_sol_stat, double *z1, double *x2, double *tau);
int main(int  argc, char *argv[]) {
	int		   num_inst, confirm;
	FILE       *ini;
	clock_t	   start_time1, end_time1;
	double     cputimeMIP1 = 0.;
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
	time_t		tt;
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
	int rk, I, S;
	double d;
	char	   instance[40];           // External file containing all the input data of the problem
	char	   instance1[40];           // External file containing all the input data of the problem
	int num_r = 0;
	double avg_obj = 0.;
	double avg_time = 0.;
	double avg_time1 = 0.;
	double avg_time2 = 0.;
	double avg_obj_iter = 0.;
	double avg_time_iter = 0.;
	double avg_time1_iter = 0.;
	double avg_time2_iter = 0.;
	double avg_iter = 0.;
	int num_fea = 0;
	double avg_col = 0.;
	int avg_col_iter = 0;
	double avg_sparse = 0.;
	double avg_sparse_iter = 0.;
	double avg_tcol = 0.;
	double avg_tcol_iter = 0.;
	for (rk = 0; rk < num_inst; rk++) {
		avg_obj_iter = 0.;
		avg_time_iter = 0.;
		avg_time1_iter = 0.;
		avg_time2_iter = 0.;
		cputimeMIP1 = 0.;
		avg_sparse_iter = 0.;
		avg_tcol_iter = 0.;
		confirm = 0;
		avg_col_iter = 0;
		/***********************Now reviewing the information for the instance**************************/
		if (fscanf(ini, "%d", &I) == 1) confirm++;
		if (fscanf(ini, "%d", &S) == 1) confirm++;
		if (fscanf(ini, "%s", &instance) == 1) confirm++;//data name
		if (fscanf(ini, "%lf", &theta) == 1) confirm++;
		if (confirm < 4)  goto TERMINATE;
		out = Open_File(outfile, "a+");
		fprintf(out, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", rk + 1, I, S, instance, inst_name, theta);
		fclose(out);

		solution = Open_File(xvalue, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", rk + 1, I, S, instance, inst_name, theta);
		fclose(solution);

		solution = Open_File(Solution, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", rk + 1, I, S, instance, inst_name, theta);
		fclose(solution);

		solution = Open_File(log1, "a+");
		fprintf(solution, "%d;\t %d\t %d\t %s\t%s\t%.1lf\n", rk + 1, I, S, instance, inst_name, theta);
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
		linklist *Ly1, *pl, *ql;
		linklist *Ly11, *pl11;
		linklist *Leta1, *ple1;
		linklist *Leta2, *ple2;
		linklist *Leta3, *ple3;
		linklist *Leta4, *ple4;
		linklist *Leta5, *ple5;
		linklist *Leta6, *ple6;

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
		int flag = 0;

		double *tau;

		int *x, alpha, *r, t, **y0, ***y1, **y2;
		int **eta1, **eta2, ***eta3, ***eta4, ***eta5, ***eta6;
		int *w1, **w2, **w3;

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


		int index1 = 0;
		for (i = 0; i < I; i++) {
			x[i] = index1;
			index1++;

		}
		int num_x_var = index1;
		index1 = 0;

		alpha = num_x_var + index1;
		index1++;

		int num_alpha_var = index1;
		index1 = 0;

		for (s = 0; s < S; s++) {
			r[s] = num_x_var + num_alpha_var + index1;
			index1++;

		}
		int num_r_var = index1;
		index1 = 0;

		t = num_r_var + num_x_var + num_alpha_var + index1;
		index1++;

		int num_t_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				y0[i][s] = num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;

			}
		}
		int num_y0_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					y1[i][j][s] = num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;
				}
			}
		}
		int num_y1_var = index1;
		index1 = 0;


		for (j = 0; j < I; j++) {
			for (s = 0; s < S; s++) {
				y2[j][s] = num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;

			}
		}
		int num_y2_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				eta1[i][s] = num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;

			}
		}
		int num_eta1_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				eta2[i][s] = num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;

			}
		}
		int num_eta2_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta3[i][j][s] = num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;

				}
			}
		}
		int num_eta3_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta4[i][j][s] = num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;

				}
			}
		}
		int num_eta4_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta5[i][j][s] = num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;

				}
			}
		}
		int num_eta5_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					eta6[i][j][s] = num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;

				}
			}
		}
		int num_eta6_var = index1;
		index1 = 0;

		for (s = 0; s < S; s++) {
			w1[s] = num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
			index1++;

		}

		int num_w1_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				w2[i][s] = num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;

			}
		}
		int num_w2_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				w3[i][s] = num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;

			}
		}
		int num_w3_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (j = 0; j < I; j++) {
				for (s = 0; s < S; s++) {
					y11[i][j][s] = num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
					index1++;

				}
			}
		}
		int num_y11_var = index1;
		index1 = 0;

		int alpha1, *r1, **y01;
		r1 = create_int_vector(S);
		y01 = create_int_matrix(I, S);

		alpha1 = num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
		index1++;

		int num_alpha1_var = index1;
		index1 = 0;

		for (s = 0; s < S; s++) {
			r1[s] = num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
			index1++;

		}
		int num_r1_var = index1;
		index1 = 0;

		for (i = 0; i < I; i++) {
			for (s = 0; s < S; s++) {
				y01[i][s] = num_r1_var + num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;
				index1++;

			}
		}
		int num_y01_var = index1;
		index1 = 0;

		int num_tau1 = num_y01_var + num_r1_var + num_alpha1_var + num_y11_var + num_w3_var + num_w2_var + num_w1_var + num_eta6_var + num_eta5_var + num_eta4_var + num_eta3_var + num_eta2_var + num_eta1_var + num_y2_var + num_y1_var + num_y0_var + num_t_var + num_r_var + num_x_var + num_alpha_var + index1;

		indexd = 0;

		tau21 = num_tau1 + indexd;
		indexd++;

		int num_tau2 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			tau31[s] = num_tau2 + num_tau1 + indexd;
			indexd++;

		}
		int num_tau3 = indexd;
		indexd = 0;


		for (s = 0; s < S; s++) {
			tau41[s] = num_tau3 + num_tau2 + num_tau1 + indexd;
			indexd++;

		}
		int num_tau4 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (j = 0; j < I; j++) {
				tau51[s][j] = num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;

			}
		}
		int num_tau5 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (j = 0; j < I; j++) {
				tau61[s][j] = num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;

			}
		}
		int num_tau6 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau71[s][i] = num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;

			}
		}
		int num_tau7 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau81[s][i] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;

			}
		}
		int num_tau8 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau91[s][i][j] = num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;

				}
			}
		}
		int num_tau9 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau101[s][i][j] = num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;

				}
			}
		}
		int num_tau10 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau111[s][i] = num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;

			}
		}
		int num_tau11 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				tau121[s][i] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
				indexd++;

			}
		}
		int num_tau12 = indexd;
		indexd = 0;

		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau131[s][i][j] = num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;

				}
			}
		}
		int num_tau13 = indexd;
		indexd = 0;


		for (s = 0; s < S; s++) {
			for (i = 0; i < I; i++) {
				for (j = 0; j < I; j++) {
					tau141[s][i][j] = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
					indexd++;

				}
			}
		}
		int num_tau14 = indexd;
		indexd = 0;

		int num_sub = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + num_tau14;
		int num_tau = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + 1 + num_tau14;
		CPXENVptr env = NULL;
		CPXLPptr  lp = NULL;
		int iter, Iter;
		iter = 0;
		Iter = 100;
		double fmin_y1 = 1000.;
		double fmin_y11 = 1000.;
		double *x2;
		x2 = create_double_vector(num_sub);
		tau = create_double_vector(num_tau);

		while (flag == 0 && iter < Iter &&avg_time1_iter<7200) {
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
			CPXsetdblparam(env, CPX_PARAM_TILIM, 7200);
			CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
			int worker_lp_sol_stat;
			double z1;
			start_time1 = clock();//CPXprimopt
			worker_lp_sol_stat = CGsolver(iter, d, c, p, uz, z, I, S, theta, Ly1, Ly11, Leta1, Leta2, Leta3, Leta4, Leta5, Leta6, env, lp, &z1, x2, tau);
			end_time1 = clock();
			cputimeMIP1 = (double)(end_time1 - start_time1) / (double)CLOCKS_PER_SEC;

			printf("status is %d\t", worker_lp_sol_stat);
			if (worker_lp_sol_stat == CPX_STAT_OPTIMAL) {



				min_y1 = 1000.;
				min_y11 = 1000.;
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

				double c_y1;
				int flag_y1 = 0;
				flag = 1;
				solution = Open_File(Solution, "a+");
				out = Open_File(log1, "a+");
				fmin_y1 = 1000.;
				fmin_y11 = 1000.;
				start_time1 = clock();

				min_eta2 = 1000.;
				min_eta1 = 1000.;

				double c_y1e, c_y2, c_y3, c_y4, c_y5, c_y6;
				for (i = 0; i < I; i++) {
					min_eta3 = 1000.;
					min_eta4 = 1000.;
					min_eta5 = 1000.;
					min_eta6 = 1000.;
					for (s = 0; s < S; s++) {
						c_y1e = -1.*(uz[i] - z[i][s])*tau[tau3[s]] + tau[tau5[s][i]] - tau[tau6[s][i]];
						c_y2 = -1.*(z[i][s])*tau[tau3[s]] - tau[tau5[s][i]] + tau[tau6[s][i]];
						for (j = 0; j < I; j++) {
							c_y3 = -1.*(uz[j] - z[j][s])*tau[tau7[s][i]] + tau[tau9[s][i][j]] - tau[tau10[s][i][j]];
							c_y4 = -1.*z[j][s] * tau[tau7[s][i]] - tau[tau9[s][i][j]] + tau[tau10[s][i][j]];
							c_y5 = -1.*(uz[j] - z[j][s])*tau[tau11[s][i]] + tau[tau13[s][i][j]] - tau[tau14[s][i][j]];
							c_y6 = -1.* z[j][s] * tau[tau11[s][i]] - tau[tau13[s][i][j]] + tau[tau14[s][i][j]];
							if (c_y3 < -0.00001) {
								ql = (linklist*)malloc(sizeof(linklist));
								ql->data = i * I*S + j * S + s;
								ple3->next = ql;                   //  利用尾差法  
								ple3 = ql;
								ple3->next = NULL;
								avg_col_iter++;
								flag = 0;
							}
							if (c_y4 < -0.00001) {
								ql = (linklist*)malloc(sizeof(linklist));
								ql->data = i * I*S + j * S + s;
								ple4->next = ql;                   //  利用尾差法  
								ple4 = ql;
								ple4->next = NULL;
								avg_col_iter++;
								flag = 0;

							}
							if (c_y5 < -0.00001) {
								ql = (linklist*)malloc(sizeof(linklist));
								ql->data = i * I*S + j * S + s;
								ple5->next = ql;                   //  利用尾差法  
								ple5 = ql;
								ple5->next = NULL;
								avg_col_iter++;
								flag = 0;
							}
							if (c_y6 < -0.00001) {
								ql = (linklist*)malloc(sizeof(linklist));
								ql->data = i * I*S + j * S + s;
								ple6->next = ql;                   //  利用尾差法  
								ple6 = ql;
								ple6->next = NULL;
								avg_col_iter++;
								flag = 0;
							}
						}
						if (c_y1e < -0.00001) {
							ql = (linklist*)malloc(sizeof(linklist));
							ql->data = i * S + s;
							ple1->next = ql;                   //  利用尾差法  
							ple1 = ql;
							ple1->next = NULL;
							avg_col_iter++;
							flag = 0;
						}
						if (c_y2 < -0.00001) {
							ql = (linklist*)malloc(sizeof(linklist));
							ql->data = i * S + s;
							ple2->next = ql;                   //  利用尾差法  
							ple2 = ql;
							ple2->next = NULL;
							avg_col_iter++;
							flag = 0;
						}
					}

				}

				int flag_y11 = 0;

				for (i = 0; i<I; i++) {
					min_y1 = 1000.;
					//	min_y11 = 1000.;
					for (j = 0; j<I; j++) {
						flag_y11 = 0;
						for (s = 0; s<S; s++) {
							c_y1 = -1.*tau[tau3[s]] * p[i] * z[j][s] - p[i] * tau[tau5[s][j]] + p[i] * tau[tau6[s][j]] + z[j][s] * tau[tau7[s][i]] + tau[tau9[s][i][j]] - tau[tau10[s][i][j]] + z[j][s] * tau[tau11[s][i]] + tau[tau13[s][i][j]] - tau[tau14[s][i][j]];
							if (c_y1<-0.00001) {
								ql = (linklist*)malloc(sizeof(linklist));
								ql->data = i * I*S + j * S + s;
								pl->next = ql;                   //  利用尾差法  
								pl = ql;
								pl->next = NULL;
								avg_col_iter++;
								flag = 0;
							}
						}

					}
				}
				for (i = 0; i<I; i++) {
					min_y11 = 1000.;
					for (j = 0; j<I; j++) {
						flag_y11 = 0;
						for (s = 0; s<S; s++) {
							//	c_y1 = tau[tau3[s]] * p[i] * z[j][s] + p[i] * tau[tau5[s][j]] - p[i] * tau[tau6[s][j]] - z[j][s] * tau[tau7[s][i]] - tau[tau9[s][i][j]] + tau[tau10[s][i][j]] - z[j][s] * tau[tau11[s][i]] - tau[tau13[s][i][j]] + tau[tau14[s][i][j]];
							//	c_y1 = -1.*tau[tau3[s]] * p[i] * z[j][s] - p[i] * tau[tau5[s][j]] + p[i] * tau[tau6[s][j]] + z[j][s] * tau[tau7[s][i]] + tau[tau9[s][i][j]] - tau[tau10[s][i][j]] + z[j][s] * tau[tau11[s][i]] + tau[tau13[s][i][j]] - tau[tau14[s][i][j]];
							double  c_y11 = tau[tau3[s]] * p[i] * z[j][s] + p[i] * tau[tau5[s][j]] - p[i] * tau[tau6[s][j]] - z[j][s] * tau[tau7[s][i]] - tau[tau9[s][i][j]] + tau[tau10[s][i][j]] - z[j][s] * tau[tau11[s][i]] - tau[tau13[s][i][j]] + tau[tau14[s][i][j]];

							if (c_y11<-0.00001) {
								ql = (linklist*)malloc(sizeof(linklist));
								ql->data = i * I*S + j * S + s;
								pl11->next = ql;                   //  利用尾差法  
								pl11 = ql;
								pl11->next = NULL;
								avg_col_iter++;
								flag = 0;
							}
						}

					}
				}
				end_time1 = clock();
				double cputimeMIPs = (double)(end_time1 - start_time1) / (double)CLOCKS_PER_SEC;
				cputimeMIPs = 0.;
				fprintf(solution, "%.6lf\t%.6lf\t%.6lf\n", cputimeMIPs, min_eta1, min_eta2);
				fclose(solution);

				int flag_y12 = 1;

				avg_time_iter += cputimeMIP1;

				if (flag == -1) {
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
					/* Init the CPLEX environment */

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
					//	tau1 = indexd;
					//	indexd++;
					double lb = -CPX_INFBOUND;
					double ub = CPX_INFBOUND;
					double cost = d;
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
						goto TERMINATE;
					}
					//	int num_tau1 = indexd;
					//	indexd = 0;

					//	tau2 = num_tau1 + indexd;
					//	indexd++;
					lb = -CPX_INFBOUND;
					ub = 0.;
					if (x2[tau21] > 0.000001) {
						lb = 0.;
					}
					status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
						goto TERMINATE;
					}
					//	int num_tau2 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						//		tau3[s] = num_tau2 + num_tau1 + indexd;
						//		indexd++;
						lb = -CPX_INFBOUND;
						ub = 0.;
						if (x2[tau31[s]] > 0.00001) {
							lb = 0.;
						}
						status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
						if (status) {
							fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
							goto TERMINATE;
						}
					}
					//	int num_tau3 = indexd;
					//	indexd = 0;


					for (s = 0; s < S; s++) {
						//		tau4[s] = num_tau3 + num_tau2 + num_tau1 + indexd;
						//		indexd++;
						lb = -CPX_INFBOUND;
						ub = 0.;
						if (x2[tau41[s]] > 0.00001) {
							lb = 0.;
						}
						status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
						if (status) {
							fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
							goto TERMINATE;
						}
					}
					//	int num_tau4 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (j = 0; j < I; j++) {
							//			tau5[s][j] = num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							//			indexd++;
							lb = -CPX_INFBOUND;
							ub = 0.;
							if (x2[tau51[s][j]] > 0.00001) {
								lb = 0.;
							}
							status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}
					//	int num_tau5 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (j = 0; j < I; j++) {
							//			tau6[s][j] = num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							//			indexd++;
							lb = -CPX_INFBOUND;
							ub = 0.;
							if (x2[tau61[s][j]] > 0.00001) {
								lb = 0.;
							}
							status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}
					//	int num_tau6 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							//			tau7[s][i] = num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							//			indexd++;
							lb = -CPX_INFBOUND;
							ub = 0.;
							if (x2[tau71[s][i]] > 0.00001) {
								lb = 0.;
							}
							status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}
					//	int num_tau7 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							//			tau8[s][i] = num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							//			indexd++;
							lb = -CPX_INFBOUND;
							ub = 0.;
							if (x2[tau81[s][i]] > 0.00001) {
								lb = 0.;
							}
							status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}
					//	int num_tau8 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								//				tau9[s][i][j] = num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
								//				indexd++;
								lb = -CPX_INFBOUND;
								ub = 0.;
								if (x2[tau91[s][i][j]] > 0.00001) {
									lb = 0.;
								}
								status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
								if (status) {
									fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
									goto TERMINATE;
								}
							}
						}
					}
					//	int num_tau9 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								//				tau10[s][i][j] = num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
								//				indexd++;
								lb = -CPX_INFBOUND;
								ub = 0.;
								if (x2[tau101[s][i][j]] > 0.00001) {
									lb = 0.;
								}
								status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
								if (status) {
									fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
									goto TERMINATE;
								}
							}
						}
					}
					//	int num_tau10 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							//			tau11[s][i] = num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							//			indexd++;
							lb = -CPX_INFBOUND;
							ub = 0.;
							cost = z[i][s];
							if (x2[tau111[s][i]] > 0.00001) {
								lb = 0.;
							}
							status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}
					//	int num_tau11 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							//			tau12[s][i] = num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
							//			indexd++;
							lb = -CPX_INFBOUND;
							ub = 0.;
							if (x2[tau121[s][i]] > 0.00001) {
								lb = 0.;
							}
							status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}
					//	int num_tau12 = indexd;
					//	indexd = 0;

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								//				tau13[s][i][j] = num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
								//				indexd++;
								lb = -CPX_INFBOUND;
								ub = 0.;
								cost = 0.;
								if (i == j) {
									cost = 1.;
								}
								if (x2[tau131[s][i][j]] > 0.00001) {
									lb = 0.;
								}
								status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
								if (status) {
									fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
									goto TERMINATE;
								}
							}
						}
					}
					//	int num_tau13 = indexd;
					//	indexd = 0;


					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								//				tau14[s][i][j] = num_tau13 + num_tau12 + num_tau11 + num_tau10 + num_tau9 + num_tau8 + num_tau7 + num_tau6 + num_tau5 + num_tau4 + num_tau3 + num_tau2 + num_tau1 + indexd;
								//				indexd++;
								lb = -CPX_INFBOUND;
								ub = 0.;
								cost = 0.;
								if (i == j) {
									cost = -1.;
								}
								if (x2[tau141[s][i][j]] > 0.00001) {
									lb = 0.;
								}
								status = CPXnewcols(env, lp, 1, NULL, &lb, &ub, NULL, NULL);
								if (status) {
									fprintf(stderr, "Error in CPXnewcols, status = %d.\n", status);
									goto TERMINATE;
								}
							}
						}
					}
					//	int num_tau14 = indexd;
					//	indexd = 0;

					CPXchgobjsen(env, lp, CPX_MAX);
					char sense;
					int nzcnt, rmatbeg, *rmatind = NULL;
					double rhs, *rmatval = NULL;
					rmatbeg = 0;

					rmatind = (int *)malloc((S + 1) * sizeof(int));
					if (rmatind == NULL) {
						fprintf(stderr, "No memory for rmatind array.\n");
						status = -1;
						goto TERMINATE;
					}
					rmatval = (double *)malloc((S + 1) * sizeof(double));
					if (rmatval == NULL) {
						fprintf(stderr, "No memory for rmatval array.\n");
						status = -1;
						goto TERMINATE;
					}

					//x
					for (i = 0; i < I; i++) {
						if (x2[x[i]] == 0.) {
							sense = 'L';
						}
						if (x2[x[i]] > 0.00001) {
							sense = 'E';
						}
						nzcnt = 0;
						rmatind[nzcnt] = tau1;
						rmatval[nzcnt++] = c[i];
						for (s = 0; s < S; s++) {
							rmatind[nzcnt] = tau11[s][i];
							rmatval[nzcnt++] = 1.;
						}
						rhs = -p[i];

						status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
							&rmatbeg, rmatind, rmatval, NULL, NULL);
						if (status) {
							fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
							goto TERMINATE;
						}
					}


					free_and_null((char **)&rmatind);
					free_and_null((char **)&rmatval);

					rmatind = (int *)malloc(1 * sizeof(int));
					if (rmatind == NULL) {
						fprintf(stderr, "No memory for rmatind array.\n");
						status = -1;
						goto TERMINATE;
					}
					rmatval = (double *)malloc(1 * sizeof(double));
					if (rmatval == NULL) {
						fprintf(stderr, "No memory for rmatval array.\n");
						status = -1;
						goto TERMINATE;
					}
					//alpha
					if (x2[alpha] == 0.) {
						sense = 'L';
					}
					if (x2[alpha] > 0.00001) {
						sense = 'E';
					}

					nzcnt = 0;
					rmatind[nzcnt] = tau2;
					rmatval[nzcnt++] = -1.;
					rhs = 1.;

					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}
					//alpha1
					if (x2[alpha1] == 0.) {
						sense = 'L';
					}
					if (x2[alpha1] > 0.00001) {
						sense = 'E';
					}
					nzcnt = 0;
					rmatind[nzcnt] = tau2;
					rmatval[nzcnt++] = 1.;
					rhs = 1.;

					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}

					free_and_null((char **)&rmatind);
					free_and_null((char **)&rmatval);

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
					//r

					for (s = 0; s < S; s++) {
						if (x2[r[s]] == 0.) {
							sense = 'L';
						}
						if (x2[r[s]] > 0.00001) {
							sense = 'E';
						}
						nzcnt = 0;
						rmatind[nzcnt] = tau2;
						rmatval[nzcnt++] = 1. / S;
						rmatind[nzcnt] = tau3[s];
						rmatval[nzcnt++] = -1.;
						rhs = 0.;

						status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
							&rmatbeg, rmatind, rmatval, NULL, NULL);
						if (status) {
							fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
							goto TERMINATE;
						}
					}

					//r1

					for (s = 0; s < S; s++) {
						if (x2[r1[s]] == 0.) {
							sense = 'L';
						}
						if (x2[r1[s]] > 0.00001) {
							sense = 'E';
						}
						nzcnt = 0;
						rmatind[nzcnt] = tau2;
						rmatval[nzcnt++] = -1. / S;
						rmatind[nzcnt] = tau3[s];
						rmatval[nzcnt++] = 1.;
						rhs = 0.;

						status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
							&rmatbeg, rmatind, rmatval, NULL, NULL);
						if (status) {
							fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
							goto TERMINATE;
						}
					}
					free_and_null((char **)&rmatind);
					free_and_null((char **)&rmatval);

					rmatind = (int *)malloc((S + 1) * sizeof(int));
					if (rmatind == NULL) {
						fprintf(stderr, "No memory for rmatind array.\n");
						status = -1;
						goto TERMINATE;
					}
					rmatval = (double *)malloc((S + 1) * sizeof(double));
					if (rmatval == NULL) {
						fprintf(stderr, "No memory for rmatval array.\n");
						status = -1;
						goto TERMINATE;
					}
					//t
					if (x2[t] > 0.00001) {
						sense = 'E';
					}
					if (x2[t] == 0.) {
						sense = 'L';
					}
					nzcnt = 0;
					rmatind[nzcnt] = tau2;
					rmatval[nzcnt++] = theta;
					for (s = 0; s < S; s++) {
						rmatind[nzcnt] = tau4[s];
						rmatval[nzcnt++] = -1.;
					}
					rhs = 0.;

					status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
						&rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status) {
						fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
						goto TERMINATE;
					}

					free_and_null((char **)&rmatind);
					free_and_null((char **)&rmatval);

					rmatind = (int *)malloc((2 * I + 1) * sizeof(int));
					if (rmatind == NULL) {
						fprintf(stderr, "No memory for rmatind array.\n");
						status = -1;
						goto TERMINATE;
					}
					rmatval = (double *)malloc((2 * I + 1) * sizeof(double));
					if (rmatval == NULL) {
						fprintf(stderr, "No memory for rmatval array.\n");
						status = -1;
						goto TERMINATE;
					}
					//w3

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							if (x2[w3[i][s]] > 0.00001) {
								sense = 'E';
							}
							if (x2[w3[i][s]] == 0.) {
								sense = 'L';
							}
							nzcnt = 0;
							rmatind[nzcnt] = tau12[s][i];
							rmatval[nzcnt++] = 1.;
							for (j = 0; j < I; j++) {
								rmatind[nzcnt] = tau13[s][i][j];
								rmatval[nzcnt++] = -1.;
							}
							for (j = 0; j < I; j++) {
								rmatind[nzcnt] = tau14[s][i][j];
								rmatval[nzcnt++] = -1.;
							}
							rhs = 0.;

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
					//y0

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							if (x2[y0[i][s]] > 0.00001) {
								sense = 'E';
							}
							if (x2[y0[i][s]] == 0.) {
								sense = 'L';
							}
							nzcnt = 0;
							rmatind[nzcnt] = tau3[s];
							rmatval[nzcnt++] = p[i];
							rmatind[nzcnt] = tau7[s][i];
							rmatval[nzcnt++] = -1.;
							rmatind[nzcnt] = tau11[s][i];
							rmatval[nzcnt++] = -1.;
							rhs = 0.;

							status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
								&rmatbeg, rmatind, rmatval, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
								goto TERMINATE;
							}
						}

					}

					//y01

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							if (x2[y01[i][s]] > 0.00001) {
								sense = 'E';
							}
							if (x2[y01[i][s]] == 0.) {
								sense = 'L';
							}
							nzcnt = 0;
							rmatind[nzcnt] = tau3[s];
							rmatval[nzcnt++] = -p[i];
							rmatind[nzcnt] = tau7[s][i];
							rmatval[nzcnt++] = 1.;
							rmatind[nzcnt] = tau11[s][i];
							rmatval[nzcnt++] = 1.;
							rhs = 0.;

							status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
								&rmatbeg, rmatind, rmatval, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
								goto TERMINATE;
							}
						}

					}

					//y2

					//		tau[tau8[s][i]] + tau[tau12[s][i]] - tau[tau4[s]] * p[i];

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							if (x2[y2[i][s]] > 0.00001) {
								sense = 'E';
							}
							if (x2[y2[i][s]] == 0.) {
								sense = 'L';
							}

							nzcnt = 0;
							rmatind[nzcnt] = tau4[s];
							rmatval[nzcnt++] = p[i];
							rmatind[nzcnt] = tau8[s][i];
							rmatval[nzcnt++] = -1.;
							rmatind[nzcnt] = tau12[s][i];
							rmatval[nzcnt++] = -1.;
							rhs = 0.;

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
					//c_y1 = -1.*tau[tau3[s]] * p[i] * z[j][s] - p[i] * tau[tau5[s][j]] + p[i] * tau[tau6[s][j]] + z[j][s] * tau[tau7[s][i]] + tau[tau9[s][i][j]] - tau[tau10[s][i][j]] + z[j][s] * tau[tau11[s][i]] + tau[tau13[s][i][j]] - tau[tau14[s][i][j]];
					rmatind = (int *)malloc(9 * sizeof(int));
					if (rmatind == NULL) {
						fprintf(stderr, "No memory for rmatind array.\n");
						status = -1;
						goto TERMINATE;
					}
					rmatval = (double *)malloc(9 * sizeof(double));
					if (rmatval == NULL) {
						fprintf(stderr, "No memory for rmatval array.\n");
						status = -1;
						goto TERMINATE;
					}
					//y1
					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								if (x2[y1[i][j][s]] > 0.00001) {
									sense = 'E';
								}
								if (x2[y1[i][j][s]] == 0.) {
									sense = 'L';
								}
								rhs = 0.;
								nzcnt = 0;
								rmatind[nzcnt] = tau3[s];
								rmatval[nzcnt++] = p[i] * z[j][s];
								rmatind[nzcnt] = tau5[s][j];
								rmatval[nzcnt++] = p[i];
								rmatind[nzcnt] = tau6[s][j];
								rmatval[nzcnt++] = -p[i];
								rmatind[nzcnt] = tau7[s][i];
								rmatval[nzcnt++] = -z[j][s];
								rmatind[nzcnt] = tau9[s][i][j];
								rmatval[nzcnt++] = -1.;
								rmatind[nzcnt] = tau10[s][i][j];
								rmatval[nzcnt++] = 1.;
								rmatind[nzcnt] = tau11[s][i];
								rmatval[nzcnt++] = -z[j][s];
								rmatind[nzcnt] = tau13[s][i][j];
								rmatval[nzcnt++] = -1.;
								rmatind[nzcnt] = tau14[s][i][j];
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

					//y11
					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								if (x2[y11[i][j][s]] > 0.00001) {
									sense = 'E';
								}
								if (x2[y11[i][j][s]] == 0.) {
									sense = 'L';
								}
								rhs = 0.;
								nzcnt = 0;
								rmatind[nzcnt] = tau3[s];
								rmatval[nzcnt++] = -p[i] * z[j][s];
								rmatind[nzcnt] = tau5[s][j];
								rmatval[nzcnt++] = -p[i];
								rmatind[nzcnt] = tau6[s][j];
								rmatval[nzcnt++] = p[i];
								rmatind[nzcnt] = tau7[s][i];
								rmatval[nzcnt++] = z[j][s];
								rmatind[nzcnt] = tau9[s][i][j];
								rmatval[nzcnt++] = 1.;
								rmatind[nzcnt] = tau10[s][i][j];
								rmatval[nzcnt++] = -1.;
								rmatind[nzcnt] = tau11[s][i];
								rmatval[nzcnt++] = z[j][s];
								rmatind[nzcnt] = tau13[s][i][j];
								rmatval[nzcnt++] = 1.;
								rmatind[nzcnt] = tau14[s][i][j];
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

					//eta1
					//	c_y1e = -1.*(uz[i] - z[i][s])*tau[tau3[s]] + tau[tau5[s][i]] - tau[tau6[s][i]];
					//	c_y2 = -1.*(z[i][s])*tau[tau3[s]] - tau[tau5[s][i]] + tau[tau6[s][i]];
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
					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							if (x2[eta1[i][s]] > 0.00001) {
								sense = 'E';
							}
							if (x2[eta1[i][s]] == 0.) {
								sense = 'L';
							}
							rhs = 0.;
							nzcnt = 0;
							rmatind[nzcnt] = tau3[s];
							rmatval[nzcnt++] = uz[i] - z[i][s];
							rmatind[nzcnt] = tau5[s][i];
							rmatval[nzcnt++] = -1.;
							rmatind[nzcnt] = tau6[s][i];
							rmatval[nzcnt++] = 1.;
							status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
								&rmatbeg, rmatind, rmatval, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}

					//eta2

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							if (x2[eta2[i][s]] > 0.00001) {
								sense = 'E';
							}
							if (x2[eta2[i][s]] == 0.) {
								sense = 'L';
							}
							rhs = 0.;
							nzcnt = 0;
							rmatind[nzcnt] = tau3[s];
							rmatval[nzcnt++] = z[i][s];
							rmatind[nzcnt] = tau5[s][i];
							rmatval[nzcnt++] = 1.;
							rmatind[nzcnt] = tau6[s][i];
							rmatval[nzcnt++] = -1.;
							status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
								&rmatbeg, rmatind, rmatval, NULL, NULL);
							if (status) {
								fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
								goto TERMINATE;
							}
						}
					}

					//	c_y3 = -1.*(uz[j] - z[j][s])*tau[tau7[s][i]] + tau[tau9[s][i][j]] - tau[tau10[s][i][j]];
					//eta3

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								if (x2[eta3[i][j][s]] > 0.00001) {
									sense = 'E';
								}
								if (x2[eta3[i][j][s]] == 0.) {
									sense = 'L';
								}
								rhs = 0.;
								nzcnt = 0;
								rmatind[nzcnt] = tau7[s][i];
								rmatval[nzcnt++] = uz[j] - z[j][s];
								rmatind[nzcnt] = tau9[s][i][j];
								rmatval[nzcnt++] = -1.;
								rmatind[nzcnt] = tau10[s][i][j];
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
					//	c_y4 = -1.*z[j][s] * tau[tau7[s][i]] - tau[tau9[s][i][j]] + tau[tau10[s][i][j]];
					//eta4

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								if (x2[eta4[i][j][s]] > 0.00001) {
									sense = 'E';
								}
								if (x2[eta4[i][j][s]] == 0.) {
									sense = 'L';
								}
								rhs = 0.;
								nzcnt = 0;
								rmatind[nzcnt] = tau7[s][i];
								rmatval[nzcnt++] = z[j][s];
								rmatind[nzcnt] = tau9[s][i][j];
								rmatval[nzcnt++] = 1.;
								rmatind[nzcnt] = tau10[s][i][j];
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
					//	c_y5 = -1.*(uz[j] - z[j][s])*tau[tau11[s][i]] + tau[tau13[s][i][j]] - tau[tau14[s][i][j]];
					//eta5

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								if (x2[eta5[i][j][s]] > 0.00001) {
									sense = 'E';
								}
								if (x2[eta5[i][j][s]] == 0.) {
									sense = 'L';
								}
								rhs = 0.;
								nzcnt = 0;
								rmatind[nzcnt] = tau11[s][i];
								rmatval[nzcnt++] = uz[j] - z[j][s];
								rmatind[nzcnt] = tau13[s][i][j];
								rmatval[nzcnt++] = -1.;
								rmatind[nzcnt] = tau14[s][i][j];
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

					//c_y6 = -1.* z[j][s] * tau[tau11[s][i]] - tau[tau13[s][i][j]] + tau[tau14[s][i][j]];
					//eta6

					for (s = 0; s < S; s++) {
						for (i = 0; i < I; i++) {
							for (j = 0; j < I; j++) {
								if (x2[eta6[i][j][s]] > 0.00001) {
									sense = 'E';
								}
								if (x2[eta6[i][j][s]] == 0.) {
									sense = 'L';
								}
								rhs = 0.;
								nzcnt = 0;
								rmatind[nzcnt] = tau11[s][i];
								rmatval[nzcnt++] = z[j][s];
								rmatind[nzcnt] = tau13[s][i][j];
								rmatval[nzcnt++] = 1.;
								rmatind[nzcnt] = tau14[s][i][j];
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

					//w1
					rmatind = (int *)malloc((2 * I + 1) * sizeof(int));
					if (rmatind == NULL) {
						fprintf(stderr, "No memory for rmatind array.\n");
						status = -1;
						goto TERMINATE;
					}
					rmatval = (double *)malloc((2 * I + 1) * sizeof(double));
					if (rmatval == NULL) {
						fprintf(stderr, "No memory for rmatval array.\n");
						status = -1;
						goto TERMINATE;
					}
					for (s = 0; s < S; s++) {
						if (x2[w1[s]] > 0.00001) {
							sense = 'E';
						}
						if (x2[w1[s]] == 0.) {
							sense = 'L';
						}
						rhs = 0.;
						nzcnt = 0;
						rmatind[nzcnt] = tau4[s];
						rmatval[nzcnt++] = 1.;
						for (j = 0; j < I; j++) {
							rmatind[nzcnt] = tau5[s][j];
							rmatval[nzcnt++] = -1.;
						}
						for (j = 0; j < I; j++) {
							rmatind[nzcnt] = tau6[s][j];
							rmatval[nzcnt++] = -1.;
						}
						status = CPXaddrows(env, lp, 0, 1, nzcnt, &rhs, &sense,
							&rmatbeg, rmatind, rmatval, NULL, NULL);
						if (status) {
							fprintf(stderr, "Error in CPXaddrows, status = %d.\n", status);
							goto TERMINATE;
						}
					}

					//w2
					for (i = 0; i < I; i++) {
						for (s = 0; s < S; s++) {
							if (x2[w2[i][s]] > 0.00001) {
								sense = 'E';
							}
							if (x2[w2[i][s]] == 0.) {
								sense = 'L';
							}
							rhs = 0.;
							nzcnt = 0;
							rmatind[nzcnt] = tau8[s][i];
							rmatval[nzcnt++] = 1.;
							//- (j = 0; j<I; j++) {
							for (j = 0; j < I; j++) {
								rmatind[nzcnt] = tau9[s][i][j];
								rmatval[nzcnt++] = -1.;
							}
							for (j = 0; j < I; j++) {
								rmatind[nzcnt] = tau10[s][i][j];
								rmatval[nzcnt++] = -1.;
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
					CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
					start_time1 = clock();
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
						cputimeMIP2 = (double)(end_time1 - start_time1) / (double)CLOCKS_PER_SEC;

						double z1;
						status = CPXgetobjval(env, lp, &z1);
						if (status) {
							fprintf(stderr, "Failed to obtain solution, status = %d.\n", status);
							goto TERMINATE;
						}
						//get optimal solution

						solution = Open_File(outfile, "a+");
						fprintf(solution, "%.6lf\t%.6lf\n", z1, cputimeMIP2);
						fclose(solution);
						flag = 1;
					}
					else {
						end_time1 = clock();
						cputimeMIP2 = (double)(end_time1 - start_time1) / (double)CLOCKS_PER_SEC;
						solution = Open_File(outfile, "a+");
						fprintf(solution, "%.6lf\n", cputimeMIP2);
						fclose(solution);
					}
				}
				//	print(Ly1);
				if (flag == 1) {
					for (i = 0; i < num_sub; i++) {
						if (x2[i] > 0.00001) {
							avg_sparse_iter = avg_sparse_iter + 1.;
						}
					}
					avg_tcol_iter = num_sub - (I*I*S * 6 + I * S * 2) + avg_col_iter;
				}
				avg_obj_iter = z1;
				avg_time1_iter += cputimeMIP1 + cputimeMIPs + cputimeMIP2;
				avg_time2_iter += cputimeMIP2;
				solution = Open_File(outfile, "a+");
				fprintf(solution, "%.4lf\t %.4lf\t %.4lf\t %.4lf\t %d\n", z1, avg_time_iter, avg_time2_iter, avg_time1_iter, avg_col_iter);
				fclose(solution);


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

		}


		if (num_r < 10) {
			if (avg_time1_iter < 7200) {
				avg_time1 += avg_time1_iter;
				avg_time += avg_time_iter;
				avg_time2 += avg_time2_iter;
				num_fea++;
				avg_obj += avg_obj_iter;
				avg_sparse += avg_sparse_iter;
				avg_tcol += avg_tcol_iter;
			}
			avg_iter += iter;
			avg_col += avg_col_iter;
			if (num_r == 9) {
				solution = Open_File(outfile, "a+");
				fprintf(solution, "\n");
				fprintf(solution, "%.4lf\t %.4lf\t %.4lf\t %.4lf\t %.4lf\t %.1lf\t %d\t %.0lf\t %.0lf\n", avg_obj / num_fea, avg_time / num_fea, avg_time1 / num_fea, avg_time2 / num_fea, avg_iter / 10., avg_col / 10., num_fea, avg_sparse*1. / num_fea, avg_tcol*1. / num_fea);
				fclose(solution);

			}
			num_r++;

		}
		if (num_r == 10) {
			num_r = 0;
			avg_obj = 0.;
			avg_time = 0.;
			avg_time1 = 0.;
			avg_time2 = 0.;
			avg_iter = 0.;
			num_fea = 0;
			avg_col = 0.;
			avg_sparse = 0.;
			avg_tcol = 0.;
		}
		/* Init the CPLEX environment */
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
		free_and_null((char **)&p);
		free_and_null((char **)&c);
		free_and_null((char **)&uz);
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
			free_and_null((char **)&y21[i]);
			free_and_null((char **)&eta1[i]);
			free_and_null((char **)&eta2[i]);
			free_and_null((char **)&w2[i]);
			free_and_null((char **)&w3[i]);
			free_and_null((char **)&z[i]);
		}
		free_and_null((char **)&y0);
		free_and_null((char **)&y2);
		free_and_null((char **)&y21);
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
void xiangjian(linklist *head1, linklist *head2)         //集合的差运算  
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
void hebing(linklist *head1, linklist *head2)     // 集合的合并  
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
