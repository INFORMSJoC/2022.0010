#ifndef _H_solveNVDROW
#define _H_solveNVDROW

typedef struct node
{
	int data;
	struct node *next;
}linklist;

double solveNVDROW(int I, int S, double *uz, double theta, double *p, double **z, double *c, double d, int Moos, double **xis_oos, double *prob, double *uu, double *prob_valid);//(int I, int S, double *ubar, double Gamma, double *bs, double **xis, double *cs, double d, int Moos, double **xis_oos) {
int ***create_int_tmatrix(int thr, int rows, int Columns);
int ****create_int_fmatrix(int fou, int thr, int rows, int Columns);
double **create_double_matrix(int rows, int Columns);
double ***create_double_tmatrix(int thr, int rows, int Columns);
void release_list(linklist *pHead);
#endif
