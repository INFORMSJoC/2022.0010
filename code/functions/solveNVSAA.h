#ifndef _H_solveNVSAA
#define _H_solveNVSAA

double solveNVSAA(int I, int S, double *p, double **z, double *c, double d, int Moos, double **xis_oos);
void free_and_null(char **ptr);
int *create_int_vector(int dim);
int **create_int_matrix(int rows, int Columns);
double *create_double_vector(int dim);
#endif
