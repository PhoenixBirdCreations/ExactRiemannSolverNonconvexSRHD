#include "read_par_ic.h"

void m_fail(char*s);
void bigToSmall(double *v, int *out, int N);
void smallToBig(double *v, int *out, int N);
void orderArray(double *v, int N);
void purgeL(double *s, int *out, int N);
void purgeR(double *s, int *out, int N);

double parabola(double *x, double *y, int j, double eval, int N);
double recta(double *x, double *y, int j, double eval, int N);
double cap(double a);


