#include "hugoniot_curves.h"

double int_e_isen(double r);
double int_R_inv(double r, double eps);
int build_rarefaction(double *v, double *P, double *r, double *vm, int s, double goal_r, double *st, double *ov, double *op, double *or, double *ovw,int first);
void find_nonl0(double ra, double ea, double rb, double eb, double *r, double *P);

void integral_interpolate(double va, double Pa, double ra, double rb, double exa, double exb, int zero_var, int N, double *xv, double *yv, int s);
void find_rgoal(double v, double P, double r, int s, double Pmiddle, double *rgoal);
