#include "auxFun.h"
double hb_Taub(double Pa, double Pb, double rb, double ha, double ra);
double F_hug_newton(double x, double Pb, double Pa, double ha, double ra);
double dF_hug_newton(double x, double Pb, double Pa, double ha, double ra);
double newton_hugoniot(double prev, double Pa, double ha, double ra, double Pb);

double jsquare(double Pa, double ra, double ha, double Pb, double rb, double hb);
double vshock(double ra, double va, double js, int s);
double vflow(double ha, double va, double Pa, double ra, double vs, double Pb, double j);

int build_shock(double *v, double *P, double *r, double *vw, int s, double goal_P, double Pb, double *st, double *ov, double *op, double *or, double *ovw, int first);
//void find_sonic(double ra, double Pa, double va, double na, double rb, double Pb, double nb, int s, double *r, double *P);

void hugoniot_interpolate(double va, double Pa, double ra, double Pgoal, double Pb, double vstart, double exa, double exb, int zero_var, int N, double *xv, double *yv,int s);
