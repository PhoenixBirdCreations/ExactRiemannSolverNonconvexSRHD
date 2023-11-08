#include "integral_curves.h"

double v_rho(double l, double D, double vr, double r);
double dv_rho(double D, double l, double vr, double r);
double eq3(double D, double vr, double Pr, double rr, double l, double r);
double deq3(double D, double vr, double Pr, double rr, double l, double r);
double delicateBisection(double D, double vr, double P, double rr, double l, double xnew, int sense);
double newton_eq(double rr, double D, double vr, double P, double l, int sense, double prev);
double eps_rho(double vr, double P, double rr, double l, double r, double v);

void store_rarefaction(double **rare, double va, double Pa, double ra, double goal_r, int s, int N);

int build_mixedcurve(double *va, double *Pa, double *ra, double *vm, int s, double rb, double *st, double *ov, double *oP, double *or, double *ovw);
               
void mixed_interpolate(double va, double Pa, double ra, double rb, double exa, double exb, int zero_var,int N, double *xv, double *yv, int s);

void findRarefactionEnd_mixedsonic(double va, double Pa, double ra, double rb, double *vout, double *Pout, double *rout, double *vm, int s);

void last_mixed(double *vr, double *Pr, double *rr, double *vwr, double *v, double *P, double *r, double *vw, double vmiddle, double Pmiddle, int s);
