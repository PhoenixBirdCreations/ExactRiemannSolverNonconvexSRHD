#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern double RL,VL,PL, RR,VR,PR;//conditions of the Riemann problem
extern double userGoalPl, userGoalPr, userGoalRl, userGoalRr; //initial goals for solution
extern double x_ini_disc, tf; //position initial discontinuity, final time of integration
extern char *wave_curves_file, *exact_solution_file;
extern FILE *fwc, *fes;

extern int sl, sr;
extern int savewc, savees;

struct Data{
    int n_wl; //number waves L (not counting the 4)
    int n_wr; //number waves R (not counting the 4)
    int *wtl; //type of waves L and the 4
    int *wtr; //type of waves R and the 4
    double **wcl; //extremes L for saving wave curves
    double **wcr; //extremes R for saving wave curves
    double **esl; //extremes L to draw exact solution
    double **esr; //extremes R to draw exact solution
}Data;

struct GGL{
    double g0;
    double g1;
    double s0;
    double r0;
};

extern struct GGL eos;

#define TOL 1.0e-12
#define MAXIT 500
#define INFTY 1e30
#define EPSILON 0.000001

#define SIGN(a) ((a) >= 0.0 ? 1 : -1)
#define MIN0(a) ((a)<0.0 ? 0 : a) 

double gmma(double r);
double dgmma(double r);
double ddgmma(double r);
double dddgmma(double r);

double nonlfactor(double css, double r);
double fundG(double r);
