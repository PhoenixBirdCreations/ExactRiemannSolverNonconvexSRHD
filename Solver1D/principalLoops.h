#include "mixed_curves.h"
void incrementSign();
void firstWave(int *cL, int *cR);
void problemType(int *pr, double *ini_dif);
int possibilityIntersection(int oc, int side, double **waveExtremesL, double **waveExtremesR, int *wtl, int *wtr, int nwl, int nwr, int liL, int liR, int zero_var, double ini_dif,double *vmiddle, double *Pmiddle);
int checkIntersection(double **waveExtremesL, double **waveExtremesR, int *wTL, int *wTR, int iL, int iR, double *vmiddle, double *Pmiddle, int);
struct Data saveInfo(int nwl, int nwr, int *wtl, int *wtr, double **wel, double **wer, double vmiddle, double Pmiddle);
struct Data findWaveExtremes();
void exportWaveCurves(struct Data sol);
void exportExactSolution(struct Data sol);



