#include "principalLoops.h"

double RL,VL,PL, RR,VR,PR; //conditions of the Riemann problem
double userGoalPl, userGoalPr, userGoalRl, userGoalRr; //initial goals for solution
double tf; //final time of integration
char *wave_curves_file, *exact_solution_file;

double x_ini_disc=0.5; //position initial discontinuity

FILE *fwc, *fes;

int sl, sr;
int savewc=0, savees=0;

struct GGL eos; //EoS


int main(){
    
    read_par_eos("eos.par");

    double *ciInfo;
    if((ciInfo=(double*)malloc(7*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    double *goalInfo;
    if((goalInfo=(double*)malloc(4*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    read_par_ic("ic.par", ciInfo, goalInfo, &exact_solution_file, &wave_curves_file);
    RL=ciInfo[0]; VL=ciInfo[1]; PL=ciInfo[2]; RR=ciInfo[3]; VR=ciInfo[4]; PR=ciInfo[5];
    tf=ciInfo[6];
    userGoalPl=goalInfo[0]; userGoalRl=goalInfo[1]; userGoalPr=goalInfo[2]; userGoalRr=goalInfo[3];

    int i;
    struct Data sol;
    fwc=fopen(wave_curves_file,"w");
    fes=fopen(exact_solution_file,"w");
    sol=findWaveExtremes();
    
    printf("Diagram v-P completed. Exporting results.\n\n");
    
    exportWaveCurves(sol);
    exportExactSolution(sol);
    fclose(fes);
    fclose(fwc);
    
    printf("Results exported. Cleaning and exiting\n\n");
 
    for(i=0; i<2*sol.n_wl; i++) {free(sol.wcl[i]);} free(sol.wcl); 
    for(i=0; i<2*sol.n_wr; i++) {free(sol.wcr[i]);} free(sol.wcr);
    for(i=0; i<2*sol.n_wl; i++) {free(sol.esl[i]);} free(sol.esl);
    for(i=0; i<2*sol.n_wr; i++) {free(sol.esr[i]);} free(sol.esr);
    free(sol.wtl); free(sol.wtr);
    
    printf("Done!\n");

    return 0;
}
