#include "principalLoops.h"

void incrementSign(){
    /*
     decide direction of variables. Use H curves bc easy, H and I are tangent at initial conditions so the result has to be equivalent although the other wave is needed
     */
    
    double g, va, Pa, ra, ha, js, vs, Pb, rb, hb, vb;
    double norm, dist0;
    
    if(SIGN(VL)!=SIGN(VR)){
        if(VL<0){ //expanding slabs
            sl=-1; sr=-1;
        }
        else{ //colliding slabs
            sl=1; sr=1;
        }
    }
    else{
        dist0=((PL-PR)*(PL-PR)+(VL-VR)*(VL-VR));
        
        //first the left side
        Pa=PL; ra=RL; va=VL;
        
        g=gmma(ra);
        ha=1.0+Pa/ra*(g/(g-1));
        
        Pb=Pa+0.1; rb=ra;   
        rb=newton_hugoniot(rb, Pa, ha, ra, Pb);
        if (rb<0) sl=-1; //there was a problem with Newton in that direction, so the valid one is the other: P decreases
        else{
            g=gmma(rb);
            hb=1.0+Pb*g/(g-1.0)/rb; if(hb<0){printf("\nERROR:Not valid hb=%f\n",hb); exit(0);}
            js=jsquare(Pa, ra, ha, Pb, rb, hb);
            vs=vshock(ra, va, js, -1);
            vb=vflow(ha,va,Pa,ra,vs,Pb,-sqrt(js)); //going to left
            norm=((Pb-PR)*(Pb-PR)+(vb-VR)*(vb-VR));
            if(norm<dist0) sl=1; else sl=-1;
        }
        
        //then the right side
        Pa=PR; ra=RR; va=VR;
        
        g=gmma(ra);
        ha=1.0+Pa/ra*(g/(g-1));
        
        Pb=Pa+0.1; rb=ra;
        rb=newton_hugoniot(rb, Pa, ha, ra, Pb);
        if (rb<0) sr=-1; //there was a problem with Newton in that direction, so the valid one is the other: P decreases
        else{
            g=gmma(rb);
            hb=1.0+Pb*g/(g-1.0)/rb; if(hb<0){printf("\nERROR:Not valid hb=%f\n",hb); exit(0);}
            js=jsquare(Pa, ra, ha, Pb, rb, hb);
            vs=vshock(ra, va, js, 1);
            vb=vflow(ha,va,Pa,ra,vs,Pb,sqrt(js)); //going to right
            norm=((Pb-PL)*(Pb-PL)+(vb-VL)*(vb-VL));
            if(norm<dist0) sr=1; else sr=-1;
        }
    }
}

void firstWave(int *cL, int *cR){
    /*
     decide type of first wave from each side
    */
    
    double h, css, g;
    
    g=gmma(RL);
    h=1.0+PL/RL*(g/(g-1));
    css=PL/(g-1)/RL/h*(g*(g-1)+RL*dgmma(RL));
    if(nonlfactor(css,RL)<0){
        if(sl<0) *cL=3;
        else *cL=1;
    }
    else{
        if(sl<0) *cL=1;
        else *cL=3;
    }
    
    g=gmma(RR);
    h=1.0+PR/RR*(g/(g-1));
    css=PR/(g-1)/RR/h*(g*(g-1)+RR*dgmma(RR));
    if(nonlfactor(css,RR)<0){
        if(sr<0) *cR=3;
        else *cR=1;
    }
    else{
        if(sr<0) *cR=1;
        else *cR=3;
    }

    printf("LEFT: first wave curve is "); if((*cL)==1) printf("an integral curve.\n"); if((*cL)==3) printf("a Hugoniot curve.\n");
    printf("      Pressure is "); if(sl<0) printf("decreasing.\n"); else printf("increasing.\n");
    printf("\n");
    printf("RIGHT: first wave curve is "); if((*cR)==1) printf("an integral curve.\n"); if((*cR)==3) printf("a Hugoniot curve.\n");
    printf("       Pressure is "); if(sr<0) printf("decreasing.\n"); else printf("increasing.\n");
}

void problemType(int *pr, double *ini_dif){
    if(PL==PR) {*pr=0; *ini_dif=VL-VR;}
    else {*pr=1; *ini_dif=PL-PR;}
}

int possibilityIntersection(int oc, int side, double **waveExtremesL, double **waveExtremesR, int *wtl, int *wtr, int nwl, int nwr, int liL, int liR, int zero_var, double ini_dif,double *vmiddle, double *Pmiddle){
    int iniL, finL, iniR, finR;
    if(side<0) {iniL=oc; finL=oc; iniR=0; finR=2*nwr-1-2*liR;}
    else {iniR=oc; finR=oc; iniL=0; finL=2*nwl-1-2*liL;}
    
    int i,j, iL, iR;
    for (i=iniL; i<=finL; i=i+2){
        for(j=iniR; j<=finR; j=j+2){
            if(((waveExtremesL[i+1][zero_var]-waveExtremesR[j+1][zero_var])*ini_dif<=0) && //change of sign in difference
                (//curves coincide in some domain
                    ((waveExtremesL[MIN0(i-1)][((zero_var+1)%2)]<=waveExtremesR[MIN0(j-1)][((zero_var+1)%2)]) && (waveExtremesR[MIN0(j-1)][((zero_var+1)%2)]<waveExtremesL[i+1][((zero_var+1)%2)] )) ||
                    ((waveExtremesR[MIN0(j-1)][((zero_var+1)%2)]<=waveExtremesL[MIN0(i-1)][((zero_var+1)%2)]) && (waveExtremesL[MIN0(i-1)][((zero_var+1)%2)]<waveExtremesR[j+1][((zero_var+1)%2)] )) ||
                    ((waveExtremesL[MIN0(i-1)][((zero_var+1)%2)]>=waveExtremesR[MIN0(j-1)][((zero_var+1)%2)]) && (waveExtremesR[MIN0(j-1)][((zero_var+1)%2)]>waveExtremesL[i+1][((zero_var+1)%2)] )) ||
                    ((waveExtremesR[MIN0(j-1)][((zero_var+1)%2)]>=waveExtremesL[MIN0(i-1)][((zero_var+1)%2)]) && (waveExtremesL[MIN0(i-1)][((zero_var+1)%2)]>waveExtremesR[j+1][((zero_var+1)%2)] )) ||
                    ((waveExtremesR[MIN0(j-1)][((zero_var+1)%2)]==waveExtremesL[MIN0(i-1)][((zero_var+1)%2)] && sl==sr))
                )
            ){
                iL=i; iR=j;
                if(checkIntersection(waveExtremesL, waveExtremesR, wtl, wtr, iL, iR, vmiddle, Pmiddle, zero_var)){
                    if (((int)(iL*0.5+1))!=nwl) wtl[(int)(iL*0.5+1)]=4;
                    if (((int)(iR*0.5+1))!=nwr) wtr[(int)(iR*0.5+1)]=4;
                    return 1;
                }
            }
        }
    }
    return 0; 
}

int checkIntersection(double **waveExtremesL, double **waveExtremesR, int *wTL, int *wTR, int iL, int iR, double *vmiddle, double *Pmiddle, int zero_var){
    double *x_vectorL, *y_vectorL, *x_vectorR, *y_vectorR, *aux; // *dd_vectorL,*dd_vectorR; 
    double *xaxis, *yLvalues, *yRvalues;
    double ext_a, ext_b, dx;
    int i, j, N;
    int located=0;
    double middlex, v1, v2;

    //Find interval of coincidence
    if(fabs(waveExtremesL[MIN0(iL-1)][((zero_var+1)%2)])>fabs(waveExtremesR[MIN0(iR-1)][((zero_var+1)%2)])){
        if(sl<0 && sr<0) {ext_a=waveExtremesR[MIN0(iR-1)][((zero_var+1)%2)]; }
        else {ext_a=waveExtremesL[MIN0(iL-1)][((zero_var+1)%2)]; }
    }
    else if (fabs(waveExtremesL[MIN0(iL-1)][((zero_var+1)%2)])<fabs(waveExtremesR[MIN0(iR-1)][((zero_var+1)%2)])){
        if(sl<0 && sr<0) { ext_a=waveExtremesL[MIN0(iL-1)][((zero_var+1)%2)];  }
        else{ ext_a=waveExtremesR[MIN0(iR-1)][((zero_var+1)%2)];  }
    }
    else{
        ext_a=waveExtremesR[MIN0(iR-1)][((zero_var+1)%2)];
    }
    if(fabs(waveExtremesL[iL+1][((zero_var+1)%2)])<fabs(waveExtremesR[iR+1][((zero_var+1)%2)])){
        if(sl<0 && sr<0) {ext_b=waveExtremesR[iR+1][((zero_var+1)%2)]; }
        else {ext_b=waveExtremesL[iL+1][((zero_var+1)%2)];  }
    }
    else if(fabs(waveExtremesL[iL+1][((zero_var+1)%2)])>fabs(waveExtremesR[iR+1][((zero_var+1)%2)])){
        if(sl<0 && sr<0) { ext_b=waveExtremesL[iL+1][((zero_var+1)%2)];  }
        else {ext_b=waveExtremesR[iR+1][((zero_var+1)%2)]; }
    }
    else {
        ext_b=waveExtremesR[iR+1][((zero_var+1)%2)];
    }
    double DI=fabs(ext_a-ext_b);
    if (DI>1) DI=0.05/DI;
    N=(int)(fmin(fmax(50.0/DI,2000),8000)); //always want same precision
    printf("\nPossible intersection of wave curves in range: %e %.e\n",ext_a,ext_b);


    //declare memory: addecuate points for calculating the curve, given number for the interval
    if ((x_vectorL=(double*)malloc(N*sizeof(double)))==NULL) m_fail("memory for x interp.");
    if ((y_vectorL=(double*)malloc(N*sizeof(double)))==NULL) m_fail("memory for y interp.");
    if ((x_vectorR=(double*)malloc(N*sizeof(double)))==NULL) m_fail("memory for x interp.");
    if ((y_vectorR=(double*)malloc(N*sizeof(double)))==NULL) m_fail("memory for y interp.");
    if ((xaxis=(double*)malloc(N*sizeof(double)))==NULL) m_fail("memory for y interp.");
    if ((yLvalues=(double*)malloc(N*sizeof(double)))==NULL) m_fail("memory for y interp.");
    if ((yRvalues=(double*)malloc(N*sizeof(double)))==NULL) m_fail("memory for y interp.");

    double v,P,r,Pb,vstart;        
    //Obtain points of the curve L
    v=waveExtremesL[iL][0]; P=waveExtremesL[iL][1]; r=waveExtremesL[iL][2]; 
    if (iL-1>=0) {Pb=waveExtremesL[iL-1][1]; vstart=waveExtremesL[iL-1][0];}
    else {Pb=P; vstart=v; }

    if (wTL[(int)(iL*0.5)]==3){
        hugoniot_interpolate(v, P, r, waveExtremesL[iL+1][1], Pb, vstart, ext_a,ext_b, zero_var,  N, x_vectorL, y_vectorL,-1);
    }
    else if (wTL[(int)(iL*0.5)]==1){
        integral_interpolate(v, P, r, waveExtremesL[iL+1][2], ext_a, ext_b, zero_var, N, x_vectorL, y_vectorL,-1);
    }
    else if (wTL[(int)(iL*0.5)]==2){
        mixed_interpolate(v, P, r,  waveExtremesL[iL][3], ext_a, ext_b, zero_var,  N, x_vectorL, y_vectorL,-1);
    }

    //Obtain points of the curve R
    v=waveExtremesR[iR][0]; P=waveExtremesR[iR][1]; r=waveExtremesR[iR][2];
    if (iR-1>=0) {Pb=waveExtremesR[iR-1][1]; vstart=waveExtremesR[iR-1][0];}
    else {Pb=P; vstart=v; }

    if (wTR[(int)(iR*0.5)]==3){
        hugoniot_interpolate(v, P, r, waveExtremesR[iR+1][1], Pb, vstart, ext_a, ext_b, zero_var, N, x_vectorR, y_vectorR, 1);
    }
    else if (wTR[(int)(iR*0.5)]==1){
        integral_interpolate(v, P, r, waveExtremesR[iR+1][2], ext_a, ext_b, zero_var, N, x_vectorR, y_vectorR,1);
    }
    else if (wTR[(int)(iR*0.5)]==2){        
        mixed_interpolate(v, P, r, waveExtremesR[iR][3],  ext_a, ext_b, zero_var, N, x_vectorR, y_vectorR,1);
    }
    
    if(zero_var==0){ //interchange memmory addresses
        aux=y_vectorL;
        y_vectorL=x_vectorL;
        x_vectorL=aux;

        aux=y_vectorR;
        y_vectorR=x_vectorR;
        x_vectorR=aux;
    }

    //Create fix abscissas, equispaced
    dx=(ext_b-ext_a)/(N-1);
    for(i=0; i<N; i++) xaxis[i]=ext_a+i*dx;

    //look for points at each side of fixes x, create parabola and interpolate
    //L
    for(i=0;i<N;i++){
        for(j=0; j<N-1; j++){
            if((x_vectorL[j]<=xaxis[i] && x_vectorL[j+1]>=xaxis[i]) || (x_vectorL[j]>=xaxis[i] && x_vectorL[j+1]<=xaxis[i])){
                yLvalues[i]=parabola(x_vectorL, y_vectorL, j, xaxis[i], N); 
                break;
            }
        }
    }
    //R
    for(i=0;i<N;i++){
        for(j=0; j<N-1; j++){
            if((x_vectorR[j]<=xaxis[i] && x_vectorR[j+1]>=xaxis[i]) || (x_vectorR[j]>=xaxis[i] && x_vectorR[j+1]<=xaxis[i])){
                yRvalues[i]=parabola(x_vectorR, y_vectorR, j, xaxis[i], N); 
                break;
            }
        }
    }

    //now we compute the difference, as new ordinate. Recycle vector yLvalues.
    for(i=0; i<N; i++) yLvalues[i]=yLvalues[i]-yRvalues[i];    
    //look for a change of sign in the differences
    for (i=0; i<N-1; i++){
        if(yLvalues[i]*yLvalues[i+1]<0){
            located=1;
            break;
        }
    }
    if(!located) {
        free(x_vectorL); free(y_vectorL); free(x_vectorR); free(y_vectorR); 
        free(xaxis); free(yLvalues);  free(yRvalues);
        printf("There was no intersection\n");
        return 0;
    } //no change of difference detected, no intersection.

    //locate intersection, save the coordinate, and guess the paired value
    middlex=recta(yLvalues, xaxis, i, 0.0, N);
    printf("Intersection of wave curves found at x coordinate %e, ", middlex);

    for(j=0; j<N; j++){
        if((x_vectorL[j]<=middlex && x_vectorL[j+1]>=middlex) || (x_vectorL[j]>=middlex && x_vectorL[j+1]<=middlex)){
            v1=parabola(x_vectorL, y_vectorL, j, middlex, N); 
            break;
        }
    }
    for(j=0; j<N; j++){
        if((x_vectorR[j]<=middlex && x_vectorR[j+1]>=middlex) || (x_vectorR[j]>=middlex && x_vectorR[j+1]<=middlex)){
            v2=parabola(x_vectorR, y_vectorR, j, middlex, N); 
            break;
        }
    }
    
    printf("where y coordinate difference is %e\n", v1-v2);
    if(zero_var==0) {(*Pmiddle)=middlex; (*vmiddle)=0.5*(v1+v2); }
    else {(*vmiddle)=middlex; (*Pmiddle)=0.5*(v1+v2);}

    free(x_vectorL); free(y_vectorL); free(x_vectorR); free(y_vectorR); 
    free(xaxis); free(yLvalues);  free(yRvalues);
    return 1;

}

struct Data saveInfo(int nwl, int nwr, int *wtl, int *wtr, double **wel, double **wer, double vmiddle, double Pmiddle){
    struct Data sol;
    int i, count;
    double v, P, r, vw, rgoal, vr, Pr, rr, vwr;
    double g, css;

    //count how many actual waves there are, save their types
    count=0;
    for(i=0; i<nwl; i++){
        if (wtl[i]!=4) count++;
        else break;
    }
    sol.n_wl=count; count=0;
    for(i=0; i<nwr; i++){
        if (wtr[i]!=4) count++;
        else break;
    }
    sol.n_wr=count;
    if(NULL==(sol.wtl=(int*)malloc((sol.n_wl+1)*sizeof(int)))) m_fail("struct array");
    if(NULL==(sol.wtr=(int*)malloc((sol.n_wr+1)*sizeof(int)))) m_fail("struct array");
    memcpy(sol.wtl, wtl, sol.n_wl*sizeof(int)); sol.wtl[sol.n_wl]=4;
    memcpy(sol.wtr, wtr, sol.n_wr*sizeof(int)); sol.wtr[sol.n_wr]=4;
    printf("LEFT waves stored:");for(i=0;i<=sol.n_wl;i++) printf(" %d",sol.wtl[i]);printf("\n");
    printf("RIGHT waves stored:");for(i=0;i<=sol.n_wr;i++) printf(" %d",sol.wtr[i]);printf("\n\n");

    //information for the wave curves
    if(NULL==(sol.wcl=(double**)malloc((sol.n_wl*2)*sizeof(double*)))) m_fail("struct array");
    if(NULL==(sol.wcr=(double**)malloc((sol.n_wr*2)*sizeof(double*)))) m_fail("struct array");
    for(i=0; i<2*sol.n_wl; i++) {
        if(NULL==(sol.wcl[i]=(double*)malloc(4*sizeof(double)))) m_fail("struct array");
        memcpy(sol.wcl[i],wel[i],4*sizeof(double));
    }
    for(i=0; i<2*sol.n_wr; i++){
        if(NULL==(sol.wcr[i]=(double*)malloc(4*sizeof(double)))) m_fail("struct array");
        memcpy(sol.wcr[i],wer[i],4*sizeof(double));
    }

    //information for the exact solution
    if(NULL==(sol.esl=(double**)malloc((sol.n_wl*2)*sizeof(double*)))) m_fail("struct array");
    if(NULL==(sol.esr=(double**)malloc((sol.n_wr*2)*sizeof(double*)))) m_fail("struct array");
    for(i=0; i<2*sol.n_wl; i++) {
        if(NULL==(sol.esl[i]=(double*)malloc(4*sizeof(double)))) m_fail("struct array");
        memcpy(sol.esl[i],wel[i],4*sizeof(double));
    }
    for(i=0; i<2*sol.n_wr; i++){
        if(NULL==(sol.esr[i]=(double*)malloc(4*sizeof(double)))) m_fail("struct array");
        memcpy(sol.esr[i],wer[i],4*sizeof(double));
    }
    


    double Pb, aux1,aux2,aux3,aux4, aux5;
    for(i=1;i<2*(sol.n_wl-1);i=i+2){
        if(sol.wtl[(int)((i-1)*0.5)]==1){
            g=gmma(sol.esl[i][2]);
            css=sol.esl[i][1]*(g*(g-1)+sol.esl[i][2]*dgmma(sol.esl[i][2]))/((g-1)*sol.esl[i][2]+sol.esl[i][1]*g);
            sol.esl[i][3]=(sol.esl[i][0]-sqrt(css))/(1-sol.esl[i][0]*sqrt(css));
        }
        else if(sol.wtl[(int)((i-1)*0.5)]==2){    //if mixed is not last wave, we have to rearrange final density of the rarefaction to draw it
            v=wel[i-1][0]; P=wel[i-1][1]; r=wel[i-1][2]; 
            findRarefactionEnd_mixedsonic(v,P,r,wel[i-1][3], &aux1, &aux2, &aux3, &aux4,-1);
            sol.esl[i-1][3]=aux3;
        }
    }
    for(i=1;i<2*sol.n_wr;i=i+2){
        if(sol.wtr[(int)((i-1)*0.5)]==1){
            g=gmma(sol.esr[i][2]);
            css=sol.esr[i][1]*(g*(g-1)+sol.esr[i][2]*dgmma(sol.esr[i][2]))/((g-1)*sol.esr[i][2]+sol.esr[i][1]*g);
            sol.esr[i][3]=(sol.esr[i][0]+sqrt(css))/(1+sol.esr[i][0]*sqrt(css));
        }
        if(sol.wtr[(int)((i-1)*0.5)]==2){
            v=wer[i-1][0]; P=wer[i-1][1]; r=wer[i-1][2]; 
            findRarefactionEnd_mixedsonic(v,P,r,wer[i-1][3], &aux1, &aux2, &aux3, &aux4,1);
            sol.esr[i-1][3]=aux3;
        }
    }

    printf("LEFT side:\n");
    //Overwrite last wave for exact solution
    v=wel[2*sol.n_wl-2][0]; P=wel[2*sol.n_wl-2][1]; r=wel[2*sol.n_wl-2][2]; 
    if(2*sol.n_wl-3>=0) Pb=wel[2*sol.n_wl-3][1]; else Pb=P;

    vw=-2.0; //now we dont check for overtaking by stack
    if(sol.wtl[sol.n_wl-1]==3){ //if it is a shock
        build_shock(&v, &P, &r, &vw, -1, Pmiddle, Pb, &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        printf("Last shock. Obtained v=%e versus vmiddle=%e\n",v,vmiddle);
        sol.esl[2*sol.n_wl-1][0]=vmiddle;
        sol.esl[2*sol.n_wl-1][1]=Pmiddle;
        sol.esl[2*sol.n_wl-1][2]=r;
        sol.esl[2*sol.n_wl-1][3]=vw;  
    }
    else if(sol.wtl[sol.n_wl-1]==1){ //if it is a rarefaction
        rgoal=wel[2*sol.n_wl-1][2];
        find_rgoal(v, P, r, -1, Pmiddle, &rgoal);
        build_rarefaction(&v, &P, &r, &vw, -1, rgoal,  &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        printf("Last rarefaction. Obtained v=%e versus vmiddle=%e, P=%e versus Pmiddle=%e\n",v,vmiddle,P, Pmiddle);
        sol.esl[2*sol.n_wl-1][0]=vmiddle;
        sol.esl[2*sol.n_wl-1][1]=Pmiddle;
        sol.esl[2*sol.n_wl-1][2]=r;
        sol.esl[2*sol.n_wl-1][3]=vw;
    }
    else {//if it is a mixed curve
        rr=wel[2*sol.n_wl-2][3];
        last_mixed(&v, &P, &r, &vw, &vr, &Pr, &rr, &vwr, vmiddle, Pmiddle, -1);
        printf("Last mixed. Obtained v=%e versus vmiddle=%e, P=%e vs Pmiddle=%e\n",vr,vmiddle,Pr,Pmiddle);
        sol.esl[2*sol.n_wl-2][3]=r; //change the info to draw the rarefaction

        sol.esl[2*sol.n_wl-1][0]=vmiddle; //mixed curve at intersection
        sol.esl[2*sol.n_wl-1][1]=Pmiddle;
        sol.esl[2*sol.n_wl-1][2]=rr;
        sol.esl[2*sol.n_wl-1][3]=vwr;
    }

    printf("RIGHT side:\n");
    v=wer[2*sol.n_wr-2][0]; P=wer[2*sol.n_wr-2][1]; r=wer[2*sol.n_wr-2][2];
    if(2*sol.n_wr-3>=0) Pb=wer[2*sol.n_wr-3][1]; else Pb=P;

    vw=-2.0;
    if(sol.wtr[sol.n_wr-1]==3){ //if it is a shock
        build_shock(&v, &P, &r, &vw, 1, Pmiddle, Pb, &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        printf("Last shock. Obtained v=%e versus vmiddle=%e\n",v,vmiddle);
        sol.esr[2*sol.n_wr-1][0]=vmiddle; 
        sol.esr[2*sol.n_wr-1][1]=Pmiddle;
        sol.esr[2*sol.n_wr-1][2]=r;
        sol.esr[2*sol.n_wr-1][3]=vw;
        
    }
    else if(sol.wtr[sol.n_wr-1]==1){ //if it is a rarefaction
        rgoal=wer[2*sol.n_wr-1][2];
        find_rgoal(v, P, r, 1, Pmiddle, &rgoal);
        build_rarefaction(&v, &P, &r, &vw, 1, rgoal,  &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        printf("Last rarefaction. Obtained v=%e versus vmiddle=%e, P=%e versus Pmiddle=%e\n",v,vmiddle,P, Pmiddle);
        sol.esr[2*sol.n_wr-1][0]=vmiddle;
        sol.esr[2*sol.n_wr-1][1]=Pmiddle;
        sol.esr[2*sol.n_wr-1][2]=r;
        sol.esr[2*sol.n_wr-1][3]=vw;
    }
    else {//if it is a mixed curve
        rr=wer[2*sol.n_wr-2][3];
        last_mixed(&v, &P, &r, &vw, &vr, &Pr, &rr, &vwr, vmiddle, Pmiddle, 1);
        printf("Last mixed. Obtained v=%e versus vmiddle=%e, P=%e vs Pmiddle=%e\n",vr,vmiddle,Pr,Pmiddle);
        sol.esr[2*sol.n_wr-2][3]=r; //overwrite end of rarefaction

        sol.esr[2*sol.n_wr-1][0]=vmiddle; //mixed wave at intersection
        sol.esr[2*sol.n_wr-1][1]=Pmiddle;
        sol.esr[2*sol.n_wr-1][2]=rr;
        sol.esr[2*sol.n_wr-1][3]=vwr;
    }

    printf("\n");
    printf("Key states spatial domain. Waves moving LEFT\n");
    printf("v            P            rho          v_wave\n");
    for(i=0; i<2*sol.n_wl; i++) printf("%e %e %e %e\n",sol.esl[i][0],sol.esl[i][1],sol.esl[i][2], sol.esl[i][3]);

    printf("\nKey states spatial domain. Waves moving RIGHT\n");
    printf("v            P            rho          v_wave\n");
    for(i=0; i<2*sol.n_wr; i++) printf("%e %e %e %e\n",sol.esr[i][0],sol.esr[i][1],sol.esr[i][2], sol.esr[i][3]);
    return sol;
}

struct Data findWaveExtremes(){
    int codeL, codeR, zero_var; //wave codes and variable to find intersection of curves
    int n_wavesL=1, n_wavesR=1; //how many waves there are from each initial position
    double ini_dif; //the initial difference between the variables that equal at intersection
    double vstack, ov, oP, or, ovw, stackable;
    int istack, index;


    double **waveExtremesL, **waveExtremesR;
    int *waveTypeL, *waveTypeR; //waves memory
    if ((waveExtremesL=(double**)malloc(2*sizeof(double*)))==NULL) m_fail("Allocating memory for waveExtremesL");
    if ((waveExtremesR=(double**)malloc(2*sizeof(double*)))==NULL) m_fail("Allocating memory for waveExtremesR");
    if ((waveExtremesL[0]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("memory for waveExtremesL i");
    if ((waveExtremesR[0]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("memory for waveExtremesR i");
    if ((waveExtremesL[1]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("memory for waveExtremesL i");
    if ((waveExtremesR[1]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("memory for waveExtremesR i");
    if ((waveTypeL=(int*)malloc(sizeof(int)))==NULL) m_fail("Allocating memory for waveTypeL");
    if ((waveTypeR=(int*)malloc(sizeof(int)))==NULL) m_fail("Allocating memory for waveTypeR");
    
    //initialization
    incrementSign(); 
    firstWave(&codeL, &codeR); 
    problemType(&zero_var, &ini_dif);
    waveExtremesL[0][0]=VL; waveExtremesL[0][1]=PL; waveExtremesL[0][2]=RL; waveExtremesL[0][3]=0; waveExtremesL[0][4]=0;
    waveExtremesR[0][0]=VR; waveExtremesR[0][1]=PR; waveExtremesR[0][2]=RR; waveExtremesR[0][3]=0; waveExtremesR[0][4]=0;
    waveTypeL[0]=codeL; waveTypeR[0]=codeR;
    
    //Variable initialization
    double v,P,r, gr, auxv,auxP,auxvm, Pb, goalPL=userGoalPl, goalRL=userGoalRl, goalPR=userGoalPr, goalRR=userGoalRr; 
    int looking_indexL,looking_indexR, bool_intersectionFound=0, Lcalculated=0, Rcalculated=0;

    int stackUsedL=0, stackUsedR=0;
    double vmiddle, Pmiddle;
    double css,g;
    int i;
    
    //calculate waves from each side until intersection
    do{
        printf("\n------------------------------------\n\n");
        
        if(codeL==0 && codeR==0){ //if both types are zero and no intersecion detected, amplify goals
            goalPL=goalPL+sl*goalPL*0.5;
            goalRL=goalRL+sl*goalRL*0.5;
            codeL=waveTypeL[n_wavesL-1];
            goalPR=goalPR+sr*goalPR*0.5;
            goalRR=goalRR+sr*goalRR*0.5;
            codeR=waveTypeR[n_wavesR-1];
        } //this way if one is finished and the other is not, there's no need to calculate more of the finished one
        
        looking_indexL=1; looking_indexR=1;
        
        //do left side
        if(codeL){
            Lcalculated=1;
            for(i=2*n_wavesL-1-2*looking_indexL; i>0; i=i-2) if (waveExtremesL[i][4]>0) break;
            if(i<0) i=0;
            istack=i;
            if(istack) {
                vstack=waveExtremesL[istack][3]; 
                ov=waveExtremesL[istack-1][0];
                oP=waveExtremesL[istack-1][1];
                or=waveExtremesL[istack-1][2];
                ovw=waveExtremesL[istack-1][3];
            }
            else{
                vstack=-2.0;
                ov=-2.0;
                oP=-2.0;
                or=-2.0;
                ovw=-2.0;
            }
            index=2*n_wavesL-1;

            v=waveExtremesL[index-1][0]; P=waveExtremesL[index-1][1]; r=waveExtremesL[index-1][2];            
            if (index-2>=0) Pb=waveExtremesL[index-2][1]; else Pb=waveExtremesL[index-1][1];
            if(codeL==3){
                printf("LEFT:calculating Hugoniot curve with origin v=%e P=%e rho=%e\n",v,P,r);
                codeL=build_shock(&v, &P, &r, &vstack, -1, goalPL, Pb, &stackable, &ov, &oP, &or, &ovw, 1); 
            }
            else if (codeL==1){
                printf("LEFT:calculating integral curve with origin v=%e P=%e rho=%e\n",v,P,r);
                codeL=build_rarefaction(&v, &P, &r, &vstack, -1, goalRL, &stackable, &ov, &oP, &or, &ovw, 1); 
            }
            else if (codeL==2){
                if(stackUsedL) {findRarefactionEnd_mixedsonic(v,P,r,waveExtremesL[index-1][3], &auxv, &auxP, &gr,&auxvm,-1);} else {gr=waveExtremesL[index-2][2];}
                waveExtremesL[index-1][3]=gr;
                printf("LEFT:calculating mixed curve with rarefaction origin v=%e P=%e rho=%e, ending at rho=%e\n",v,P,r,gr);
                codeL=build_mixedcurve(&v, &P, &r, &vstack, -1, gr, &stackable, &ov, &oP, &or, &ovw);
            }
            
            waveExtremesL[index][0]=v; waveExtremesL[index][1]=P; waveExtremesL[index][2]=r; waveExtremesL[index][3]=vstack;
            
            if (fabs(stackable)>0.6) {stackUsedL=0; waveExtremesL[index][4]=stackable;}
            else{
                waveExtremesL[index][4]=-1.0; waveExtremesL[istack][4]=-1.0;
                stackUsedL=1;
                codeL=waveTypeL[(int)((istack-1)*0.5)];
                if (codeL==1) {if(n_wavesL==2) codeL=3; else codeL=2;} //a rarefaction does not overtake, a related mixed does
            }
            if (n_wavesL==1 && waveTypeL[0]==1) {g=gmma(waveExtremesL[index-1][2]); css=waveExtremesL[index-1][1]*(g*(g-1)+waveExtremesL[index-1][2]*dgmma(waveExtremesL[index-1][2]))/((g-1)*waveExtremesL[index-1][2]+waveExtremesL[index-1][1]*g);waveExtremesL[index][4]=1.0; waveExtremesL[index][3]=(waveExtremesL[index-1][0]-sqrt(css))/(1-waveExtremesL[index-1][0]*sqrt(css));}
            printf("done LEFT curve, ended v=%e P=%e rho=%e ",v,P,r); 
            if(codeL==1) printf("wave speed=%f. Next wave type: integral curve.\n", ovw);
            else if(codeL==2) printf("wave speed=%f. Next wave type: mixed curve.\n", vstack); 
            else if (codeL==3) printf("wave speed=%f. Next wave type: Hugoniot curve.\n", vstack);   
            else printf("Next wave type: none for the moment (curve not terminated).\n");
        }
        
        if (codeL) { //if the goal was not reached and there is another wave afterwards
            n_wavesL++;
            if((waveTypeL=(int*)realloc(waveTypeL,n_wavesL*sizeof(int)))==NULL) m_fail("realloc waveTypeL");
            waveTypeL[n_wavesL-1]=codeL;
            if((waveExtremesL=(double**)realloc(waveExtremesL,(n_wavesL*2)*sizeof(double*)))==NULL) m_fail("realloc waveExtremesL");
            if((waveExtremesL[index+1]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("extremes entry after realloc L");
            waveExtremesL[index+1][0]=ov; waveExtremesL[index+1][1]=oP; waveExtremesL[index+1][2]=or; waveExtremesL[index+1][3]=ovw;
            if((waveExtremesL[index+2]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("extremes entry after realloc L");
        }
        else { 
            looking_indexL=0;
        }
        if(Lcalculated){
            if (possibilityIntersection(index-1, -1, waveExtremesL, waveExtremesR, waveTypeL, waveTypeR, n_wavesL, n_wavesR, looking_indexL, looking_indexR, zero_var, ini_dif, &vmiddle, &Pmiddle))
                break;
        }
        else {
           if (possibilityIntersection(2*n_wavesL-2, -1, waveExtremesL, waveExtremesR, waveTypeL, waveTypeR, n_wavesL, n_wavesR, looking_indexL, looking_indexR, zero_var, ini_dif, &vmiddle, &Pmiddle))
                break; 
        }
        Lcalculated=0;

        printf("\n");
        //do right side
        if(codeR ){
            Rcalculated=1;
            for(i=2*n_wavesR-1-2*looking_indexR; i>0; i=i-2) if (waveExtremesR[i][4]>0) break;
            if(i<0) i=0;
            istack=i;
            if(istack) {
                vstack=waveExtremesR[istack][3]; 
                ov=waveExtremesR[istack-1][0];
                oP=waveExtremesR[istack-1][1];
                or=waveExtremesR[istack-1][2];
                ovw=waveExtremesR[istack-1][3];
            }
            else{
                vstack=-2.0;
                ov=-2.0;
                oP=-2.0;
                or=-2.0;
                ovw=-2.0;
            }
            index=2*n_wavesR-1;

            v=waveExtremesR[index-1][0]; P=waveExtremesR[index-1][1]; r=waveExtremesR[index-1][2];
            if (index-2>=0) Pb=waveExtremesR[index-2][1]; else Pb=P;
            
            if(codeR==3){
                printf("RIGHT:calculating Hugoniot curve with origin v=%e P=%e rho=%e\n",v,P,r);
                codeR=build_shock(&v, &P, &r, &vstack, 1, goalPR, Pb, &stackable, &ov, &oP, &or, &ovw, 1); 
            }
            else if (codeR==1){
                printf("RIGHT:calculating integral curve with origin v=%e P=%e rho=%e\n",v,P,r);
                codeR=build_rarefaction(&v, &P, &r, &vstack, 1, goalRR, &stackable, &ov, &oP, &or, &ovw, 1); 
            }
            else if (codeR==2){
                if(stackUsedR) {findRarefactionEnd_mixedsonic(v,P,r,waveExtremesR[index-1][3], & auxv, &auxP, &gr,&auxvm,1);} else {gr=waveExtremesR[index-2][2];}
                waveExtremesR[index-1][3]=gr;
                printf("RIGHT:calculating mixed curve with rarefaction origin v=%e P=%e rho=%e, ending at rho=%e\n",v,P,r,gr);
                codeR=build_mixedcurve(&v, &P, &r, &vstack, 1, gr, &stackable, &ov, &oP, &or, &ovw);
            }
            
            waveExtremesR[index][0]=v; waveExtremesR[index][1]=P; waveExtremesR[index][2]=r; waveExtremesR[index][3]=vstack;
            if (fabs(stackable)>0.6) {stackUsedR=0; waveExtremesR[index][4]=stackable;}
            else{
                stackUsedR=1;
                waveExtremesR[index][4]=-1.0; 
                waveExtremesR[istack][4]=-1.0;
                codeR=waveTypeR[(int)((istack-1)*0.5)];
                if (codeR==1) {if(n_wavesR==2) codeR=3; else codeR=2;}
            }      
            if (n_wavesR==1 && waveTypeR[0]==1) {g=gmma(waveExtremesR[index-1][2]); css=waveExtremesR[index-1][1]*(g*(g-1)+waveExtremesR[index-1][2]*dgmma(waveExtremesR[index-1][2]))/((g-1)*waveExtremesR[index-1][2]+waveExtremesR[index-1][1]*g);waveExtremesR[index][4]=1.0; waveExtremesR[index][3]=(waveExtremesR[index-1][0]+sqrt(css))/(1+waveExtremesR[index-1][0]*sqrt(css));}
            printf("done RIGHT curve, ended v=%e P=%e rho=%e ",v,P,r); 
            if(codeR==1) printf("wave speed=%f. Next wave type: integral curve.\n", ovw);
            else if(codeR==2) printf("wave speed=%f. Next wave type: mixed curve.\n", vstack); 
            else if (codeR==3) printf("wave speed=%f. Next wave type: Hugoniot curve.\n", vstack);   
            else printf("Next wave type: none for the moment (curve not terminated).\n");
        }
        
        if (codeR) { //if the goal was not reached and there is another wave afterwards
            n_wavesR++;
            if((waveTypeR=(int*)realloc(waveTypeR,n_wavesR*sizeof(int)))==NULL) m_fail("realloc waveTypeR");
            waveTypeR[n_wavesR-1]=codeR;
            if((waveExtremesR=(double**)realloc(waveExtremesR,(n_wavesR*2)*sizeof(double*)))==NULL) m_fail("realloc waveExtremesR");
            if((waveExtremesR[index+1]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("extremes entry after realloc R");
            waveExtremesR[index+1][0]=ov; waveExtremesR[index+1][1]=oP; waveExtremesR[index+1][2]=or; waveExtremesR[index+1][3]=ovw;
            if((waveExtremesR[index+2]=(double*)malloc(5*sizeof(double)))==NULL) m_fail("extremes entry after realloc R");        }
        else { 
            looking_indexR=0;
        }

      if(Rcalculated){
            if (possibilityIntersection(index-1, 1, waveExtremesL, waveExtremesR, waveTypeL, waveTypeR, n_wavesL, n_wavesR, looking_indexL, looking_indexR, zero_var, ini_dif, &vmiddle, &Pmiddle))
                break;
      }
      else{
          if (possibilityIntersection(2*n_wavesR-2, 1, waveExtremesL, waveExtremesR, waveTypeL, waveTypeR, n_wavesL, n_wavesR, looking_indexL, looking_indexR, zero_var, ini_dif, &vmiddle, &Pmiddle))
                break;
      }
        Rcalculated=0;

    } while (!bool_intersectionFound); //looop to calculate wave curves until intersection
    
    //put solution together for output    
    struct Data solution;
    printf("\n------------------------------------\n");
    printf("-----------Saving solution----------\n");
    printf("------------------------------------\n");
    printf("Wave codes: 1:integral curve, 2:mixed curve, 3:Hugoniot curve (4: end of sequence)\n\n");

    solution=saveInfo(n_wavesL, n_wavesR, waveTypeL, waveTypeR, waveExtremesL, waveExtremesR, vmiddle, Pmiddle);
    printf("\nSolution information saved\n\n");
    free(waveTypeL); free(waveTypeR); for(codeL=0; codeL<2*n_wavesL; codeL++) free(waveExtremesL[codeL]);
    free(waveExtremesL); for(codeL=0; codeL<n_wavesR*2; codeL++) free(waveExtremesR[codeL]); free(waveExtremesR);
    //return solution
    return solution;
}

void exportWaveCurves(struct Data sol){
    int i, code;
    char letter;
    double v, P, r, vwave=-2.0, Pb, aux1, aux2, aux3, aux4, aux5;
    savewc=1;
    fprintf(fwc, "#1-v 2-P 3-rho 4-v_wave 5-lambda 6-nonl_factor\n");

    for(i=0; i<sol.n_wl; i++){
        vwave=-2.0;
        code=sol.wtl[i];
        if(code==1) letter='I'; else if (code==2) letter='M'; else if (code==3) letter='H'; else letter='W';
        fprintf(fwc, "\n\n\"L-%c\"\n",letter);

        v=sol.wcl[i*2][0]; P=sol.wcl[i*2][1];  r=sol.wcl[i*2][2];
        if(code==1){
            build_rarefaction(&v, &P, &r, &vwave, -1, sol.wcl[i*2+1][2], &aux1, &aux2, &aux3, &aux4, &aux5, 0); 
        }
        if(code==3){
            if(i*2-1>0) Pb=sol.wcl[i*2-1][1]; else Pb=P;
            build_shock(&v, &P, &r, &vwave, -1, sol.wcl[i*2+1][1], Pb, &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        }
        if(code==2){
            build_mixedcurve(&v, &P, &r, &vwave, -1, sol.wcl[i*2][3], &aux1, &aux2, &aux3, &aux4, &aux5);
        }
            
    }
    
    for(i=0; i<sol.n_wr; i++){
        vwave=-2.0;
        code=sol.wtr[i];
        if(code==1) letter='I'; else if (code==2) letter='M'; else if (code==3) letter='H'; else letter='W';
        fprintf(fwc, "\n\n\"R-%c\"\n",letter);

        v=sol.wcr[i*2][0]; P=sol.wcr[i*2][1];  r=sol.wcr[i*2][2];
        if(code==1){
            build_rarefaction(&v, &P, &r, &vwave, 1, sol.wcr[i*2+1][2], &aux1, &aux2, &aux3, &aux4, &aux5, 0); 
        }
        if(code==3){
            if(i*2-1>0) Pb=sol.wcr[i*2-1][1]; else Pb=P;
            build_shock(&v, &P, &r, &vwave, 1, sol.wcr[i*2+1][1], Pb, &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        }
        if(code==2){
            build_mixedcurve(&v, &P, &r, &vwave, 1, sol.wcr[i*2][3], &aux1, &aux2, &aux3, &aux4, &aux5);
        }
    }
    savewc=0;
}

void exportExactSolution(struct Data sol){
    int i;
    double *wavespeedsL, *wavespeedsR;
    int *orderL, *orderR;
    if(NULL==(wavespeedsL=(double*)malloc(sol.n_wl*sizeof(double)))) m_fail("exact");
    if(NULL==(wavespeedsR=(double*)malloc(sol.n_wr*sizeof(double)))) m_fail("exact");
    if(NULL==(orderL=(int*)malloc(sol.n_wl*sizeof(int)))) m_fail("exact");
    if(NULL==(orderR=(int*)malloc(sol.n_wr*sizeof(int)))) m_fail("exact");
    
    for(i=1; i<2*sol.n_wl; i=i+2) wavespeedsL[(int)((i-1)*0.5)]=sol.esl[i][3];
    purgeL(wavespeedsL, orderL, sol.n_wl); //select which waves appear in the spatial domain
    for(i=1; i<2*sol.n_wr; i=i+2) wavespeedsR[(int)((i-1)*0.5)]=sol.esr[i][3];
    purgeR(wavespeedsR, orderR, sol.n_wr); 

    printf("Index of waves appearing in the spatial domain:\n");
    printf("Left: ");for(i=0; i<sol.n_wl; i++) printf("%d ",orderL[i]);printf("\n");
    printf("Right: ");for(i=0; i<sol.n_wr; i++) printf("%d ",orderR[i]);printf("\n");

    double v, P, r, css, g, vw,  prevv,prevP,prevr,  aux1,aux2,aux3,aux4,aux5, l0;
    savees=1;
    fprintf(fes, "#1-x 2-v 3-P 4-rho\n");

    g=gmma(RL);
    css=PL*(g*(g-1)+RL*dgmma(RL))/((RL*(g-1))*(1+PL/RL*g/(g-1)));
    prevv=VL; prevP=PL; prevr=RL;
    fprintf(fes, "0.00000000 %.8e %.8e %.8e \n", VL, PL, RL);

    for(i=0; i<sol.n_wl; i++){
        if (orderL[i]<0) break;

        if(sol.wtl[orderL[i]]==3){
            v=sol.esl[orderL[i]*2+1][0]; P=sol.esl[orderL[i]*2+1][1]; r=sol.esl[orderL[i]*2+1][2]; vw=sol.esl[orderL[i]*2+1][3];
            g=gmma(r);
            css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
            fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw-EPSILON,
                 prevv, prevP, prevr);
            fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw,
                 v, P, r);
        }
        else if(sol.wtl[orderL[i]]==2){
            v=sol.esl[orderL[i]*2][0]; P=sol.esl[orderL[i]*2][1]; r=sol.esl[orderL[i]*2][2];
            g=gmma(r);
            css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
            l0=(v-sqrt(css))/(1-v*sqrt(css));
            if (l0<sol.esl[orderL[i]*2+1][3]){
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*l0-EPSILON,
                 prevv, prevP, prevr); 
                build_rarefaction(&v,&P, &r, &vw, -1, sol.esl[orderL[i]*2][3], &aux1, &aux2, &aux3, &aux4, &aux5, 0);
                v=sol.esl[orderL[i]*2+1][0]; P=sol.esl[orderL[i]*2+1][1]; r=sol.esl[orderL[i]*2+1][2]; vw=sol.esl[orderL[i]*2+1][3];
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw,
                 v, P, r);
            }
            else {
                v=sol.esl[orderL[i]*2+1][0]; P=sol.esl[orderL[i]*2+1][1]; r=sol.esl[orderL[i]*2+1][2]; vw=sol.esl[orderL[i]*2+1][3];
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw-EPSILON,
                 prevv, prevP, prevr);
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw,
                 v, P, r);
            }
        }
        else{
            v=sol.esl[orderL[i]*2][0]; P=sol.esl[orderL[i]*2][1]; r=sol.esl[orderL[i]*2][2];
            g=gmma(r);
            css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
            l0=(v-sqrt(css))/(1-v*sqrt(css));
            fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*l0-EPSILON,
                 prevv, prevP, prevr);
            build_rarefaction(&v,&P, &r, &vw, -1, sol.esl[orderL[i]*2+1][2], &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        }
        prevv=v; prevP=P; prevr=r;  
    }

    //contact discontinuity from the left
    v=sol.esl[2*orderL[i-1]+1][0]; P=sol.esl[orderL[i-1]*2+1][1]; r=sol.esl[orderL[i-1]*2+1][2]; //por asegurar consistencia decimales
    g=gmma(r);
    css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
    double fx,fv,fP,fr;
    fx=x_ini_disc+tf*v; fv=v; fP=P; fr=r;
    fprintf(fes, "%8f %.8e %.8e %.8e\n", x_ini_disc+tf*v, v, P, r);

    fprintf(fes,"\n\n");

    g=gmma(RR);
    css=PR*(g*(g-1)+RR*dgmma(RR))/((RR*(g-1))*(1+PR/RR*g/(g-1)));
    prevv=VR; prevP=PR; prevr=RR;
    fprintf(fes, "1.00000000 %.8e %.8e %.8e\n", VR, PR, RR);
    for(i=0; i<sol.n_wr; i++){
        if (orderR[i]<0) break;

        if(sol.wtr[orderR[i]]==3){
            v=sol.esr[orderR[i]*2+1][0]; P=sol.esr[orderR[i]*2+1][1]; r=sol.esr[orderR[i]*2+1][2]; vw=sol.esr[orderR[i]*2+1][3];
            g=gmma(r);
            css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
            fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw+EPSILON,
                 prevv, prevP, prevr);
            fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw,
                 v, P, r);
        }
        else if(sol.wtr[orderR[i]]==2){
            v=sol.esr[orderR[i]*2][0]; P=sol.esr[orderR[i]*2][1]; r=sol.esr[orderR[i]*2][2];
            g=gmma(r);
            css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
            l0=(v+sqrt(css))/(1+v*sqrt(css));
            if (l0>sol.esr[orderR[i]*2+1][3]){
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*l0+EPSILON,
                 prevv, prevP, prevr); 
                build_rarefaction(&v,&P, &r, &vw, 1, sol.esr[orderR[i]*2][3], &aux1, &aux2, &aux3, &aux4, &aux5, 0);
                v=sol.esr[orderR[i]*2+1][0]; P=sol.esr[orderR[i]*2+1][1]; r=sol.esr[orderR[i]*2+1][2]; vw=sol.esr[orderR[i]*2+1][3];
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw,
                 v, P, r);
            }
            else {
                v=sol.esr[orderR[i]*2+1][0]; P=sol.esr[orderR[i]*2+1][1]; r=sol.esr[orderR[i]*2+1][2]; vw=sol.esr[orderR[i]*2+1][3];
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw+EPSILON,
                 prevv, prevP, prevr);
                fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*vw,
                 v, P, r);
            }
        }
        else{
            v=sol.esr[orderR[i]*2][0]; P=sol.esr[orderR[i]*2][1]; r=sol.esr[orderR[i]*2][2];
            g=gmma(r);
            css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
            l0=(v+sqrt(css))/(1+v*sqrt(css));
            fprintf(fes, "%.8e %.8e %.8e %.8e\n", x_ini_disc+tf*l0+EPSILON,
                 prevv, prevP, prevr);
            build_rarefaction(&v,&P, &r, &vw, 1, sol.esr[orderR[i]*2+1][2], &aux1, &aux2, &aux3, &aux4, &aux5, 0);
        }
        prevv=v; prevP=P; prevr=r; 
    }
    
    //contact discontinuty from the right
    v=sol.esr[orderR[i-1]*2+1][0]; P=sol.esr[orderR[i-1]*2+1][1]; r=sol.esr[orderR[i-1]*2+1][2]; //por asegurar consistencia decimales
    g=gmma(r);
    css=P*(g*(g-1)+r*dgmma(r))/((r*(g-1))*(1+P/r*g/(g-1)));
    fprintf(fes, "%8e %.8e %.8e %.8e\n", x_ini_disc+tf*v, v, P, r);

    savees=0;

    fprintf(fes, "\n\n");
    fprintf(fes, "%8e %.8e %.8e %.8e\n", x_ini_disc+tf*v, v, P, r);
    fprintf(fes, "%8e %.8e %.8e %.8e\n", fx, fv, fP, fr);

    free(wavespeedsL); free(wavespeedsR); free(orderL); free(orderR);
}
