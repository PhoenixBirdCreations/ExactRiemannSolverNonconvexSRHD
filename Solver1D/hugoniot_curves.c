#include "hugoniot_curves.h"

double hb_Taub(double Pa, double Pb, double rb, double ha, double ra){
    double root=(Pb-Pa)*(Pb-Pa)/(rb*rb)+4*(ha*ha+ha/ra*(Pb-Pa));
    root=sqrt(root);
    double num=(Pb-Pa)/rb+root;
    num*=0.5;
    return num;
}

double F_hug_newton(double x, double Pb, double Pa, double ha, double ra){
/* 
 function for Newton method. It is the definition of the pressure using the EOS, introducing the enthalpy obtained from the Taub adiabat.
*/
    double g=gmma(x);
    return -Pb+(g-1)*x/g*(hb_Taub(Pa,Pb,x,ha,ra)-1);
}

double dF_hug_newton(double x, double Pb, double Pa, double ha, double ra){
/*
 derivative of the function for Newton method
*/
    double g=gmma(x);
    double dg=dgmma(x);
    return (g*(g-1)+dg*x)/g/g*(-1+0.5*((Pb-Pa)/x+sqrt((Pb-Pa)*(Pb-Pa)/x/x+4*(ha*ha+ha/ra*(Pb-Pa)))))+(g-1)*x/g*0.5*((Pa-Pb)/x/x+0.5*(-2*(Pb-Pa)*(Pb-Pa)/x/x/x))/sqrt((Pb-Pa)*(Pb-Pa)/x/x+4*(ha*ha+ha/ra*(Pb-Pa)));
}

double newton_hugoniot(double prev, double Pa, double ha, double ra, double Pb){
/*
 Newton method to find the density corresponding to some pressure in a Hugoniot curve 
*/
    int i;
    double x,xnew,xmemory;
    x=prev;
    for(i=0;i<MAXIT;i++){
        xnew=x-F_hug_newton(x,Pb,Pa,ha,ra)/dF_hug_newton(x,Pb,Pa,ha,ra);
        if(fabs(F_hug_newton(xnew,Pb,Pa,ha,ra))<TOL) break;
        xmemory=x;
        x=xnew;
    }
    if (i>=MAXIT) {
        printf("\nERROR:newton_hugoniot did not converge F=%.15f x_k+1=%f x_k=%f\n",F_hug_newton(xnew,Pb,Pa,ha,ra), x, xmemory);
        return -1.0;
    }
    return xnew;
}

double jsquare(double Pa, double ra, double ha, double Pb, double rb, double hb){
/*
 square of the mass flux invariant across the shock
*/
    double j;
    j=rb*ra*(Pa-Pb)/(hb*ra-ha*rb);
    if(j<0) {printf("\nERROR:jsquare is negative!! Pa=%f Pb=%f ra=%f rb=%f\n",Pa,Pb,ra,rb); exit(1);}
    return j;
}

double vshock(double ra, double va, double js, int s){
/*
 shock speed, s=-1 to left, s=1 to right
*/
    double wa=1.0/sqrt(1-va*va);
    return (ra*ra*wa*wa*va+s*js*sqrt(1+ra*ra/js))/(ra*ra*wa*wa+js);
}

double vflow(double ha, double va, double Pa, double ra, double vs, double Pb, double j){
/* 
 flow speed given by the Rankine-Hugoniot jump conditions
*/
    double ws=1.0/sqrt(1-vs*vs);
    double wa=1.0/sqrt(1-va*va);
    double dem;
    dem=ha*wa+(Pb-Pa)*(ws*va/j+1.0/(ra*wa));
    return (ha*wa*va+ws*(Pb-Pa)/j)/dem;
}

int build_shock(double *v, double *P, double *r, double *vw, int s, double goal_P, double Pb, double *stackable, double *ov, double *oP, double *or, double *ovw, int first){
    //Initialization
    double va, Pa, ra, g, ha;
    va=*v; Pa=*P; ra=*r;
    g=gmma(ra);
    ha=1.0+Pa/ra*g/(g-1);
    double vs, vb, rb, hb, js, css, dP;
    rb=ra;
    double vs_m, r_m, vdif=0.0, difStack=0.0;
    vs_m=0.0;
    int count=0,sonic=0;
    int sense;
    if(s<0) sense=sl; else sense=sr;
    if (first) dP=sense*0.0005; else dP=(goal_P-Pb)/5000;
    //Calculate shock
    while(sense*Pb<sense*goal_P){
        count++;
        if(count==1) difStack=vs-(*vw);
        Pb=Pb+dP; if(fabs(Pb-goal_P)<TOL || sense*Pb>sense*goal_P) Pb=goal_P;
        rb=newton_hugoniot(rb, Pa, ha, ra, Pb);
        g=gmma(rb);
        hb=1.0+Pb*g/(g-1.0)/rb; if(hb<0){printf("\nERROR:Not valid hb=%f\n",hb); exit(0);}
        js=jsquare(Pa, ra, ha, Pb, rb, hb); 
        vs=vshock(ra, va, js, s);
        vb=vflow(ha,va,Pa,ra,vs,Pb,s*sqrt(js));
        
        //check monotony v_shock       
        if ((vs-(*vw))*difStack<0){ //termination by stack has preference
            if(dP<TOL) {sonic=2; printf("--Hugoniot curve reached by stack\n");break;}
            else {Pb-=dP; dP*=0.5;}
        }
        else if((vs_m-vs)*vdif<0){  //then sonic point
            if(dP<TOL) {sonic=1; printf("--Hugoniot curve terminated for admissibility. v_s^(i-1)-v_s^(i):%e, vs'=%e\n",vs_m-vs,(vs_m-vs)/(rb-r_m));break;}
            else {Pb-=dP; dP*=0.5;}
        }
        else {  //else continue
            if(count>1) vdif=vs_m-vs;

            //prepare next iteration
            vs_m=vs; r_m=rb;

            //exporting output wave curves
            if (savewc){
                css=Pb*(g*(g-1)+rb*dgmma(rb))/((g-1)*rb+Pb*g);
                fprintf(fwc, "%.8e %.8e %.8e %.8e %.8e %.8e\n", vb, Pb, rb, vs,(vb+s*sqrt(css))/(1+s*vb*sqrt(css)) ,nonlfactor(css,rb));
            }
        }
    }
    
    (*r)=rb; (*P)=Pb; (*v)=vb; 
    
    if(fabs(vs-(*vw))<1.0e-5) sonic=2;

    (*vw)=vs;

    if(sonic==1){ //breaks because becomes sonic. Up to stack, and origin of next
        (*stackable)=1.0;
        (*ov)=vb; (*oP)=Pb; (*or)=rb; (*ovw)=vs;
        return 1; 
    }
    else if (sonic==2){ //breaks because overtaken by stack
        (*stackable)=-0.5;
        return 1;
    }
    else{
        (*stackable)=-1.0;
        (*ov)=vb; (*oP)=Pb; (*or)=rb; (*ovw)=vs;
        return 0;
    }
}

void hugoniot_interpolate(double va, double Pa, double ra, double Pgoal, double Pb, double vstart, double exa, double exb, int zero_var, int N, double *xv, double *yv, int s){
    int entering=0, exiting=0;                          
    double Penter, Pexit, difv=0.0;
    
    if (zero_var==0){ //gas slabs, extremes in P
        if (fabs(Pb-exa)<TOL){ //=, but floating point tolerance
            Penter=Pb;
            entering=1;
        }
    }
    else{
        if (fabs(vstart-exa)<TOL){ //=, but floating point tolerance
            Penter=Pb;
            entering=1;
            difv=exb-vstart;
        }
        else difv=exa-vstart;
    }
    
    int sense;
    if(s<0) sense=sl; else sense=sr;
    
    double g, ha;
    g=gmma(ra);
    ha=1+Pa/ra*g/(g-1);
    
    double vs, vb, rb, hb, js, dP, Pstart=Pb, P_m;
    rb=ra;
    int i;
    dP=(Pgoal-Pstart)/(N-1);
    P_m=Pstart;
    
    //Calculate shock
    for(i=1;i<N;i++){
        Pb=Pstart+i*dP;
        rb=newton_hugoniot(rb, Pa, ha, ra, Pb);
        g=gmma(rb);
        hb=1.0+Pb*g/(g-1.0)/rb; if(hb<0){printf("\nERROR:Not valid hb=%f\n",hb); exit(0);}
        js=jsquare(Pa, ra, ha, Pb, rb, hb); if (js<0) {printf("\nERROR:jsquared is negative\n"); exit(0);}
        vs=vshock(ra, va, js, s);
        vb=vflow(ha,va,Pa,ra,vs,Pb,s*sqrt(js));
        if (!entering){ //interval of intersection NOT reached yet
            if(zero_var==1){ //case blast wave, looking at pressure
                if(difv*(exa-vb)<0){ //interval reached
                    entering=1;
                    Penter=P_m; 
                }
            }
            else{
                if(sense*Pb>sense*exa){ //interval reached
                    entering=1;
                    Penter=P_m;
                }
            }
        }
        else if (entering && !exiting){
            if(zero_var==1){ //case blast wave, looking at pressure
                if(difv*(exb-vb)<0){ //interval reached
                    exiting=1;
                    Pexit=Pb;
                }
            }
            else{
                if(sense*Pb>sense*exb){ //interval reached
                    exiting=1;
                    Pexit=Pb;
                }
            }
        }
        P_m=Pb;
        if (entering && exiting) break;

    }
    if(!exiting) Pexit=Pb;
    dP=(Pexit-Penter)/(N-1);
    //Calculate shock
    for(i=0;i<N;i++){
        Pb=Penter+i*dP;
        rb=newton_hugoniot(rb, Pa, ha, ra, Pb);
        g=gmma(rb);
        hb=1.0+Pb*g/(g-1.0)/rb; if(hb<0){printf("\nERROR:Not valid hb=%f\n",hb); exit(0);}
        js=jsquare(Pa, ra, ha, Pb, rb, hb); if (js<0) {printf("\nERROR:jsquared is negative\n"); exit(0);}
        vs=vshock(ra, va, js, s);
        vb=vflow(ha,va,Pa,ra,vs,Pb,s*sqrt(js));
        xv[i]=vb; yv[i]=Pb;
    } 

}







