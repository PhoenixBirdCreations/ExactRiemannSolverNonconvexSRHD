#include "integral_curves.h"

double int_e_isen(double r){
/* 
 the function to integrate to obtain internal energy e, isentropic flow dU=-PdV 
*/   
    return (gmma(r)-1)/r;
}

double int_R_inv(double r, double eps){
/* 
 the function to integrate to calculate the Riemann invariant between two points 
*/    
    double g=gmma(r);
    return sqrt(eps*(g*(g-1)+r*dgmma(r))/(1+eps*g))/r;
}

int build_rarefaction(double *v, double *P, double *r, double *vm, int s, double goal_r, double *stackable, double *ov, double *oP, double *or, double *ovw, int first){
/* 
 double integration algorithm that describes the rarefaction. Taking the initial point _a, it creates an array of desired densities to traverse. Then it calculates the velocity of the fluid at that point using the Riemann invariant. To do so, it also calculates the internal energy of the state with the desired density using that the flow is isentropic during rarefactions. The function ends successfully if the goal density is reached. If the convexity fails, it calculates the exact point and return the state in it.
 v, P, r are input (initial states) and output (nonl=0 state). s is input, the direction of the wave: -1 left, 1 right. goal_r is input desired final density. N is input number of points for discretization.
 return 0=goal reached, return=1 nonl=0
 */
    //Initialization
    double va, Pa, ra, ea;
    va=(*v); Pa=(*P); ra=(*r);
    ea=Pa/((gmma(ra)-1)*ra);
    
    double *R, *eps;
    R=(double*) malloc(7*sizeof(double)); if(!R) m_fail("build_rarefaction, R");
    eps=(double*) malloc(4*sizeof(double)); if(!eps) m_fail("build_rarefaction, eps");
    
    int count,j, undef=0;
    double drho;
    int sense;
    if(s<0) sense=sl; else sense=sr;
    if (first) drho=sense*0.0005; else drho=(goal_r-ra)/10000;
    eps[0]=ea;
    
    
    double X, vb, Pb, css, l, nu, g;
    double nu_m, e_m, r_m;
    
    e_m=ea; r_m=ra;
    nu_m=nonlfactor(ea*(gmma(ra)*(gmma(ra)-1)+ra*dgmma(ra))/(1+ea*gmma(ra)),ra);
    R[0]=ra;
    //start calculus of rarefaction
    while(sense*R[0]<sense*goal_r){
        count++;
        //middle points for integration
        R[6]=R[0]+drho; if (fabs(R[6]-goal_r)<TOL) R[6]=goal_r;
        R[1]=(5*R[0]+R[6])/6.0; R[2]=(2*R[0]+R[6])/3.0; R[3]=R[2];
        R[4]=(R[0]+2*R[6])/3.0; R[5]=(R[0]+R[6])*0.5;
        //calculate e at new densities
        for(j=1;j<4;j++){
            eps[j]=eps[0]*exp((R[2*j]-R[0])/6.0*(int_e_isen(R[0])+4*int_e_isen(R[2*j-1])+int_e_isen(R[2*j])));
        }
        
        //check convexity
        g=gmma(R[6]);
        Pb=(g-1)*R[6]*eps[3];
        css=Pb*(g*(g-1)+R[6]*dgmma(R[6]))/((g-1)*R[6]+Pb*g);
        nu=nonlfactor(css,R[6]);
        
        if(first && nu*nu_m<0.0 && !undef) {
            find_nonl0(r_m, e_m, R[6], eps[3], r, P);
            drho=(*r)-R[0];
            undef=1;
        }
        
        else{
            //integral of Riemann invariant 
            X=(R[6]-R[0])/8.0*(int_R_inv(R[0],eps[0])+3*int_R_inv(R[2],eps[1])+3*int_R_inv(R[4],eps[2])+int_R_inv(R[6],eps[3]));
            
            //use to calculate flow speed
            vb=((1+va)/(1-va)*exp(s*2*X)-1)/((1+va)/(1-va)*exp(s*2*X)+1); 
            
            //prepare next iteration and save
            eps[0]=eps[3];
            va=vb;
            nu_m=nu; e_m=eps[3]; r_m=R[6];
            R[0]=R[6];
            
            if (savewc){
                l=(vb+s*sqrt(css))/(1+s*vb*sqrt(css));
                fprintf(fwc, "%.8e %.8e %.8e %.8e %.8e %.8e\n", vb, Pb, R[6],l, l,nonlfactor(css,R[6]));
            }
            if (savees){
                l=(vb+s*sqrt(css))/(1+s*vb*sqrt(css));
                fprintf(fes, "%.8e %.8e %.8e %.8e\n",x_ini_disc+tf*l, vb, Pb, R[6]);
            }

            if(undef) {printf("--Integral curve terminated at nonlinearity term: %e\n",nu); break;}
        }
    }
    l=(vb+s*sqrt(css))/(1+s*vb*sqrt(css));
    (*or)=ra; (*oP)=Pa;  (*ov)=(*v);
    (*r)=R[6];
    (*P)=Pb; 
    (*v)=vb;
    (*vm)=l; 
    (*ovw)=R[6];
    (*stackable)=-1;
    free(eps); free(R);
    if (undef){
        return 2;
    }
    else{
        return 0;   
    } 

}

void find_nonl0(double ra, double ea, double rb, double eb, double *r, double *P){
    double na, nx, f_ea, f_ra, x, eps;
    f_ea=ea; f_ra=ra;
    na=nonlfactor(ea*(gmma(ra)*(gmma(ra)-1)+ra*dgmma(ra))/(1+ea*gmma(ra)),ra); 
    int count=0;
    while(count<MAXIT){
        count++;
        x=0.5*(ra+rb);
        eps=f_ea*exp((x-f_ra)/6.0*(int_e_isen(f_ra)+4.0*int_e_isen(0.5*(f_ra+x))+int_e_isen(x)));
        nx=nonlfactor(eps*(gmma(x)*(gmma(x)-1)+x*dgmma(x))/(1+eps*gmma(x)),x);
        if(fabs(nx)<TOL){
            (*r)=x;
            (*P)=eps*x*(gmma(x)-1);
            break;
        }
        else if (na*nx<0){
            rb=x;
            eb=eps;
        }
        else{
            ra=x;
            ea=eps;
            na=nonlfactor(ea*(gmma(ra)*(gmma(ra)-1)+ra*dgmma(ra))/(1+ea*gmma(ra)),ra);
        }
    }   
    if (count==MAXIT) {printf("\nERROR:Bisection looking for nu=0 did not converge\n"); exit(2);}
}

void integral_interpolate(double va, double Pa, double ra, double rb, double exa, double exb, int zero_var, int N, double *xv, double *yv, int s){
    int entering=0, exiting=0;
    double venter, Penter, renter, rexit, difv=0.0;
    
    if (zero_var==0){ //gas slabs, extremes in P
        if (fabs(Pa-exa)<TOL){ //=, but floating point tolerance
            venter=va; Penter=Pa; renter=ra;
            entering=1;
        }
    }
    else{
        if (fabs(va-exa)<TOL){ //=, but floating point tolerance
            venter=va; Penter=Pa; renter=ra;
            entering=1;
            difv=exb-va;
        }
        else difv=exa-va;
    }

    double *R, *eps;
    R=(double*) malloc(7*sizeof(double)); if(!R) m_fail("build_rarefaction, R");
    eps=(double*) malloc(4*sizeof(double)); if(!eps) m_fail("build_rarefaction, eps");
    int i,j;
    int sense;
    if (s<0) sense=sl;
    else sense=sr;
    double drho=(rb-ra)/(N-1);
    eps[0]=Pa/((gmma(ra)-1)*ra);
    double X, vb, Pb, P_m;
    R[0]=ra;
    //start calculus of rarefaction
    for(i=0;i<N-1;i++){
        //middle points for integration
        R[6]=R[0]+drho;
        R[1]=(5*R[0]+R[6])/6.0; R[2]=(2*R[0]+R[6])/3.0; R[3]=R[2];
        R[4]=(R[0]+2*R[6])/3.0; R[5]=(R[0]+R[6])*0.5;
        
        //calculate e at new densities
        for(j=1;j<4;j++){
            eps[j]=eps[0]*exp((R[2*j]-R[0])/6.0*(int_e_isen(R[0])+4*int_e_isen(R[2*j-1])+int_e_isen(R[2*j])));
        }
        X=(R[6]-R[0])/8.0*(int_R_inv(R[0],eps[0])+3*int_R_inv(R[2],eps[1])+3*int_R_inv(R[4],eps[2])+int_R_inv(R[6],eps[3]));
        vb=((1+va)/(1-va)*exp(s*2*X)-1)/((1+va)/(1-va)*exp(s*2*X)+1); 
        Pb=(gmma(R[6])-1)*R[6]*eps[3];

        if (!entering){ //interval of intersection NOT reached yet
            if(zero_var==1){ //case blast wave, looking at pressure
                if(difv*(exa-vb)<0){ //interval reached
                    entering=1;
                    venter=va; Penter=P_m; renter=R[0];
                }
            }
            else{
                if(sense*Pb>sense*exa){ //interval reached
                    entering=1;
                    venter=va; Penter=P_m; renter=R[0];
                }
            }
        }
        else if (entering && !exiting){
            if(zero_var==1){ //case blast wave, looking at pressure
                if(difv*(exb-vb)<0){ //interval reached
                    exiting=1;
                    rexit=R[6];
                }
            }
            else{
                if(sense*Pb>sense*exb){ //interval reached
                    exiting=1;
                    rexit=R[6];
                }
            }
        }
        
        //prepare next iteration
        eps[0]=eps[3];
        va=vb;
        P_m=Pb;
        R[0]=R[6];
        if (entering && exiting) break;
    }
    if (!exiting) rexit=R[6];
    
    xv[0]=venter; yv[0]=Penter;
    drho=(rexit-renter)/(N-1); 
    R[0]=renter; eps[0]=Penter/(renter*(gmma(renter)-1)); va=venter;
    for(i=0;i<N-1;i++){
        //middle points for integration
        R[6]=R[0]+drho;
        R[1]=(5*R[0]+R[6])/6.0; R[2]=(2*R[0]+R[6])/3.0; R[3]=R[2];
        R[4]=(R[0]+2*R[6])/3.0; R[5]=(R[0]+R[6])*0.5;
        //calculate e at new densities
        for(j=1;j<4;j++){
            eps[j]=eps[0]*exp((R[2*j]-R[0])/6.0*(int_e_isen(R[0])+4*int_e_isen(R[2*j-1])+int_e_isen(R[2*j])));
        }
        X=(R[6]-R[0])/8.0*(int_R_inv(R[0],eps[0])+3*int_R_inv(R[2],eps[1])+3*int_R_inv(R[4],eps[2])+int_R_inv(R[6],eps[3]));
        vb=((1+va)/(1-va)*exp(s*2*X)-1)/((1+va)/(1-va)*exp(s*2*X)+1); 
        Pb=(gmma(R[6])-1)*R[6]*eps[3];

        //prepare next iteration
        eps[0]=eps[3];
        va=vb;
        P_m=Pb;
        R[0]=R[6];
        xv[i+1]=vb; yv[i+1]=Pb;
    }
    
    free(eps); free(R);
}

void find_rgoal(double v, double P, double r, int s, double Pmiddle, double *rgoal){
    //Initialization
    double va, Pa, ra, ea;
    va=v; Pa=P; ra=r;
    ea=Pa/((gmma(ra)-1)*ra);
    
    double *R, *eps;
    R=(double*) malloc(7*sizeof(double)); if(!R) m_fail("build_rarefaction, R");
    eps=(double*) malloc(4*sizeof(double)); if(!eps) m_fail("build_rarefaction, eps");
    
    int j;
    double drho;
    int sense;
    if(s<0) sense=sl; else sense=sr;
    drho=(*rgoal-ra)/1000;
    eps[0]=ea;
    
    double X, vb, Pb, g;
    
    R[0]=ra;
    while(sense*R[0]<sense*(*rgoal)){
        R[6]=R[0]+drho;
        R[1]=(5*R[0]+R[6])/6.0; R[2]=(2*R[0]+R[6])/3.0; R[3]=R[2];
        R[4]=(R[0]+2*R[6])/3.0; R[5]=(R[0]+R[6])*0.5;
        
        for(j=1;j<4;j++){
            eps[j]=eps[0]*exp((R[2*j]-R[0])/6.0*(int_e_isen(R[0])+4*int_e_isen(R[2*j-1])+int_e_isen(R[2*j])));
        }
        
        g=gmma(R[6]);
        Pb=(g-1)*R[6]*eps[3];
        if(fabs(Pb-Pmiddle)<TOL) {*rgoal=R[6]; break;}
        if (sense*Pb>sense*Pmiddle){ //goes beyond
            Pb=(gmma(R[4])-1)*R[4]*eps[2];
            if (sense*Pb<sense*Pmiddle){
                if(fabs(Pb-Pmiddle)<TOL) {*rgoal=R[4]; break;}
                drho=R[4]-R[0];
            }
            else if (sense*(gmma(R[2])-1)*R[2]*eps[1]<sense*Pmiddle){
                if(fabs((gmma(R[2])-1)*R[2]*eps[1]-Pmiddle)<TOL) {*rgoal=R[2]; break;}
                drho=R[2]-R[0];
            }
            else drho=0.5*(R[2]-R[0]);
        }
        else{
            //integral of Riemann invariant 
            X=(R[6]-R[0])/8.0*(int_R_inv(R[0],eps[0])+3*int_R_inv(R[2],eps[1])+3*int_R_inv(R[4],eps[2])+int_R_inv(R[6],eps[3]));
            
            //use to calculate flow speed
            vb=((1+va)/(1-va)*exp(s*2*X)-1)/((1+va)/(1-va)*exp(s*2*X)+1); 
            
            //prepare next iteration
            eps[0]=eps[3];
            va=vb;
            R[0]=R[6];
        }
        
    }
    
    free(eps); free(R);
}
