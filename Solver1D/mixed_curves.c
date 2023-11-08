#include "mixed_curves.h"

double v_rho(double l, double D, double vr, double r){
    double M=D*(vr-l);
    return (l*r*r+M*sqrt(r*r-l*l*r*r+M*M))/(r*r+M*M);
}

double dv_rho(double D, double l, double vr, double r){
    double M=D*(vr-l), raiz=sqrt(1-l*l+M*M/(r*r));
    return -M*(M*M*(l*l+1)+r*r*(1-l*l)-2*M*l*r*raiz)/(raiz*(M*M+r*r)*(M*M+r*r));
}

double eq3(double D, double vr, double Pr, double rr, double l, double r){
    double v=v_rho(l,D,vr,r);
    double w=1.0/sqrt(1-v*v);
    double g, hr;
    g=gmma(r);
    hr=1+Pr/rr*gmma(rr)/(gmma(rr)-1);
    return w*(1-l*v)*(g*Pr+(g-1)*r)+hr/sqrt(1-vr*vr)*(g*w*D*(vr-l)*(vr-v)-(g-1)*r*(1-l*vr));
}

double deq3(double D, double vr, double Pr, double rr, double l, double r){
    double v, w, dw, dv, g, dg;
    v=v_rho(l,D,vr,r); dv=dv_rho(D,l,vr,r);
    w=1.0/sqrt(1-v*v); dw=w*w*w*v*dv;
    g=gmma(r); dg=dgmma(r);
    double hr=1+Pr/rr*gmma(rr)/(gmma(rr)-1);

    return (g*Pr+(g-1)*r)*(dw*(1-l*v)-w*l*dv)+w*(1-l*v)*(dg*Pr+g-1+r*dg)+hr/sqrt(1-vr*vr)*(D*(vr-l)*((vr-v)*(dg*w+g*dw)-g*w*dv)-(1-l*vr)*(g-1+dg*r));
    
}

double delicateBisection(double D, double vr, double P, double rr, double l, double xnew, int sense){
 double ba, bb, bx, sa, sb, sx;
    double step=0.01;
    int count=0;
    sa=eq3(D,vr,P,rr,l,cap(xnew+sense*0.01));
    sa=-sa;

    bx=xnew;
    do{
        ba=bx;
        bx=bx-sense*step;
        sb=eq3(D,vr,P,rr,l,bx);
        if(sense*bx<=sense*rr){
            step/=100;
            if(step<1.0e-14){printf("\nERROR:step too small in mixed curve\n"); exit(1);}
            bx=xnew;
        }  
        else if(sa*sb<0 && eq3(D,vr,P,rr,l,ba)*sb>0) {step/=100; bx=xnew; sb=sa;}
    }while(sa*sb>0);
    bb=bx; 
    sa=eq3(D,vr,P,rr,l,ba);
    sb=eq3(D,vr,P,rr,l,bb); 
    count=0;
    
    if(sa*sb>0) {
        if (fabs(sa)<TOL && fabs(sb)<TOL) return ba;
        printf("\nERROR:Same sign extremes in bisection for mixed curves\n"); exit(1);}
    if (fabs(sa)<TOL) {
        bx=ba-sense*step/10;
        sx=eq3(D,vr,P,rr,l,bx);
        while(sx*sb>0){
            step*=0.1;
            bx=ba-sense*step/10;
            sx=eq3(D,vr,P,rr,l,bx);
        }
        ba=bx; sa=sx;
    }
    do{
        count++;
        bx=0.5*(ba+bb);
        sx=eq3(D,vr,P,rr,l,bx);
        if(fabs(sx)<TOL) {break;}
        else{
            if(sx*sa<=0) bb=bx;
            else {sa=sx; ba=bx;}
        }
    }while(count<2*MAXIT);
    if(count>=2*MAXIT) {
        if (fabs(sx)>1.0e-10) {
            printf("\nERROR:Bisection for mixed curves did not converge\n"); 
            exit(1);
        }
    }
    return bx;
}

double newton_eq(double rr, double D, double vr, double P, double l, int sense, double prev){
    double F, dF, x, xnew, alpha;
    //printf("\ndouble D=%.15f, vr=%.15f, Pr=%.15f, l=%.15f, rr=%.15f;\n",D,vr,P,l,rr);
    xnew=rr;
    alpha=deq3(D,vr,P,rr,l,rr+sense*EPSILON);
    int count=0, outercount=1;
    do{
        if(sense>0) {if(prev<0) xnew=rr*0.5*(outercount+3); else xnew=prev*0.5*(outercount+2);}
        else {if (outercount<5){ if(prev<0) xnew=rr*(1-0.1*outercount); else xnew=prev*(1-0.1*outercount);}else {if(prev<0) xnew=rr/(outercount+1); else xnew=prev/(outercount+1);}}
        outercount++;
        count=0;
        do{
            x=xnew;
            F=eq3(D,vr,P,rr,l,x);
            if(fabs(F)<TOL) {break;} //tolerance for zero
            dF=deq3(D,vr,P,rr,l,x); 
            xnew=fabs(x-F/dF); //because a negative value is not physical
            count++;
        }while(fabs(x-xnew)>1e-14 && count<MAXIT); //Newton convergence

        if(count>=MAXIT) { //Newton get lost in steep derivatives
            if(fabs(eq3(D,vr,P,rr,l,xnew))>1e-10 && fabs(xnew-rr)<0.0001)
                return -1.0;
        }
        
        else if(prev<0){ //for the first point of the system, make sure we are not in trivial solution
            if (fabs(rr-xnew)>0.00001)
                break;
        }

        else if (outercount>100){ //too many x0 tried. There's no other solution than the trivial
            return -1.0;
        }

        else if (fabs(prev-xnew)>0.5){ //went to a root very far away
            xnew=rr;
        }

    }while(fabs(rr-xnew)<fabs(prev-xnew)); //try with other starting point if wrong basin of attraction

    if(xnew<0){printf("\nERROR:negative density in Newton mixed curves\n"); exit(1);}

    if(fabs(F)>1.0e-6) return -1.0;
    if(deq3(D,vr,P,rr,l,xnew)*alpha>0) { //if the curve has three solutions, and we have the wrong one
        xnew=delicateBisection(D,vr,P,rr,l,xnew,sense);
    }
    if (prev>0 && sense*prev>sense*xnew) { //something weird, wrong direction
        return -1.0;
    }
    return xnew;
}

double eps_rho(double vr, double P, double rr, double l, double r, double v){
    double hr=1+P/rr*gmma(rr)/(gmma(rr)-1);
    double fw=sqrt(1-v*v)/sqrt(1-vr*vr);
    return (fw*hr*((1-l*vr)/(1-l*v))-1)/gmma(r);
}

void store_rarefaction(double **rare, double va, double Pa, double ra, double goal_r, int s, int N){
    int i,j;
    double X, vb=va, Pb=Pa, drho; 
    double  *R, *eps;
    R=(double*) malloc(7*sizeof(double));
    eps=(double*) malloc(4*sizeof(double));
    drho=(goal_r-ra)/N;
    eps[0]=Pa/((gmma(ra)-1)*ra);
    //calculate the rarefaction using N points and store it
    for (i=0;i<N;i++){
        R[0]=ra+i*drho; R[6]=ra+(i+1)*drho;
        R[1]=(5*R[0]+R[6])/6.0; R[2]=(2*R[0]+R[6])/3.0; R[3]=R[2];
        R[4]=(R[0]+2*R[6])/3.0; R[5]=(R[0]+R[6])*0.5;
        rare[i][0]=vb; rare[i][1]=Pb; rare[i][2]=R[0];
    
        for(j=1;j<4;j++){
            eps[j]=eps[0]*exp((R[2*j]-R[0])/6.0*(int_e_isen(R[0])+4*int_e_isen(R[2*j-1])+int_e_isen(R[2*j])));
        }
        X=(R[6]-R[0])/8.0*(int_R_inv(R[0],eps[0])+3*int_R_inv(R[2],eps[1])+3*int_R_inv(R[4],eps[2])+int_R_inv(R[6],eps[3]));
        
        vb=((1+va)/(1-va)*exp(s*2*X)-1)/((1+va)/(1-va)*exp(s*2*X)+1);
        Pb=(gmma(R[6])-1)*R[6]*eps[3];
        
        eps[0]=eps[3];
        va=vb;
    }
    rare[N][0]=vb; rare[N][1]=Pb; rare[N][2]=R[6];
    free(R); free(eps);
}

int build_mixedcurve(double *vini, double *Pini, double *rini, double *vm, int s, double rb, double *st, double *ov, double *oP, double *or, double *ovw){
    double **rare;
    int i, N=1000, sonic=0;
    double va,Pa,ra, linw, difStack;
    va=(*vini); Pa=(*Pini); ra=(*rini); //difStack=va-(*vm);
    rare=(double**) malloc((N+1)*sizeof(double*));
    for(i=0;i<=N;i++) rare[i]=(double*) malloc(3*sizeof(double));
    
    store_rarefaction(rare, va, Pa, ra, rb, s, N );
    //prepare for solving the system of equations
    double vr, Pr, rr, wr, D, g, css, l, m_dif=0.0;
    Pr=rare[N-1][1]; vr=rare[N-1][0]; rr=rare[N-1][2];
    g=gmma(rr);        
    css=Pr*(g*(g-1)+rr*dgmma(rr))/((g-1)*rr+Pr*g);
    l=(vr+s*sqrt(css))/(1+s*vr*sqrt(css));
    difStack=l-(*vm);

    double v,r, P;
    double prev=-2;
    int count, sense; 
    if(s<0) sense=sl; else sense=sr;
    for(i=N-1;i>=0;i--){ //the solution for the last point of the rarefaction is itself 
        Pr=rare[i][1]; vr=rare[i][0]; rr=rare[i][2];
        g=gmma(rr);        
        css=Pr*(g*(g-1)+rr*dgmma(rr))/((g-1)*rr+Pr*g);
        l=(vr+s*sqrt(css))/(1+s*vr*sqrt(css));
        wr=1.0/sqrt(1-vr*vr);
        D=rr*wr; 
        r=newton_eq(rr, D, vr, Pr, l, sense, prev);
        if (r<0){
            if(fabs(rare[i][2]-rare[i+1][2])<1.0e-12){
               if(fabs(linw-l)<1.0e-5){
                sonic=1; r=prev;
                break;  
               }
               else{
                printf("\nERROR:mixed curve without solution and can't refine more. Last time solved l*=%.9f l=%.9f\n",l,linw);
                exit(1);
               }
            } 
            //refining rarefaction
            store_rarefaction(rare, rare[i][0], rare[i][1], rare[i][2], rare[i+1][2], s, N );
            i=N;
        }
        else{
            v=v_rho(l,D,vr,r);
            g=gmma(r);
            P=(g-1)*r*eps_rho(vr,Pr,rr,l,r,v); 
            css=P*(g*(g-1)+r*dgmma(r))/((g-1)*r+P*g);
            linw=(v+s*sqrt(css))/(1+s*v*sqrt(css));

            if((l-(*vm))*difStack<0){
                if(fabs(l-(*vm))<TOL){
                    sonic=2;
                    break;
                }
                else{
                    if (fabs(rare[i][2]- rare[i+1][2])<1.0e-12){
                        if(fabs(v-(*vm))<1.0e-5){
                            sonic=2;
                            break;
                        }
                    }
                    else{
                        store_rarefaction(rare, rare[i][0], rare[i][1], rare[i][2], rare[i+1][2], s, N );
                        i=N;
                    }
                }

            }

            if(fabs(linw-l)<1.0e-6 && i<N-2){ //depending on the problem, it cannot handle more accuracy
                printf("--Mixed curve ending at lambda-lambda*= %e\n",linw-l);
                sonic=1; 
                break;
            }

            else if((l-linw)*m_dif<0){
                if (fabs(rare[i][2]- rare[i+1][2])<1.0e-12){
                    if(fabs(linw-l)<1.0e-5){
                        printf("--Mixed curve ending at lambda-lambda*= %e\n",linw-l);
                        sonic=1;
                        break;
                    }
                    else{
                        printf("\nERROR:mixed curve can't refine more end point, lambda-lambda*= %e\n",linw-l);
                        exit(1); 
                    }
                }
                store_rarefaction(rare, rare[i][0], rare[i][1], rare[i][2], rare[i+1][2], s, N );
                i=N;
            }
            else m_dif=(l-linw);

            prev=r;

            //Guardamos el resultado
            if (savewc){          
                fprintf(fwc, "%.8e %.8e %.8e %.8e %.8e %.8e\n", v, P, r, l,linw, nonlfactor(css,r) );    
            }
        }
    
    }

    if (savewc){
        v=v_rho(l, D, vr, r);
        P=(gmma(r)-1)*r*eps_rho(vr,Pr,rr,l,r,v);            
        fprintf(fwc, "%.8e %.8e %.8e %.8e %.8e %.8e\n", v, P, r, l,linw, nonlfactor(css,r) );    
    
    }
    for(count=0;count<=N;count++) {free(rare[count]);} free(rare);

    (*vini)=v;
    (*Pini)=P;
    (*rini)=r;
    
    if(fabs(l-(*vm))<1.0e-5) sonic=2; //has done all rarefaction, coincides with stack speed

    (*vm)=l;

    if (sonic==1){ //breaks because sonic
        printf("--Mixed curve terminated\n");
        (*st)=1.0;
        (*ov)=v;
        (*oP)=P;
        (*or)=r;
        (*ovw)=l;
        return 1;
    }
    else if (sonic==2){  //overtaken by stack
    printf("--Mixed curve reaches stack speed\n");
        (*st)=-0.5;
        return 1;
    }
    else{ //covers all rarefaction
        printf("--Mixed uses all integral curve\n");
        (*st)=-1.0;
        (*ov)=v;
        (*oP)=P;
        (*or)=r;
        (*ovw)=l;
        return 3;
    }
}

void mixed_interpolate(double va, double Pa, double ra, double rb, double exa, double exb, int zero_var, int N, double *xv, double *yv, int s){
    double **rare;
    int i;
    
    rare=(double**) malloc((N+1)*sizeof(double*));
    for(i=0;i<N+1;i++) rare[i]=(double*) malloc(3*sizeof(double)); 
    store_rarefaction(rare, va, Pa, ra, rb, s, N );


    int entering=0, exiting=0;
    double renter, vexit, Pexit, rexit, difv;
    if (zero_var==0){ //gas slabs, extremes in P
        if (fabs(rare[N][1]-exa)<TOL){ //=, but floating point tolerance
            renter=rare[N][2];
            entering=1;
        }
    }
    else{
        if (fabs(rare[N][0]-exa)<TOL){ //=, but floating point tolerance
            renter=rare[N][2];
            entering=1;
            difv=exb-rare[N][0];
        }
        else difv=exa-rare[N][0];
    }

    //prepare for solving the system of equations
    double vr, Pr, rr, wr, D, g, css, l, m_dif=0.0, linw;
    double v,r,P, prev=-2;
    int sense; 
    if(s<0) sense=sl; else sense=sr;
    for(i=N-1;i>=0;i--){ 
        Pr=rare[i][1]; vr=rare[i][0]; rr=rare[i][2];
        g=gmma(rr);        
        css=Pr*(g*(g-1)+rr*dgmma(rr))/((g-1)*rr+Pr*g);
        l=(vr+s*sqrt(css))/(1+s*vr*sqrt(css));
        wr=1.0/sqrt(1-vr*vr);
        D=rr*wr;

        r=newton_eq(rr, D, vr, Pr, l, sense, prev); 
        if (r<0){
            if(fabs(rare[i][2]-rare[i+1][2])<1.0e-12){
                printf("\nERROR:No solution and can't refine more. Last time solved l*=%.9f l=%.9f\n",l,(v+s*sqrt(css))/(1+s*v*sqrt(css)));
                exit(1);
            } 
            store_rarefaction(rare, rare[i][0], rare[i][1], rare[i][2], rare[i+1][2], s, N );
            i=N;
        }
        else{
            prev=r;
            v=v_rho(l,D,vr,r);
            g=gmma(r);
            P=(g-1)*r*eps_rho(vr,Pr,rr,l,r,v); 

            if (!entering){ //interval of intersection NOT reached yet
                if(zero_var==1){ //case blast wave, looking at pressure
                    if(difv*(exa-v)<0){ //interval reached
                        entering=1;
                        renter=rr;
                    }
            }
            else{
                if(sense*P>sense*exa){ //interval reached
                    entering=1;
                    renter=rr;
                }
            }
            }
            else if (entering && !exiting){
                if(zero_var==1){ //case blast wave, looking at pressure
                    if(difv*(exb-v)<0){ //interval reached
                        exiting=1;
                        vexit=vr; Pexit=Pr; rexit=rr;
                    }
                }
                else{
                    if(sense*P>sense*exb){ //interval reached
                        exiting=1;
                        vexit=vr; Pexit=Pr; rexit=rr;
                    }
                }
            }

            if (entering && exiting) break;

            css=P*(g*(g-1)+r*dgmma(r))/((g-1)*r+P*g);
            linw=(v+s*sqrt(css))/(1+s*v*sqrt(css));
            
            if(fabs(linw-l)<1.0e-10){
                break;
            }

            else if((l-linw)*m_dif<0){
                if (fabs(rare[i][2]- rare[i+1][2])<1.0e-12){
                    if(fabs(linw-l)<1.0e-5){
                        break;
                    }
                    else{
                        printf("\nERROR:interpolating mixed curve can't refine more sonic value, l-l*=%e\n",l-linw);
                        exit(1); 
                    }
                }
                store_rarefaction(rare, rare[i][0], rare[i][1], rare[i][2], rare[i+1][2], s, N );
                i=N;
            }

            else m_dif=(l-linw);  
        }      
    }
    if (!exiting) {vexit=vr; Pexit=Pr; rexit=rr;}
    store_rarefaction(rare, vexit, Pexit, rexit, renter, s, N ); 
    prev=-2;
    for(i=N-1;i>=0;i--){
        Pr=rare[i][1]; vr=rare[i][0]; rr=rare[i][2];
        g=gmma(rr);        
        css=Pr*(g*(g-1)+rr*dgmma(rr))/((g-1)*rr+Pr*g);
        l=(vr+s*sqrt(css))/(1+s*vr*sqrt(css));
        wr=1.0/sqrt(1-vr*vr);
        D=rr*wr;

        r=newton_eq(rr, D, vr, Pr, l, sense, prev); 
        v=v_rho(l,D,vr,r);
        g=gmma(r);
        P=(g-1)*r*eps_rho(vr,Pr,rr,l,r,v); 
        xv[N-1-i]=v; yv[N-1-i]=P;
        prev=r;
    }

for(i=0;i<=N;i++) {free(rare[i]);} free(rare);
}

void findRarefactionEnd_mixedsonic(double va, double Pa, double ra, double rb, double *vout, double *Pout, double *rout, double *vm, int s){
    double **rare;
    int i, N=1000;
    
    rare=(double**) malloc((N+1)*sizeof(double*));
    for(i=0;i<=N;i++) rare[i]=(double*) malloc(3*sizeof(double));
    
    store_rarefaction(rare, va, Pa, ra, rb, s, N );

    //prepare for solving the system of equations
    double vr, Pr, rr, wr, D, g, css, l, m_dif=0.0, linw;
    double v,r, P, prev=-2;
    int count, sense; 
    if(s<0) sense=sl; else sense=sr;
    for(i=N-1;i>=0;i--){ 
        Pr=rare[i][1]; vr=rare[i][0]; rr=rare[i][2];
        g=gmma(rr);        
        css=Pr*(g*(g-1)+rr*dgmma(rr))/((g-1)*rr+Pr*g);
        l=(vr+s*sqrt(css))/(1+s*vr*sqrt(css));
        wr=1.0/sqrt(1-vr*vr);
        D=rr*wr; 

        r=newton_eq(rr, D, vr, Pr, l, sense, prev);
        if (r<0){
            if(fabs(rare[i][2]-rare[i+1][2])<1.0e-12){
                printf("\nERROR:mixed curve without solution and can't refine more. Last time solved l-l*=%e\n",linw-l);
                exit(1);
            } 
            store_rarefaction(rare, rare[i][0], rare[i][1], rare[i][2], rare[i+1][2], s, N );
            i=N;
        }
        else{
            prev=r;
            v=v_rho(l,D,vr,r);
            g=gmma(r);
            P=(g-1)*r*eps_rho(vr,Pr,rr,l,r,v); 
            css=P*(g*(g-1)+r*dgmma(r))/((g-1)*r+P*g);
            linw=(v+s*sqrt(css))/(1+s*v*sqrt(css));
            
            if(fabs(linw-l)<1.0e-6){ //no admite mayor tol según el problema
                break;
            }

            else if((l-linw)*m_dif<0){
                if (fabs(rare[i][2]- rare[i+1][2])<TOL){
                    printf("\nERROR:mixed curve without solution and can't refine more. Last time solved l-l*=%e\n",linw-l);
                    exit(1); 
                }
                store_rarefaction(rare, rare[i][0], rare[i][1], rare[i][2], rare[i+1][2], s, N );
                i=N;
            }

            else m_dif=(l-linw);    

            if (savewc){
                v=v_rho(l, D, vr, r);
                P=(gmma(r)-1)*r*eps_rho(vr,Pr,rr,l,r,v);            
                fprintf(fwc, "%.8e %.8e %.8e\n", v, P, r);    
            }
        }
    }
    if(i<0) i=0;
    (*Pout)=rare[i][1]; (*vout)=rare[i][0]; (*rout)=rare[i][2]; (*vm)=(v+s*sqrt(css))/(1+s*v*sqrt(css));
    if (savewc){
        v=v_rho(l, D, vr, r);
        P=(gmma(r)-1)*r*eps_rho(vr,Pr,rr,l,r,v);            
        fprintf(fwc, "%.8e %.8e %.8e\n", v, P, r);    
    }
    for(count=0;count<=N;count++) {free(rare[count]);} free(rare);
}

void last_mixed(double *vr, double *Pr, double *rr, double *vwr, double *v, double *P, double *r, double *vw, double vmiddle, double Pmiddle, int s){
    double **rare;
    int i, N=1000;
    double va, Pa, ra, rb;
    va=*vr; Pa=*Pr; ra=*rr; rb=*r;
    rare=(double**) malloc((N+1)*sizeof(double*));
    for(i=0;i<N+1;i++) rare[i]=(double*) malloc(3*sizeof(double)); 
    store_rarefaction(rare, va, Pa, ra, rb, s, N );
    //prepare for solving the system of equations
    double cvr, cPr, crr, wr, D, g, css, l, difP=0.0;
    double cv,cr,cP, prev=-2;
    int sense; 
    if(s<0) sense=sl; else sense=sr;

    for(i=N-1;i>=0;i--){ //the solution for the last point of the rarefaction is itself (?)
        cPr=rare[i][1]; cvr=rare[i][0]; crr=rare[i][2];
        g=gmma(crr);        
        css=cPr*(g*(g-1)+crr*dgmma(crr))/((g-1)*crr+cPr*g);
        l=(cvr+s*sqrt(css))/(1+s*cvr*sqrt(css));
        wr=1.0/sqrt(1-cvr*cvr);
        D=crr*wr;

        cr=newton_eq(crr, D, cvr, cPr, l, sense, prev);
        prev=cr;
        cv=v_rho(l,D,cvr,cr);
        g=gmma(cr);
        cP=(g-1)*cr*eps_rho(cvr,cPr,crr,l,cr,cv); 

        if(difP*(Pmiddle-cP)<0){
            if(fabs(Pmiddle-cP)<TOL){
                break;
            }
            else{
                if(fabs(rare[i][2]-rare[i+1][2])<1.0e-12) break; //can't refine more
                store_rarefaction(rare,rare[i][0],rare[i][1],rare[i][2], rare[i+1][2],s,N);
                i=N;
            }
        }
        else{
            difP=Pmiddle-cP;
        }

    }

    (*vr)=cvr; (*Pr)=cPr; (*rr)=crr; (*vwr)=l;
    (*v)=cv; (*P)=cP; (*r)=cr; (*vw)=l;

    for(i=0;i<=N;i++) {free(rare[i]);} free(rare);
}

