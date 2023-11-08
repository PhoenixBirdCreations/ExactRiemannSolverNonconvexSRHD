#include "auxFun.h"


void m_fail(char*s){
    printf("Error allocating memory: %s\n",s);
    exit(1);
}

void bigToSmall(double *v, int *out, int N){
    int i,j;
    double *aux, var;
    aux=(double*)malloc(N*sizeof(double));
    memcpy(aux, v, N*sizeof(double));
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if(aux[i]<aux[j]){
                var=aux[i];
                aux[i]=aux[j];
                aux[j]=var;
            }
        }
    }
    for(i=0;i<N;i++){
        for(j=N-1;j>=0;j--){
            if (v[j]==aux[i]) {out[i]=j; break;}
        }
    }
    free(aux);
}

void smallToBig(double *v, int *out, int N){
    int i,j;
    double *aux, var;
    aux=(double*)malloc(N*sizeof(double));
    memcpy(aux, v, N*sizeof(double));
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if(aux[i]>aux[j]){
                var=aux[i];
                aux[i]=aux[j];
                aux[j]=var;
            }
        }
    }
    for(i=0;i<N;i++){
        for(j=N-1;j>=0;j--){
            if (v[j]==aux[i]) {out[i]=j; break;}
        }
    }
    free(aux);
}

void orderArray(double *v, int N){
    double *aux;
    if((aux=(double*)malloc(N*sizeof(double)))==NULL) {m_fail("ordering array");}
    memcpy(aux,v,N*sizeof(double));
    int i;
    for(i=0;i<N;i++) v[i]=aux[N-1-i];
    free(aux);
}

void purgeL(double *s, int *out, int N){
   int i,j,k=0;
   for (i=0; i<N; i++){ //look at all the waves
       for(j=i+1;j<N;j++){//Q1.is there a faster one than me later?
           if (s[j]<s[i]) {//A1.yes
               i=j-1; //move to that wave
               break;
            } 
        }
        if(j>=N){//A1.no (ended the search)
            out[k]=i;
            k++;
        }
    }
   if(k!=N){ //draw less waves than the initial amount
        out[k]=-1;
    }
}

void purgeR(double *s,  int *out, int N){
   int i,j,k=0;
   for (i=0; i<N; i++){ //look at all the waves
       for(j=i+1;j<N;j++){//Q1.is there a faster one than me later?
           if (s[j]>s[i]) {//A1.yes
               i=j-1; //move to that wave
               break;
            } 
        }
        if(j>=N){//A1.no (ended the search)
            out[k]=i;
            k++;
        }
    }
   if(k!=N){ //draw less waves than the initial amount
        out[k]=-1;
    }
}

/* interpolation */
int minsign(double x, double y){
    if (fabs(x)<=fabs(y)) return SIGN(x);
    else return SIGN(y);
}
double powereno(double x, double y){
    if (fabs(x)<TOL && fabs(y)<TOL) 
        return 0.0;
    else 
        return minsign(x,y)*fmin(fabs(x),fabs(y))*(x*x+y*y+2*fmax(fabs(x),fabs(y))*fmax(fabs(x),fabs(y)))/((fabs(x)+fabs(y))*(fabs(x)+fabs(y)));
}
double powermod (double x, double y){
    if (fabs(x)<TOL && fabs(y)<TOL) 
        return 0.0;
    int sx, sy; double power;
    sx=SIGN(x);
    sy=SIGN(y);
    power=fmin(fabs(x),fabs(y))*(x*x+y*y+2*fmax(fabs(x),fabs(y))*fmax(fabs(x),fabs(y)))/((fabs(x)+fabs(y))*(fabs(x)+fabs(y)));
    return (sx+sy)*0.5*power;
}

double minmod(double x, double y){
    return (SIGN(x)+SIGN(y))*0.5*fmin(fabs(x),fabs(y));
}

/*MUSCL*/
double parabola(double *x, double *y, int j, double eval, int N){
    double mm,mp,m;
    if (j>0){
        mm=(y[j]-y[j-1])/(x[j]-x[j-1]);
        if(j<N-1){
            mp=(y[j+1]-y[j])/(x[j+1]-x[j]);
        }
        else mp=0;
    }
    else{
        mm=0;
        mp=(y[j+1]-y[j])/(x[j+1]-x[j]);
    }
    m=minmod(mm,mp);
    return y[j]+m*(eval-x[j]);
}

double recta(double *x, double *y, int j, double eval, int N){
    double a, b;
    if(j==N-1){
      b=(y[j-1]-y[j])/(x[j-1]-x[j]);  
    }
    else 
        b=(y[j]-y[j+1])/(x[j]-x[j+1]);
    a=y[j]-b*x[j];
    return a+b*eval;
}
 
double cap(double a){
    if(a<0) return 0.01;
    else return a;
}
