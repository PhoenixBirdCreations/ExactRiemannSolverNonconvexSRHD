#include "EOS_GGL.h"


double gmma(double r){
    return eos.g0+(eos.g1-eos.g0)*exp(-(r-eos.r0)*(r-eos.r0)/eos.s0/eos.s0);
}

double dgmma(double r){
    return -2.0*(r-eos.r0)/eos.s0/eos.s0*(eos.g1-eos.g0)*exp(-(r-eos.r0)*(r-eos.r0)/eos.s0/eos.s0);
}

double ddgmma(double rho){
  return (2*(rho-eos.r0)*(rho-eos.r0)/(eos.s0*eos.s0)-1)*(2*(eos.g1-eos.g0))/(eos.s0*eos.s0)*exp(-(rho-eos.r0)*(rho-eos.r0)/(eos.s0*eos.s0));
}

double dddgmma(double r){
  double s0=eos.s0, r0=eos.r0;
  return (12*(r-r0)-8*(r-r0)*(r-r0)*(r-r0)/s0/s0)*(eos.g1-eos.g0)*exp(-(r-r0)*(r-r0)/(s0*s0))/(s0*s0*s0*s0);
}

double nonlfactor( double css, double r){
    double ga=gmma(r);
    double dga=dgmma(r);
    double G=0.5*(1+ga+r*(2*ga*dga+r*ddgmma(r))/(ga*(ga-1)+r*dga));
    return sqrt(css)/r*(G-1.5*css);
}

double fundG(double r){
    double ga=gmma(r);
    double dga=dgmma(r);
    return 0.5*(1+ga+r*(2*ga*dga+r*ddgmma(r))/(ga*(ga-1)+r*dga));
}




