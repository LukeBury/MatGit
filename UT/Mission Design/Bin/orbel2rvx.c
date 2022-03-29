/* function [r, v] = orbel2rvx(mu,a,e,i,o,w,TA) */

#include <math.h>
#include "mex.h"



void orbel2rvx(double mu,double a, double e,
        double i,double o,double w,double TA,
        double *r, double *v) 
{
#define pi 3.14159265358979323846264338327    
double deg2rad,rmag,vmag,FPA,taw;
double coso,sino,cosi,sini,costaw,sintaw;

deg2rad=pi/180;
i=i*deg2rad;
o=o*deg2rad;
w=w*deg2rad;
TA=TA*deg2rad;

rmag = a*(1-e*e)/(1+e*cos(TA));
vmag = sqrt(mu*(2/rmag-1/a));
FPA = atan2(e*sin(TA),1+e*cos(TA));
taw=TA+w;
coso=cos(o);
sino=sin(o);
cosi=cos(i);
sini=sin(i);
costaw=cos(taw);
sintaw=sin(taw);
    r[0] = rmag*(costaw*coso - sintaw*cosi*sino);
    r[1] = rmag*(costaw*sino + sintaw*cosi*coso);
    r[2] = rmag*(sintaw*sini);
    taw=taw-FPA;
    costaw=cos(taw);
    sintaw=sin(taw);
    v[0] = vmag*(-sintaw*coso - costaw*cosi*sino);
    v[1] = vmag*(-sintaw*sino + costaw*cosi*coso);
    v[2] = vmag*(costaw*sini);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    double *mu,*a,*e,*i,*o,*w,*TA;
    double *r,*v;
    mu=mxGetPr(prhs[0]);
    a=mxGetPr(prhs[1]);
    e=mxGetPr(prhs[2]);
    i=mxGetPr(prhs[3]);
    o=mxGetPr(prhs[4]);
    w=mxGetPr(prhs[5]);
    TA=mxGetPr(prhs[6]);
    plhs[0]=mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(3,1,mxREAL);
    r=mxGetPr(plhs[0]);
    v=mxGetPr(plhs[1]);  
    orbel2rvx(*mu,*a,*e,*i,*o,*w,*TA,r,v);

}