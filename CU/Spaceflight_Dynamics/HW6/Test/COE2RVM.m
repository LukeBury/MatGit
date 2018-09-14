function [Rijk, Vijk] = COE2RVM(a,e,i,O,w,M,mu)

E=Mean2E(M,e);
nu=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
if nu<0
    nu=nu+2*pi;
end

R3O=[cos(-O),sin(-O),0;-sin(-O),cos(-O),0;0,0,1];
R1i=[1,0,0;0,cos(-i),sin(-i);0,-sin(-i),cos(-i)];
R3w=[cos(-w),sin(-w),0;-sin(-w),cos(-w),0;0,0,1];
ROI=R3O*R1i*R3w;

p = a*(1 - e^2);
Rpqw=[p*cos(nu)/(1+e*cos(nu));p*sin(nu)/(1+e*cos(nu));0];
Vpqw=[-sqrt(mu/p)*sin(nu);sqrt(mu/p)*(e+cos(nu));0];
Rijk=ROI*Rpqw;
Vijk=ROI*Vpqw;