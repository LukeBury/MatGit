function [rf, vf] = fandgdta(r0,v0,dta,mu)
%Functions determines the state [rf vf] propagated by a change in 
%true anomaly from an initial state [r0 v0]

%INPUTS
% r0: initial radius  (3x1)
% v0: initial velocity (3x1)
% dta: change in true anomaly (radians)
% mu: gravitational parameter

%OUTPUTS
% rf: final radius   (3x1)
% vf: final velocity (3x1)

%First Determine f and g functions
h=[r0(2)*v0(3) - r0(3)*v0(2);...%cross product 
   r0(3)*v0(1) - r0(1)*v0(3);...
   r0(1)*v0(2) - r0(2)*v0(1)];
h2=(h'*h);
h=sqrt(h2);
r0mag=sqrt(r0'*r0);
vr0=(v0'*r0)/r0mag;
sdta=sin(dta);
cdta=cos(dta);

r=h2/mu/(1+(h2/(mu*r0mag)-1)*cdta-h*vr0*sdta/mu);
f=1-mu*r*(1-cdta)/h2;
g=r*r0mag*sdta/h;

%Determine fdot and gdot functions
fdot=mu/h*(vr0/h*(1-cdta)-sdta/r0mag);
gdot=1-mu*r0mag/h2*(1-cdta);
rf=f*r0+g*v0;
vf=fdot*r0+gdot*v0;
end