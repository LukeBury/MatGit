function [tap]=kepdt(coep,JD,mu)
%Determines true anomaly at specified julian date from classical orbital 
%elements and reference date for those particular elements

% INPUTS
%coep: classical orbital elements of parking orbit
%JD: julian date
%mu: gravitational parameter of central body
%tol: tolerance of convergence

% OUTPUTS
% tap: true anomaly (deg)

ap = coep(1);
ep = coep(2);
ip = coep(3); %(deg)
op = coep(4); %(deg)
wp = coep(5); %(deg)
tap = coep(6)*pi/180; %(rad)
JD0=coep(7); %(days)


p=ap*(1-ep^2);
dt=(JD-JD0)*86400; %(sec)
r=p/(1+ep*cos(tap));
if ep~=0
   anomaly=ta2a(ep,tap); %(rad)
else
   anomaly=tap;  %(rad)
end

if ep<1.0
   n=sqrt(mu/ap^3); %mean motion
   M0=anomaly-ep*sin(anomaly); %(rad)
   M=M0+n*dt;
   anomaly=KepEqn(mu,M,ep,p,dt); 
elseif ep==1.0
   M0=anomaly+anomaly^3/3;
   anomaly=KepEqn(mu,M0,ep,p,dt); %(rad)
elseif ep>1.0
   n=sqrt(-mu/ap^3); %mean motion
   M0=ep*sinh(anomaly)-anomaly; %(rad)
   M=M0+n*dt;
   anomaly=KepEqn(mu,M,ep,p,dt); %(rad)
end
if ep~=0
   tap = a2ta(anomaly,ep,p,r)*180/pi; %(deg)
else
   tap = anomaly*180/pi; %(deg)
end
end
