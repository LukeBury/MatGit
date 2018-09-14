function [rP,vP]=ephstate(type,event,planet,JD,coep,rnb)
% type: Ephemeris type (1: DE405 ephemeris, 2: 2nd Order Ephemeris
% event: Defines if departure node is a planet or not
% planet: Departure planet number, found in "pheph.M"
% JD: Departure julian date
% coep: Classical Orbital Elements for non-planetary body
% rnb: non-body position vector (typically non-zero if used)

muj = 1.26712764e+008;      %km3/s2 Gravitational Parameter of Sun
dxtol=1E-7;
Rj=71492; 

if planet==0
    if (event==40)||(event==41)
    rP(1:3,1)=rnb.*Rj;
    vP(1:3,1)=0;
    elseif (event==17)||(event==27)
%     coep(6)=kepdt(coep(1:7),JD,mu);   
    [rP,vP]=orbel2rv(muj,coep(1),coep(2),coep(3),coep(4),coep(5),coep(6));
    [rP,vP]=kepuv(muj,(JD-coep(7))*86400,rP,vP,0,dxtol);
    end
else
    
    [rP,vP]=pleph2o(planet,JD);
     	
end
end