function [ta] = a2ta(anomaly,e,p,r)
%Determines true anomaly from respective anomaly (eccentric, parabolic,
%hyperbolic).

%anomaly: eccentric, parabolic, hyperbolic (rad)
%e: eccentricity
%p: semiparameter [km or ER]
%r: radius (magnitude) [km or ER]

%ta: true anomaly (rad)

if e < 1.0 %elliptical orbit
    E=anomaly;
%  ta=atan2((sin(E)*sqrt(1-e^2))/(1-e*cos(E)),(cos(E)-e)/(1-e*cos(E)))*180/pi;   
   ta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
elseif e >1.0 %hyperbolic orbit
   H=anomaly;
   ta=atan2((-sinh(H)*sqrt(e^2-1))/(1-e*cosh(H)),(cosh(H)-e)/(1-e*cosh(H)));
elseif e==1.0 %parabolic orbit
    B=anomaly; 
    ta=atan2(p*B,(p-r));
end
end
