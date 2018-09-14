function [anomaly] = ta2a(e,ta)
%Determines anomaly (eccentric, parabolic,hyperbolic)from respective true 
%anomaly .
%anomaly: eccentric, parabolic, or hyperbolic(rad)
%e: eccentricity
%ta: true anomaly (rad)

if e < 1.0 %for elliptical orbit
%   anomaly=atan2((sin(ta)*sqrt(1-e^2))/(1+e*cos(ta)),(e+cos(ta))/(1+e*cos(ta)));
    anomaly=2*atan(sqrt((1-e)/(1+e))*tan(ta/2));
elseif e >1.0 %for hyperbolic orbit
    anomaly=asinh(sin(ta)*sqrt(e^2-1)/(1+e*cos(ta)));
elseif e==1.0 %for parabolic orbit
    anomaly=atan(ta/2)*pi/180;
end
end