function [dvd]= park2hypoptOPT(vinfdmag,coepd,mud)
% warning off
% Determines the hyperbolic departure trajectory from a parking orbit of
% orbital elements a, e, and i

% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory 

%------INPUTS-----

% vinfdmag: vinf magnitude
% coepd(6): classical orbital elements of parking departure orbit 
% mud:  gravitational parameter of departure planet
 
%------OUTPUTS-----
% dvd: Delta-V required for hyperbolic depature trajectory planet frame

ap = coepd(1);
ep = coepd(2);
% ip = coepd(3);
% op = coepd(4);
% wp = coepd(5);
% tap = 0; %true anomaly; set to zero b/c departure at periapsis

% rppdmag = ap*(1-ep); %periapsis radius (km)
% rapdmag = ap*(1+ep); %apoapsis radius (km)
vcmag=(mud/(ap*(1-ep))); %velocity^2 of circular orbit with same periapsis 
% dvd=sqrt(2*vcmag+vinfdmag*vinfdmag)-sqrt(2*vcmag*rapdmag/(rppdmag+rapdmag)); %scalar delta-V 
dvd=abs(sqrt(2*vcmag+vinfdmag*vinfdmag)-sqrt(vcmag*(1+ep))); %scalar delta-V 

end