function [dva]= hyp2parkoptOPT(vinfamag,coepa,mua)
% warning off
% Determines the hyperbolic arrival trajectory to a parking orbit of
% orbital elements a & e assuming i >= declination of vinf
% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory
%------INPUTS-----

% vinfamag: vinf magnitude
% coepa(6): classical orbital elements of parking departure orbit 
% mua:  gravitational parameter of departure planet
% rsoi: radius of sphere of influence of planet
 
%------OUTPUTS-----
% dva: Delta-V required for hyperbolic depature trajectory planet frame

ap = coepa(1);
ep = coepa(2);
% ip = coepa(3);
% op = coepa(4);
% wp = coepa(5);
% tap = 0; %true anomaly; set to zero b/c arrival at periapsis
% rppamag = ap*(1-ep);%periapsis radius (km)
% rapamag = ap*(1+ep);%apoapsis radius (km)
vcmag=(mua/(ap*(1-ep))); %velocity^2 of circular orbit with same periapsis 
% dva=sqrt(2*vcmag+vinfamag*vinfamag)-sqrt(2*vcmag*rapamag/(rppamag+rapamag)); %scalar delta-V 
dva=abs(sqrt(2*vcmag+vinfamag*vinfamag)-sqrt(vcmag*(1+ep))); %scalar delta-V 

end
