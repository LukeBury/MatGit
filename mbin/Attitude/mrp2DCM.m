function [ DCM ] = mrp2DCM( s )
%%% Maps MRPs to DCM
%%% Inputs:
%       1) s - [3x1] MRP vector
%%% Outputs: 
%       1) DCM - [3x3] rotation matrix corresponding to CRP input
%=========================================================================
if size(s) == [1,3]
    s = s';
    
elseif size(s) ~= [3,1]
    warning('MRP input error')
end

sSq = norm(s)^2;

%%% Direction Cosine Matrix
DCM = zeros(3,3);
DCM(1,1) = 4*(s(1)^2-s(2)^2-s(3)^2) + (1-sSq)^2;
DCM(1,2) = 8*s(1)*s(2) + 4*s(3)*(1-sSq);
DCM(1,3) = 8*s(1)*s(3) - 4*s(2)*(1-sSq);
DCM(2,1) = 8*s(2)*s(1) - 4*s(3)*(1-sSq);
DCM(2,2) = 4*(-s(1)^2+s(2)^2-s(3)^2) + (1-sSq)^2;
DCM(2,3) = 8*s(2)*s(3) + 4*s(1)*(1-sSq);
DCM(3,1) = 8*s(3)*s(1) + 4*s(2)*(1-sSq);
DCM(3,2) = 8*s(3)*s(2) - 4*s(1)*(1-sSq);
DCM(3,3) = 4*(-s(1)^2-s(2)^2+s(3)^2) + (1-sSq)^2;
  
DCM = DCM./(1+sSq)^2;
end

