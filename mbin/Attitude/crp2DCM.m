function [ DCM ] = crp2DCM( q )
%%% Maps CRPs to DCM
%%% Inputs:
%       1) q - [3x1] CRP vector
%%% Outputs: 
%       1) DCM - [3x3] rotation matrix corresponding to CRP input
%=========================================================================
if size(q) == [1,3]
    q = q';
    
elseif size(q) ~= [3,1]
    warning('CRP input error')
end

%%% Direction Cosine Matrix
DCM = [1+q(1)^2-q(2)^2-q(3)^2, 2*(q(1)*q(2)+q(3)), 2*(q(1)*q(3)-q(2));...
       2*(q(2)*q(1)-q(3)), 1-q(1)^2+q(2)^2-q(3)^2, 2*(q(2)*q(3)+q(1));...
       2*(q(3)*q(1)+q(2)), 2*(q(3)*q(2)-q(1)), 1-q(1)^2-q(2)^2+q(3)^2];
DCM = DCM / (1 + q'*q);    

end

