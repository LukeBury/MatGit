function [ Quat ] = DCM2Quat( C )
%%% Maps MRPs to DCM
%%% Inputs:
%       1) DCM - [3x3] rotation matrix
%%% Outputs: 
%       1) Quat - [4x1] quaternion vector
%=========================================================================

B0 = sqrt(.25*(1 + trace(C)));
B1 = sqrt(.25*(1 + 2*C(1,1) - trace(C)));
B2 = sqrt(.25*(1 + 2*C(2,2) - trace(C)));
B3 = sqrt(.25*(1 + 2*C(3,3) - trace(C)));

Quat = [B0; B1; B2; B3];

end

