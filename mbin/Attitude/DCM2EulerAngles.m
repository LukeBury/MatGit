function [ EA ] = DCM2Quat( C , Seq)
%%% Maps MRPs to DCM
%%% Inputs:
%       1) C - [3x3] rotation matrix
%       2) Seq - [str] rotation sequence
%%% Outputs: 
%       1) EA = [3x1] Euler Angles
%=========================================================================

if strcmp(Seq,'321') == 1
    ea1 = atan(C(1,2)/C(1,1));
    ea2 = -asin(C(1,3));
    ea3 = atan(C(2,3)/C(3,3));
end

EA = [ea1; ea2; ea3];

end

