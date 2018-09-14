function [ DCM ] = quat2DCM( B )
%%% Maps quaternion (Euler Parameters) to DCM
%%% Inputs:
%       1) B - [4x1] unit quaternion
%%% Outputs: 
%       1) DCM - [3x3] rotation matrix corresponding to quaternion input
%=========================================================================
%%% Ensuring unit length
B = B./norm(B);

%%% Direction Cosine Matrix
DCM = [B(1)^2+B(2)^2-B(3)^2-B(4)^2, 2*(B(2)*B(3)+B(1)*B(4)), 2*(B(2)*B(4)-B(1)*B(3));...
       2*(B(2)*B(3)-B(1)*B(4)), B(1)^2-B(2)^2+B(3)^2-B(4)^2, 2*(B(3)*B(4)+B(1)*B(2));...
       2*(B(2)*B(4)+B(1)*B(3)), 2*(B(3)*B(4)-B(1)*B(2)), B(1)^2-B(2)^2-B(3)^2+B(4)^2];

end

