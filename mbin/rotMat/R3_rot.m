function [vector_out] = R3_rot(vector_in,theta_rad)
%%% Descrip: Writes body-vectors in inertial frame: v_N = R3 * v_B ...
%%% so [sqrt(2)/2, sqrt(2)/2, 0] = R3(pi/4) * [1, 0, 0]
%%% ***************************************************
%%% NEED -THETA FOR EULER ANGLE ROTATIONS
%%% ***************************************************
%%% Inputs:
% 1) 3x1 or 1x3 matrix to be rotated on Z axis
% 2) Rotation angle in radians
%%% Outputs:
% 1) Rotated vector with same dimensions as input
% =========================================================================
R = [cos(theta_rad) -sin(theta_rad) 0; sin(theta_rad) cos(theta_rad) 0; 0 0 1]; % Rot Matrix
if any(size(vector_in) == [3 1]) % If vector_in is column vector
    vector_out = R*vector_in; % Return rotated column vector

elseif any(size(vector_in) == [1 3]) % If vector_in is row vector
    vector_in = vector_in'; % Make vector_in a column vector
    vector_out = R*vector_in; % Rotated column vector
    vector_out = vector_out'; % Return vector_out as a row vector

else
    warning('Wrong dimensions fed to R3')
    return
end
end