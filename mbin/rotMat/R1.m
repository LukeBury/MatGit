function [m2] = R1(m1,w)
%%% ***************************************************
%%% NEED (-)THETA FOR EULER ANGLE ROTATIONS
%%% ***************************************************
%%% Given:
% 1) 3x1 or 1x3 matrix to be rotated on X axis
% 2) Rotation angle in radians
%%% Outputs:
% 1) Rotated vector with same dimensions as input
% =========================================================================
R = [1 0 0; 0 cos(w) -sin(w); 0 sin(w) cos(w)]; % Rot Matrix
if any(size(m1) == [3 1]) % If m1 is column vector
    m2 = R*m1; % Return rotated column vector

elseif any(size(m1) == [1 3]) % If m1 is row vector
    m1 = m1'; % Make m1 a column vector
    m2 = R*m1; % Rotated column vector
    m2 = m2'; % Return m2 as a row vector

else
    warning('Wrong dimensions fed to R1')
    return
end
end
