function [AB] = getDCM(A_N, B_N)
%%% Description
%       Given vectors describing two different sets of reference frames, will return a
%       DCM to rotate vectors between the two frames
%       
% --------------------------------------------------------------
%%% Inputs
%       A_N - [3x3] A-frame as described in the N-frame. Matrix of unit vectors,
%             such that: A_n = [Ax_N, Ay_N, Az_N], 
%             where, Ax_N is a [3x1] describing the position of the 'x'
%             or 'first' axis in the N frame, so 
%             Ax_N = [Axx_N, Axy_N, Axz_N]'.
%       B_N - [3x3] B-frame as described in the N-frame. Matrix of unit vectors,
%             such that: B_n = [Bx_N, By_N, Bz_N], 
%             where, Bx_N is a [3x1] describing the position of the 'x'
%             or 'first' axis in the N frame, so 
%             Bx_N = [Bxx_N, Bxy_N, Bxz_N]'.
%
% *** Note: The 'N' frame can just be the A frame (or B frame). In the case
% where the A frame and N frame are the same, A_N = eye(3);
% --------------------------------------------------------------
%%% Outputs
%       AB - [3x3] matrix that rotates a [3x1] vector from B frame to A
%       frame, such that vector_in_A = AB * vector_in_B
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Check orthonormal
% ------------------------------------
%%% Ensure orthogonal columns:
tol = 9e-16;
if (abs(dot(A_N(:,1),A_N(:,2))) > tol) || (abs(dot(A_N(:,1),A_N(:,3))) > tol) || (abs(dot(A_N(:,2),A_N(:,3))) > tol)
    warning('A not orthogonal')
    return
end

if (abs((dot(B_N(:,1),B_N(:,2)) > tol))) || (abs((dot(B_N(:,1),B_N(:,3)) > tol))) || (abs((dot(B_N(:,2),B_N(:,3)) > tol)))
    warning('B not orthogonal')
    return
end

%%% Enforce unit length
A_N(:,1) = A_N(:,1) ./ norm(A_N(:,1));
A_N(:,2) = A_N(:,2) ./ norm(A_N(:,2));
A_N(:,3) = A_N(:,3) ./ norm(A_N(:,3));

B_N(:,1) = B_N(:,1) ./ norm(B_N(:,1));
B_N(:,2) = B_N(:,2) ./ norm(B_N(:,2));
B_N(:,3) = B_N(:,3) ./ norm(B_N(:,3));

% ------------------------------------
%%% Create Rotation Matrix
% ------------------------------------
AB = A_N' * B_N;

end