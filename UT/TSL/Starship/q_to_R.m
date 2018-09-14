function R = q_to_R(q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION: q_to_R
%BY: Sebastian Munoz 
%CREATED: Feb 11, 2009
%
%This function transforms a unit quaternion into a direction cosine matrix.
%
%
%INPUTS:
%     q ; 4x1 unit quaternion
%     
%OUTPUTS:
%     R ;  3x3 Direction Cosine Matrix 
%
%EXAMPLE FUNCTION CALL:
%   N/A
%
%REFERENCES:
% Kuipers, J.B., "Quaternions and Rotation Sequences: A Primer with
% Applications to Orbits, Aerospace, and Virtual Reality", Princeton
% University Press, 1999. 
%   
%
%
%REVISION HISTORY:
%      date         modified by      comments on revision
%  ------------    -------------    --------------------------------------
%   02/11/2009          SMT          function created
%   02/17/2009          JAC          fixed q_norm call for warning check
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = [1; 0; 1; 0];
eps = 1e-8;
R = zeros(3,3);
%if (length(q) == 4)
%    if (q_norm(q) > 1 + eps || q_norm(q) < 1 - eps)
%        warning('Quaternion is not of unity norm');
%    end
    q1 = q(1,1);
    q2 = q(2,1);
    q3 = q(3,1);
    q4 = q(4,1);
    
    R(1,1) = 2*q4^2 - 1 + 2*q1^2;
    R(1,2) = 2*q1*q2 + 2*q4*q3;
    R(1,3) = 2*q1*q3 - 2*q4*q2;
    R(2,1) = 2*q1*q2 - 2*q4*q3;
    R(2,2) = 2*q4^2 - 1 + 2*q2^2;
    R(2,3) = 2*q2*q3 + 2*q4*q1;
    R(3,1) = 2*q1*q3 + 2*q4*q2;
    R(3,2) = 2*q2*q3 - 2*q4*q1;
    R(3,3) = 2*q4^2 - 1 + 2*q3^2;
%else
%    error('Quaternion is not the right length');
2*q4^2 - 1 + 2*q1^2
end