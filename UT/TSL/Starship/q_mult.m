function r = q_mult(p,q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION: q_mult
%BY: Sebastian Munoz
%CREATED: Feb 10, 2009
%
%This function returns the product of two quaternions p and q 
%
%
%INPUTS:
%     p ;  4x1 quaternion p
%     q ;  4x1 quaternion r
%
%OUTPUTS:
%     r ;  4x1 output quaternion r
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
%   02/10/2009          SMT          function created
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = zeros(4,1);

p = [1;0;1;0];
q = [0;2;0;1];

if (length(p) == 4 && length(q) == 4)
    p1 = p(1,1);
    p2 = p(2,1);
    p3 = p(3,1);
    p4 = p(4,1);
    
    q1 = q(1,1);
    q2 = q(2,1);
    q3 = q(3,1);
    q4 = q(4,1);
    
    r(1,1) = p4*q1 + p1*q4 + p2*q3 - p3*q2;
    r(2,1) = p4*q2 - p1*q3 + p2*q4 + p3*q1;
    r(3,1) = p4*q3 + p1*q2 - p2*q1 + p3*q4;
    r(4,1) = p4*q4 - p1*q1 - p2*q2 - p3*q3;
    
else
    error('Quaternion is not the right length');
end