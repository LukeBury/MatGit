function q = R_to_q(R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION: q_to_R
%BY: Sebastian Munoz 
%CREATED: Feb 11, 2009
%
%This function transforms a direction cosine matrix into a unit quaternion.
%
%
%INPUTS:
%     R ;  3x3 Direction Cosine Matrix
%     
%OUTPUTS:
%     q ; 4x1 unit quaternion 
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
%   02/17/2009          JAC          correction to q4 calculation
%   02/18/2009          JAC          added successive rotation check to
%                                      avoid singularities when q4 is near 0
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=[1 0 2;0 -1 0; 2 0 1];
eps = 1e-8;
q = zeros(4,1);
[n,o] = size(R);
if (n == 3 && o == 3)
    detR = det(R);
    if (detR > 1 + eps || detR < 1 - eps)
        warning('Direction Cosine Matrix does not have a determinant equal to 1');
    end
    m11 = R(1,1);
    m12 = R(1,2);
    m13 = R(1,3);
    m21 = R(2,1);
    m22 = R(2,2);
    m23 = R(2,3);
    m31 = R(3,1);
    m32 = R(3,2);
    m33 = R(3,3);
    
    q4 = sqrt(m11 + m22 + m33 + 1)/2;
    
    RotIndex = 1;
    LoopCheck = q4;
    qfix = [0;0;0;1];
    
    while LoopCheck < 0.01  %Need to rotate solution to avoid singularity
        
       if RotIndex <= 3    %Loop through the three primary axes
           
          %Perform rotation about axis indicated by RotIndex
             qfix = zeros(4,1);
             qfix(RotIndex) = 1;
             Rfix = q_to_R(qfix);
             Rcurrent = Rfix*R;
             
          %Recompute matrix elements   
             m11 = Rcurrent(1,1);
             m12 = Rcurrent(1,2);
             m13 = Rcurrent(1,3);
             m21 = Rcurrent(2,1);
             m22 = Rcurrent(2,2);
             m23 = Rcurrent(2,3);
             m31 = Rcurrent(3,1);
             m32 = Rcurrent(3,2);
             m33 = Rcurrent(3,3);
             
          %Recompute q4
             q4 = sqrt(m11 + m22 + m33 + 1)/2;
             LoopCheck = q4;
             
       else
           
          LoopCheck = 1;
          
       end
       
       RotIndex = RotIndex + 1;
       
    end
    
    q(1,1) = (m23 - m32)/(4*q4);
    q(2,1) = (m31 - m13)/(4*q4);
    q(3,1) = (m12 - m21)/(4*q4);
    q(4,1) = q4;
    
   
    %q = q_mult(q,qfix);
    
else
    error('Direction Cosine Matrix is not the right size');
end