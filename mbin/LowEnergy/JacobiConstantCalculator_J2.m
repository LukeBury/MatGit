function [JCs] = JacobiConstantCalculator_J2(u,rBCR_n,vBCR_n, R1_n, R2_n, J21, J22)
%%% Inputs:
% 1) Mass ratio of system [1x1]
% 2) Positions in normalized BCR frame ............[nx3]
% 3) Velocities in normalized BCR frame ...........[nx3]
% 4) Normalized radius of primary 
% 5) Normalized radius of secondary
% 6) J2 value of primary
% 7) J2 value of secondary
%%% Outputs:
% 1) Jacobi constants ............................[nx1]
% ========================================================================
if size(rBCR_n) ~= size(vBCR_n)
    warning('Mismatch of input dimensions')
    return
end

rB1 = [-u, 0, 0];
rB2 = [1-u, 0, 0];

JCs = zeros(size(rBCR_n,1),1);
for k = 1:size(rBCR_n,1)
    x = rBCR_n(k,1);
    y = rBCR_n(k,2);
    z = rBCR_n(k,3);
    dx = vBCR_n(k,1);
    dy = vBCR_n(k,2);
    dz = vBCR_n(k,3);
    r1 = norm(rB1 - [x,y,z]);
    r2 = norm(rB2 - [x,y,z]);
    JCs(k) = x^2 + y^2 + 2*(1-u)/r1 + 2*u/r2 - dx^2 - dy^2 - dz^2 + R1_n^2*J21*(1-u)*(r1^2-3*z^2)/(r1^5) + R2_n^2*J22*u*(r2^2-3*z^2)/(r2^5);
end


% rB1 = [-u, 0, 0];
% rB2 = [1-u, 0, 0];
% JCs = zeros(size(rBCR_n,1),1);
% for k = 1:size(rBCR_n,1)
%     r1 = norm(rB1 - rBCR_n(k,:));
%     r2 = norm(rB2 - rBCR_n(k,:));
%     JCs(k) = rBCR_n(k,1)^2 + rBCR_n(k,2)^2 + 2*(1-u)/r1 + 2*u/r2 - vBCR_n(k,1)^2 - vBCR_n(k,2)^2 - vBCR_n(k,3)^2 + R1_n^2*J21*(1-u)*(r1^2-3*rBCR_n(k,3)^2)/(r1^5) + R2_n^2*J22*u*(r2^2-3*rBCR_n(k,3)^2)/(r2^5);
% end
end