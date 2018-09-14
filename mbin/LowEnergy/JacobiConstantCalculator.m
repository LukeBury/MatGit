function [JCs] = JacobiConstantCalculator(u,rBCR_n,vBCR_n)
%%% Inputs:
% 1) mass ratio of system [1x1]
% 2) Positions in normalized BCR frame ............[nx3]
% 3) Velocities in normalized BCR frame ...........[nx3]
%%% Outputs:
% 1) Jacobi constants ............................[nx1]
% ========================================================================
if size(rBCR_n,1) ~= size(vBCR_n,1)
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
    JCs(k) = x^2 + y^2 + 2*(1-u)/r1 + 2*u/r2 - dx^2 - dy^2 - dz^2;
end
end