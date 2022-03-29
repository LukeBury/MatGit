function [L2FlyoverVelocity_mps] = JC_2_L2FlyoverVelocity(JC_sc, prms, L2_BCR, vNorm)
%%% Description
%       Calculates a Jacobi Constant from an L2 flyover velocity       
%       (0 => JC_L2)
%       
% --------------------------------------------------------------
%%% Inputs
%       JC_2_L2FlyoverVelocity - [1x1] Jacobi constant to conver to L2
%           Flyover Velocity
%       prms     - [struct] (u)Mass ratio of CR3BP system, (n) normalized
%                            mean motion
%       L2_BCR - [1x3] Position coordinates of L2 in normalized CR3BP
%       vNorm  - [1x1] Velocity normalizing factor (km/s)
% --------------------------------------------------------------
%%% Outputs
%       L2FlyoverVelocity_mps - [1x1] L2 flyover velocity (m/s)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Calculate Jacobi constant of L2
% [JC_L2] = JacobiConstantCalculator(prms, L2_BCR, [0,0,0]);
JC_L2 = getJacobiConstant_ZH([L2_BCR,0,0,0], prms);

%%% Calculate difference between sc JC and JC_L2
dJC_L2 = JC_L2 - JC_sc;

%%% Convert the difference in JC to difference in velocity
L2FlyoverVelocity_kps = sqrt(dJC_L2)*vNorm;
L2FlyoverVelocity_mps = L2FlyoverVelocity_kps * 1000;

end