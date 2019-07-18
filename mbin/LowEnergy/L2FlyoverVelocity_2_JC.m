function [JC_sc] = L2FlyoverVelocity_2_JC(L2FlyoverVelocity_mps, mu, L2_BCR, vNorm)
%%% Description
%       Calculates a Jacobi Constant from an L2 flyover velocity       
%       (0 => JC_L2)
%       
% --------------------------------------------------------------
%%% Inputs
%       L2FlyoverVelocity_mps - [1x1] Desired L2 flyover velocity (m/s) 
%       mu - [1x1] Mass ratio of CR3BP system 
%       L2_BCR - [1x3] Position coordinates of L2 in normalized CR3BP
%       vNorm - [1x1] Velocity normalizing factor (km/s)
% --------------------------------------------------------------
%%% Outputs
%       JC_sc - [1x1] Jacobi constant corresponding to L2 flyover velocity
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Convert L2 flyover velocity to km/s
L2FlyoverVelocity_kps = L2FlyoverVelocity_mps/1000;

%%% Convert the L2 flyover velocity to a component of Jacobi constant
dJC_L2 = (L2FlyoverVelocity_kps/vNorm)^2;

%%% Calculate Jacobi constant of L2
[JC_L2] = JacobiConstantCalculator(mu, L2_BCR, [0,0,0]);

%%% Calculate Jacobi constant corresponding to L2 flyover velocity
JC_sc = JC_L2-dJC_L2;

end