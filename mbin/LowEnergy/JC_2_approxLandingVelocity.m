function [landingVelocity_mps] = JC_2_approxLandingVelocity(JC_sc, prms, vNorm)
%%% Description
%       Based on the Jacobi Constant of the spacecraft, returns a
%       landing-velocity if the spacecraft were to land at [R, 0, 0] on the
%       secondary body
%       
% ------------------------------------------------------------------------
%%% Inputs
%       JC_sc - [1x1] Jacobi constant of spacecraft
%       prms - [struct] fields: u (mass ratio), R2 (radius of secondary)
%       vNorm  - [1x1] Velocity normalizing factor (km/s)
% ------------------------------------------------------------------------
%%% Outputs
%       landingVelocity_mps - [1x1] Landing velocity on postive-x surface
%                             (mps)
% ------------------------------------------------------------------------
% Created: 10/16/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Calculate Jacobi constant of positive-x surface
x_posXSurface = 1 - prms.u + prms.R2;
[JC_posXSurface] = getJacobiConstant_ZH([x_posXSurface, 0, 0, 0, 0, 0], prms);

%%% Calculate difference between sc JC and JC_L2
dJC = JC_sc - JC_posXSurface;

%%% Make sure the surface is accessible
if dJC > 0
    warning('Surface not reachable')
    landingVelocity_mps = NaN;
    return
end

%%% Convert the difference in JC to difference in velocity
landingVelocity_mps = sqrt(abs(dJC)) * vNorm * 1000;
end