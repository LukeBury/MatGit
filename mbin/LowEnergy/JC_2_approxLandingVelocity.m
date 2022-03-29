function [landingVelocity_mps] = JC_2_approxLandingVelocity(JC_sc, prms, vNorm)
%%% Description
%       Based on the Jacobi Constant of the spacecraft, returns a
%       landing-velocity if the spacecraft were to land at [R, 0, 0] on the
%       secondary body
%       
% ------------------------------------------------------------------------
%%% Inputs
%       JC_sc  - [nx1] Jacobi constant of spacecraft
%       prms   - [struct] fields: u (mass ratio), R2 (radius of secondary),
%                               n (normalized mean motion)
%       vNorm  - [1x1] Velocity normalizing factor (km/s)
% ------------------------------------------------------------------------
%%% Outputs
%       landingVelocity_mps - [nx1] Landing velocity on postive-x surface
%                             (mps)
% ------------------------------------------------------------------------
% Created: 10/16/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Calculate Jacobi constant of positive-x surface
x_posXSurface = 1 - prms.u + prms.R2;
[JC_posXSurface] = getJacobiConstant_ZH([x_posXSurface, 0, 0, 0, 0, 0], prms);

landingVelocity_mps = zeros(size(JC_sc));
for JC_i = 1:length(JC_sc)
    %%% Calculate difference between sc JC and JC_L2
    dJC = JC_sc(JC_i) - JC_posXSurface;

    %%% Make sure the surface is accessible
    if dJC > 0
        warning('Surface not reachable')
        landingVelocity_mps(JC_i) = NaN;
        return
    end

    %%% Convert the difference in JC to difference in velocity
    landingVelocity_mps(JC_i) = sqrt(abs(dJC)) * vNorm * 1000;
end
end