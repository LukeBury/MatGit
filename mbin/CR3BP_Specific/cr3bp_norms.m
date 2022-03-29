function [rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, G)
%%% Description
% Determine the normalizing constants for classical CR3BP system
%       
% ------------------------------------------------------------------------
%%% Inputs
% primary - [struct] struct of primary body parameters
% secondary - [struct] struct of secondary body parameters
% G - [1x1] Gravitational constant (km^3 * kg^-1 * s^-2)
% ------------------------------------------------------------------------
%%% Outputs
% rNorm - [1x1] distance-normalizing constant (km)
% tNorm - [1x1] time-normalizing constant (sec)
% vNorm - [1x1] velocity-normalizing constant (km/sec)
% ------------------------------------------------------------------------
% Created: 11/05/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Determine theoretical time period of circular-restricted system
Tp = 2*pi*sqrt((secondary.a^3) / (G*(secondary.mass + primary.mass))); % sec

%%% Set distance-normalizing constant
rNorm = secondary.a;         % n <-> km

%%% Set time-normalizing constant
tNorm = Tp / (2*pi); % n <-> sec

%%% Set velocity-normalizing constant
vNorm = rNorm / tNorm;

end