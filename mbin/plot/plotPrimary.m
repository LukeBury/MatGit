function plotPrimary(primary,secondary)
%%% Description
%       Plots primary body in normalized, 3D coordinates
%       
% ------------------------------------------------------------------------
%%% Inputs
%       primary   - struct with fields .R (radius in km) and .img (image 
%                   file)
%       secondary - struct with fields .MR (mass ratio) and .a (semi-major
%                   axis in km)
% ------------------------------------------------------------------------
%%% Outputs
%
% ------------------------------------------------------------------------
% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km

%%% Plot Primary
plotBodyTexture3(primary.R/rNorm,[-secondary.MR,0,0],primary.img);

end % function