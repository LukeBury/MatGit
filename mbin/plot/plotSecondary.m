function plotSecondary(secondary)
%%% Description
%       Plots secondary body in normalized, 3D coordinates
%       
% ------------------------------------------------------------------------
%%% Inputs
%       secondary - struct with fields .MR (mass ratio), .R_n (normalized
%                   radius), and .img (image file)
% ------------------------------------------------------------------------
%%% Outputs
%
% ------------------------------------------------------------------------
% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Plot secondary
plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img);

end % function