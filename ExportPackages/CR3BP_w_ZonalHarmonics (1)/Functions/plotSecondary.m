function plotSecondary(secondary, alpha)
%%% Description
%       Plots secondary body in normalized, 3D coordinates
%       
% ------------------------------------------------------------------------
%%% Inputs
%       secondary - struct with fields .MR (mass ratio), .R_n (normalized
%                   radius), and .img (image file)
%       alpha     - alpha value for body (1 = solid, 0 = invisible)
% ------------------------------------------------------------------------
%%% Outputs
%
% ------------------------------------------------------------------------
% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Plot secondary
if nargin == 1
    if isfield(secondary,'img')
        plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img);
    elseif isfield(secondary,'color')
        plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color,1)
    else
        plotBody3(secondary.R_n,[1-secondary.MR,0,0],[0.196078431372549,0.392156862745098,0.784313725490196],1)
    end
elseif nargin == 2
    if isfield(secondary,'img')
        plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img, alpha)
    elseif isfield(secondary,'color')
        plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color,alpha)
    else
        plotBody3(secondary.R_n,[1-secondary.MR,0,0],[0.196078431372549,0.392156862745098,0.784313725490196],alpha)
    end
end

end % function