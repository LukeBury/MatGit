function [v_BCR] = v_BodyCI2BodyCR(v_BCI, r_BCI, times, w )
%%% Description
%       Convert velocities in inertial frame to velocities in rotating
%       frame with a common origin
%       
% ------------------------------------------------------------------------
%%% Inputs
%       v_BCI - [nx3] matrix of velocity values in Body-Centered-Inertial
%               frame
%       r_BCI - [nx3] matrix of position values in Body-Centered-Inertial
%               frame
%       times - [nx1] vector of times associated with positions (sec)
%       w     - angular velocity of the system (rad/s) 
% ------------------------------------------------------------------------
%%% Outputs
%       v_BCR - [nx3] matrix of velocity values in Body-Centered-Rotating
%               frame
% ------------------------------------------------------------------------
% Created: 07/16/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================

%%% Dimension checks
if size(v_BCI,2) ~= 3
    warning('Wrong dimensions for input matrix')
    return
elseif size(v_BCI,1) ~= length(times)
    warning('Dimensions of position matrix and times don''t agree')
    return
elseif size(v_BCI) ~= size(r_BCI)
    warning('Dimensions of v_BCR and r_BCI don''t agree')
    return
end


%%% Preallocate output vector
v_BCR = zeros(size(v_BCI));

%%% Rotate BCR into BCI
for kk = 1:length(times)
    % Create rotation angle (rad)
    theta = times(kk)*w;
    
    % Transport Theorem
    v_BCR(kk,:) = R3(v_BCI(kk,:) - cross([0,0,w],r_BCI(kk,:)), -theta);
end

end % function