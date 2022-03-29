function [X_BCI] = X_BCR2BCI(X_BCR, times, w )
%%% Description
%       Convert body-centered-inertial state to body-centered-rotating 
%       state
%       
% ------------------------------------------------------------------------
%%% Inputs
%       X_BCR - [nx3] matrix of state values in Body-Centered-Rotating
%               frame
%       times - [nx1] vector of times associated with positions (sec)
%       w     - angular velocity of the system (rad/s) 
% ------------------------------------------------------------------------
%%% Outputs
%       X_BCI - [nx3] matrix of state values in Body-Centered-Inertial
%               frame
% ------------------------------------------------------------------------
% Created: 07/17/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Call other functions to create r_BCI and v_BCI, then combine
% -------------------------------------------------
%%% Seperate input data
r_BCR = X_BCR(:,1:3);
v_BCR = X_BCR(:,4:6);

%%% Call to functions
[ r_BCI ] = r_BCR2BCI( r_BCR, times, w );
[v_BCI] = v_BodyCR2BodyCI(v_BCR, r_BCI, times, w );

%%% Combine results
X_BCI = [r_BCI, v_BCI];
end % function
