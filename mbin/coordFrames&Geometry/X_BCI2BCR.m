function [X_BCR] = X_BCI2BCR(X_BCI, times, w )
%%% Description
%       Convert body-centered-inertial state to body-centered-rotating 
%       state
%       
% ------------------------------------------------------------------------
%%% Inputs
%       X_BCI - [nx3] matrix of state values in Body-Centered-Inertial
%               frame
%       times - [nx1] vector of times associated with positions (sec)
%       w     - angular velocity of the system (rad/s) 
% ------------------------------------------------------------------------
%%% Outputs
%       X_BCR - [nx3] matrix of state values in Body-Centered-Rotating
%               frame
% ------------------------------------------------------------------------
% Created: 07/17/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Call other functions to create r_BCR and v_BCR, then combine
% -------------------------------------------------
%%% Seperate input data
r_BCI = X_BCI(:,1:3);
v_BCI = X_BCI(:,4:6);

%%% Call to functions
[ r_BCR ] = r_BCI2BCR( r_BCI, times, w );
[v_BCR] = v_BodyCI2BodyCR(v_BCI, r_BCI, times, w );

%%% Combine results
X_BCR = [r_BCR, v_BCR];

end % function