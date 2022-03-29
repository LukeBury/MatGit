function [X_SCI] = X_BaCR2SCI(X_BaCR, times, prms)
%%% Description
% Converting a state from bary-centered rotating coordinates 
% (standard CR3BP) to secondary-centered inertial
%       
% ------------------------------------------------------------------------
%%% Inputs
% X_BaCR - [nx6] Array of state values in bary-centered-rotating frame
% times  - [nx1] array of times associated with the states
% prms   - [structure] Required fields: u, n
% ------------------------------------------------------------------------
%%% Outputs
% X_SCI - [nx6] Array of state values in secondary-centered inertial frame
% ------------------------------------------------------------------------
% Created: 8/18/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Convert BaCR states to BaCI (rotating to inertial)
X_BaCI = X_BaCR2BaCI(X_BaCR, times, prms);

%%% Build state vector history for secondary body in BaCR (note: its state is
%%% constant here)
X_secondary_BaCR = repmat([1-prms.u, 0, 0, 0, 0, 0],length(times),1);

%%% Get the inertial state vector for the secondary body in bary-cenetered coordinates
X_secondary_BaCI = X_BaCR2BaCI(X_secondary_BaCR,times, prms);

%%% Calculate inertial state with respect to seconary body by differencing the two
X_SCI = X_BaCI - X_secondary_BaCI;

end % function