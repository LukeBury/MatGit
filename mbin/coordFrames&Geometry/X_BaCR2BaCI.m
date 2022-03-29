function [X_BaCI] = X_BaCR2BaCI(X_BaCR, times, prms)
%%% Description
% Converts a bary-centered rotating 3-Body state to a bary-centered
% inertial state
% NOTE: ALL NORMALIZED UNITS
% ------------------------------------------------------------------------
%%% Inputs
% X_BaCR - [nx6] Array of state values in bary-centered-rotating frame
% times  - [nx1] array of times associated with the states
% prms   - [structure] Required fields: u, n
% ------------------------------------------------------------------------
%%% Outputs
% X_BaCI - [nx6] Array of state values in bary-centered-inertial frame
% ------------------------------------------------------------------------
% Created: 8/18/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Convert bary-centered rotating state to bary-centered inertial
% -------------------------------------------------
%%% Preallocate
X_BaCI = zeros(size(X_BaCR));


%%% If there's only 1 state (1 time)
if length(times) == 1
    %%% Compute r_BaCI
    X_BaCI(1:3) = R3_rot(X_BaCR(1:3), times.*prms.n);

    %%% Compute v_BaCI
    X_BaCI(4:6) = R3_rot(X_BaCR(4:6), times.*prms.n) + cross([0,0,1], X_BaCI(1:3));

% If there are more than 1 states
else
    for X_i = 1:length(times)
%         % Compute r_BaCI
%         X_BaCI(1:3,X_i) = R3_rot(X_BaCR(1:3,X_i), times(X_i).*prms.n);
% 
%         % Compute v_BaCI
%         X_BaCI(4:6,X_i) = R3_rot(X_BaCR(4:6,X_i), times(X_i).*prms.n) + cross([0,0,1], X_BaCI(1:3,X_i));
        % Compute r_BaCI
        X_BaCI(X_i, 1:3) = R3_rot(X_BaCR(X_i, 1:3), times(X_i).*prms.n);

        % Compute v_BaCI
        X_BaCI(X_i, 4:6) = R3_rot(X_BaCR(X_i, 4:6), times(X_i).*prms.n) + cross([0,0,1], X_BaCI(X_i, 1:3));
    end
end


% # -------------------------------------------------
% # Convert bary-centered rotating state to bary-centered inertial
% # -------------------------------------------------
% # Preallocate
% X_BaCI = zeros(size(X_BaCR))
% 
% # If there's only 1 state (1 time)
% if length(times) == 1
%     ### Compute r_BaCI
%     X_BaCI[1:3] = Rot3(X_BaCR[1:3], times)
% 
%     ### Compute v_BaCI
%     X_BaCI[4:6] = Rot3(X_BaCR[4:6], times) + cross([0,0,1], X_BaCI[1:3])
% 
% # If there are more than 1 states
% else
%     for X_i in 1:length(times)
%         # Compute r_BaCI
%         X_BaCI[1:3,X_i] = Rot3(X_BaCR[1:3,X_i], times[X_i])
% 
%         # Compute v_BaCI
%         X_BaCI[4:6,X_i] = Rot3(X_BaCR[4:6,X_i], times[X_i]) + cross([0,0,1], X_BaCI[1:3,X_i])
%     end
% end








% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
% 
% --------------------------


end % function