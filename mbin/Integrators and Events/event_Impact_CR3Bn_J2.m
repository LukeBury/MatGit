function [value, isterminal, direction] = event_Impact_CR3Bn_J2(t,X,u,R1,R2,J21,J22)
%%% Designed for normalized CR3BP with J2
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system
%          R1 - radius of primary body
%          R2 - radius of secondary body
%          J21 - J2 of primary body
%          J22 - J2 of secondary body

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position

%%% Distances to primary (1) and secondary (2) bodies
r2 = sqrt((x+u-1)^2 + y^2 + z^2);

value = r2 - R2; % When the surface is impacted
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end