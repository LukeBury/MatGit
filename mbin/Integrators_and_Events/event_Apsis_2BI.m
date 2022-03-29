function [value, isterminal, direction] = event_Apsis_2BI(t,X,u)
%%% Designed for standard normalized CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - time vector
%          X - state [6x1]
%          u - gravitational parameter of primary

%%% Unpack the barycentric state vector
r = X(1:3); % Position
v = X(4:6); % Velocity

%%% When velocity is perpendicular to position, dot(r,v) = 0
value = dot(r,v);

isterminal = 1; % stops the integration
direction = 0; % 0 - all directions
end