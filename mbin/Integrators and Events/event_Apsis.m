function [value, isterminal, direction] = event_Apsis(t,X,u,R2)
%%% Designed for standard normalized CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system
%          R2 - radius of secondary body

%%% Unpack the barycentric state vector
r = X(1:3); % Position
v = X(4:6); % Velocity

%%% Centering position on secondary
r(1) = r(1) - (1-u);

%%% When velocity is perpendicular to position, dot(r,v) = 0
value = dot(r,v);

isterminal = 0; % stops the integration
direction = 0; % 0 - all directions
end