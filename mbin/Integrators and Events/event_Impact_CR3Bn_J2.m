function [value, isterminal, direction] = event_Impact_CR3Bn_J2(t,X,u,R1,R2,J21,J22)
%%% Designed for CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R2_n)
prms.u = u;
prms.R2_n = R2;
%%% Distances to primary (1) and secondary (2) bodies
r2 = sqrt((X(1)+prms.u-1)^2 + X(2)^2 + X(3)^2);

value = r2 - prms.R2_n; % When the surface is impacted
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end