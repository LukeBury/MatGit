function [value, isterminal, direction] = event_ImpactEscape_CR3Bn(t,X,prms)
%%% Designed for standard normalized CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R2_n, L1x, L2x)

%%% Distances to primary (1) and secondary (2) bodies
r2 = sqrt((X(1)+prms.u-1)^2 + X(2)^2 + X(3)^2);

value = [r2 - prms.R2_n, X(1)-prms.L1x, prms.L2x-X(1)]; % When the surface is impacted
isterminal = [1, 1, 1]; % stops the integration
direction = [-1, -1, -1]; % negative direction only
end