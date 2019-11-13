function [value, isterminal, direction] = event_Impact_CR3Bn(t,X,prms)
%%% Designed for CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R2_n)

%%% Distances to primary (1) and secondary (2) bodies
r2 = sqrt((X(1)+prms.u-1)^2 + X(2)^2 + X(3)^2);

if isfield(prms,'R2') == 0
    warning('Need field of prms.R2')
end
value = r2 - prms.R2; % When the surface is impacted
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end