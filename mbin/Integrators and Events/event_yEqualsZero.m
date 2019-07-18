function [value, isterminal, direction] = event_yEqualsZero(t,X,prms)
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R2_n)

value      = X(2); % When 'y' == 0
isterminal = 0; % stops the integration
direction  = 0; % negative direction only
end