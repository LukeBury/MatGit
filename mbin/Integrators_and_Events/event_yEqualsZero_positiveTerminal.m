function [value, isterminal, direction] = event_yEqualsZero_positiveTerminal(t,X,prms)
value      = X(2); % When 'y' == 0
isterminal = 1; % stops the integration
direction  = 1; % negative direction only
end

