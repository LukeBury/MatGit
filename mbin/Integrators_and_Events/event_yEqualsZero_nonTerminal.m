function [value, isterminal, direction] = event_yEqualsZero_nonTerminal(t,X,prms)
value      = X(2); % When 'y' == 0
isterminal = 0; % stops the integration
direction  = 0; % negative direction only
end

