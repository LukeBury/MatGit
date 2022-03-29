function [value, isterminal, direction] = event_zEqualsZero_nonTerminal(t,X,prms)
value      = X(3); % When 'z' == 0
isterminal = 0; % stops the integration
direction  = 0; % negative direction only
end

