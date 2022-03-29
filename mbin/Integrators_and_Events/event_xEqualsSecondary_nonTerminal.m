function [value, isterminal, direction] = event_xEqualsSecondary_nonTerminal(t,X,prms)
value      = X(1) - (1-prms.u); % When 'x' == 1-mu
isterminal = 0; % stops the integration
direction  = 0; % negative direction only
end

