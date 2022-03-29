function [value, isterminal, direction] = event_yEqualsZero(t,X,prms)
value      = X(2); % When 'y' == 0
isterminal = 1; % stops the integration
direction  = 0; % negative direction only
end

