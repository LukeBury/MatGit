function [value, isterminal, direction] = event_zEqualsZero(t,X,prms)
value      = X(3); % When 'z' == 0
isterminal = 1; % stops the integration
direction  = 1; % negative direction only
end

