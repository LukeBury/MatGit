function [value, isterminal, direction] = event_zEqualsZero_inNeck(t,X,prms)
value      = [X(3), X(1) - prms.xStop]; % When 'z' == 0
isterminal = [0, 1]; % stops the integration
direction  = [0, 0]; % negative direction only
end

