function [value, isterminal, direction] = impact_FP_nn(t,X,rad2,uJ,uEur,aEur,wEur)
%%% Event function watching for when "value" = 0
% r2 = sqrt((X(1)-rB2(1))^2 + (X(2)-rB2(2))^2 + (X(3)-rB2(3))^2);
rEurCI = norm(X(1:3));
value = rEurCI - rad2; % When the surface is impacted
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end