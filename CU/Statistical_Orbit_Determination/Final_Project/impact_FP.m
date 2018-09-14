function [value, isterminal, direction] = impact_FP(t,X,u,rB1_n,rB2_n,rad2_n,tNorm,rNorm,vNorm,JD0,uJ,uEur)
%%% Event function watching for when "value" = 0
r2 = sqrt((X(1)-rB2_n(1))^2 + (X(2)-rB2_n(2))^2 + (X(3)-rB2_n(3))^2);

value = r2 - rad2_n; % When the surface is impacted
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end