function [X_p] = rk4_MRPw(X,dt,I,u)
%%%
%%% Inputs:
%           1) X  - State, [6x1]
%           2) dt - step size, [1x1]
%           3) I  - Inertia tensor, [3x3]
%           4) u  - control input, [3x1]
%%% Outputs:
%           1) X_p - State at time t+dt, [nx1]
% ========================================================================
%%% EOM
dXdt = @(s,w,I,u) [0.25 * [1-norm(s)^2+2*s(1)^2, 2*(s(1)*s(2)-s(3)), 2*(s(1)*s(3)+s(2))] * w;... % dMRP-1
                 0.25 * [2*(s(2)*s(1)+s(3)), 1-norm(s)^2+2*s(2)^2, 2*(s(2)*s(3)-s(1))] * w;...   % dMRP-2
                 0.25 * [2*(s(3)*s(1)-s(2)), 2*(s(3)*s(2)+s(1)), 1-norm(s)^2+2*s(3)^2] * w;...   % dMRP-3
                 -(I(3,3)-I(2,2))*w(2)*w(3)/I(1,1) + u(1)/I(1,1);...                                    % dw-1
                 -(I(1,1)-I(3,3))*w(1)*w(3)/I(2,2) + u(2)/I(2,2);...                                    % dw-2
                 -(I(2,2)-I(1,1))*w(1)*w(2)/I(3,3) + u(3)/I(3,3)];                                      % dw-3

%%% Runge Kutta Time!
k1 = dXdt(X(1:3),X(4:6),I,u);
X1 = X + k1.*(dt/2);
k2 = dXdt(X1(1:3),X1(4:6),I,u);
X2 = X + k2.*(dt/2);
k3 = dXdt(X2(1:3),X2(4:6),I,u);
X3 = X + k3.*dt;
k4 = dXdt(X3(1:3),X3(4:6),I,u);

X_p = X + (k1 + 2*k2 + 2*k3 + k4)*(dt/6);

end









