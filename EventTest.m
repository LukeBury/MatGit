clear
clc
close all

r0 = [0;0];
v0 = [10;30];
X0 = [r0; v0];
ti = 0;
tf = 10;
time = ti:tf;
dummyVar1 = 10;
dummyVar2 = 22;
tol = 1E-6;
% options1 = odeset('Events',@(y)impactEvent(t,Y,E_radius,uE,uJ,J_pos,nE));%,'RelTol',tol,'AbsTol',tol);%%%%%%%%%%%%%%%%%
options = odeset('Events',@myEvent,'RelTol',tol,'AbsTol',tol);

[Times,States] = ode45(@EventTestIntegrator,time,X0,options,dummyVar2);

plot(States(:,1),States(:,2))

function [ dY ] = EventTestIntegrator(t,Y,dummyVar2)
dY = zeros(4,1);

%%% Unpack the state vector
y = Y(1:2);
dy = Y(3:4);
dummyVar2

% %%% Ensure we haven't crashed into Europa
% if sqrt(rH(1)^2 + rH(2)^2 + rH(3)^2) < E_radius-1
%     sqrt(rH(1)^2 + rH(2)^2 + rH(3)^2);
%     if T_impact == 0
%         T_final = t;
%         T_impact = t;
%     end
%     if t < T_impact
%         T_final = t;
%         T_impact = t;
%     end
%     return
% end

%%% Dynamics
ddy = [0; -9.8];
% ddy = (-uJ/(norm(yE)^3))*yE; % Jupiter pull, km/s^2

%%% Output the derivative of the state
dY(1:2) = dy;
dY(3:4) = ddy;
end

function [value, isterminal, direction] = myEvent(t,Y,dummyVar1)
%%% Unpack thestate vector
y = Y(1:2);
dummyVar1
value = y(2);
isterminal = 1; % stops the integration
direction = -1; % negative direction only

end