function [dX] = integrator_FP_nn(t,X,rad2,uJ,uEur,aEur,wEur)
% ------------------------------------------------------------------------
%%% Dynamics
% ------------------------------------------------------------------------
%%% Unpack the Hopper state vector
rsE = X(1:3); % Hopper Position
vsE = X(4:6); % Hopper Velocity

%%% Find current rotation angle
th = t*wEur(3); % rad/s

%%% Finding Jupiter Position
rJupiter = R3([-aEur;0;0],th); % km
rEJ = -rJupiter; % km

%%% Finding spacecraft position relative to jupiter
rsJ = rEJ + rsE; % km

a2B = -uEur*rsE/(norm(rsE)^3); % km/s^2
a3B = -uJ*(rsJ/(norm(rsJ)^3) - (rEJ)/(norm(rEJ)^3)); % km/s^2

%%% Equations of Motion (CR3BP)
ddx = a2B(1) + a3B(1);
ddy = a2B(2) + a3B(2);
ddz = a2B(3) + a3B(3);

%%% Storing EOMs
ddrH = [ddx; ddy; ddz];

%%% Output the derivative of the state
dX = zeros(6+6^2,1);
dX(1:3) = vsE; % km/s
dX(4:6) = ddrH; % km/s^2

% ------------------------------------------------------------------------
%%% State Transition Matrix
% ------------------------------------------------------------------------
%%% Assigning necessary vector components
xsE = X(1); ysE = X(2); zsE = X(3); 
xsJ = rsJ(1); ysJ = rsJ(2); zsJ = rsJ(3);

%%% Build A matrix and evaluate at current state
A = zeros(6,6);
A(1:3,4:6) = eye(3,3);
% Asym(4,1)
A(4,1) = (3*uEur*xsE*abs(xsE)*sign(xsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) - uJ*(1/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(3/2) - (3*abs(xsE - xsJ)*sign(xsE - xsJ)*(xsE - xsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2)) - uEur/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(3/2);
% Asym(4,2)
A(4,2) = (3*uEur*xsE*abs(ysE)*sign(ysE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(ysE - ysJ)*sign(ysE - ysJ)*(xsE - xsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2);
% Asym(4,3)
A(4,3) = (3*uEur*xsE*abs(zsE)*sign(zsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(zsE - zsJ)*sign(zsE - zsJ)*(xsE - xsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2);
% Asym(5,1)
A(5,1) = (3*uEur*ysE*abs(xsE)*sign(xsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(xsE - xsJ)*sign(xsE - xsJ)*(ysE - ysJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2);
% Asym(5,2)
A(5,2) = (3*uEur*ysE*abs(ysE)*sign(ysE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) - uJ*(1/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(3/2) - (3*abs(ysE - ysJ)*sign(ysE - ysJ)*(ysE - ysJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2)) - uEur/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(3/2);
% Asym(5,3)
A(5,3) = (3*uEur*ysE*abs(zsE)*sign(zsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(zsE - zsJ)*sign(zsE - zsJ)*(ysE - ysJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2);
% Asym(6,1)
A(6,1) = (3*uEur*zsE*abs(xsE)*sign(xsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(xsE - xsJ)*sign(xsE - xsJ)*(zsE - zsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2);
% Asym(6,2)
A(6,2) = (3*uEur*zsE*abs(ysE)*sign(ysE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(ysE - ysJ)*sign(ysE - ysJ)*(zsE - zsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2);
% Asym(6,3)
A(6,3) = (3*uEur*zsE*abs(zsE)*sign(zsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) - uJ*(1/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(3/2) - (3*abs(zsE - zsJ)*sign(zsE - zsJ)*(zsE - zsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2)) - uEur/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(3/2);

%%% Reshape (n^2,1) stm to (n,n)
stm = reshape(X(7:end),6,6);

%%% Calculate new STM
stm_dot = A*stm;

%%% Output the derivative of the state
dX(7:end) = reshape(stm_dot,36,1);
end

