function [dX_n] = integrator_FP(t_n,X_n,u,rB1_n,rB2_n,rad2_n,tNorm,rNorm,vNorm,JD0, uJ, uEur)
% ------------------------------------------------------------------------
%%% Dynamics
% ------------------------------------------------------------------------
%%% Unpack the Hopper state vector
rP_n = X_n(1:3); % Hopper Position
vP_n = X_n(4:6); % Hopper Velocity

%%% Assigning Variables
x_n = rP_n(1); y_n = rP_n(2); z_n = rP_n(3);
dx_n = vP_n(1); dy_n = vP_n(2); dz_n = vP_n(3);

%%% Creating hopper distances to bodies
r1_n = sqrt((rP_n(1)-rB1_n(1))^2 + (rP_n(2)-rB1_n(2))^2 + (rP_n(3)-rB1_n(3))^2);
r2_n = sqrt((rP_n(1)-rB2_n(1))^2 + (rP_n(2)-rB2_n(2))^2 + (rP_n(3)-rB2_n(3))^2);

%%% Equations of Motion (CR3BP)
ddx_CR3BP_n = 2*dy_n + x_n - (1-u)*(x_n+u)/(r1_n^3) - u*(x_n+u-1)/(r2_n^3);
ddy_CR3BP_n = -2*dx_n + y_n -((1-u)/(r1_n^3) + u/(r2_n^3))*y_n;
ddz_CR3BP_n = -((1-u)/(r1_n^3) + u/(r2_n^3))*z_n;

%%% Storing EOMs
ddxP_n = [ddx_CR3BP_n; ddy_CR3BP_n; ddz_CR3BP_n];

% ------------------------------------------------------------------------
%%% Converting States to Inertial non-normalized for STM purposes
% ------------------------------------------------------------------------
%%% Barycenter Inertial State (BCI)
wEur = 2.047200349303344e-05; % rad/s
t = t_n*tNorm; % sec
th = t * wEur; % rad
state_BCI = zeros(6,1);
state_BCI(1:3) = R3(X_n(1:3),th).*rNorm; % km
state_BCI(4:6) = R3(X_n(4:6),th).*vNorm + cross([0;0;wEur],state_BCI(1:3)); % km/s

% jupiter_BCI = R3([-u; 0; 0],th).*rNorm; % km
% europa_BCI = R3([1-u; 0; 0],th).*rNorm; % km
jupiter_BCI = R3(rB1_n',th).*rNorm; % km
europa_BCI = R3(rB2_n',th).*rNorm; % km

rsE = state_BCI(1:3) - europa_BCI; % km
rsJ = state_BCI(1:3) - jupiter_BCI; % km

xsE = rsE(1); ysE = rsE(2); zsE = rsE(3); % km
xsJ = rsJ(1); ysJ = rsJ(2); zsJ = rsJ(3); % km

% fprintf('------------\n')
% t_n*tNorm
% rEJ = rsJ-rsE;
% rsE
% a2B = -uEur*rsE/(norm(rsE)^3) % km/s^2
% scJ = -uJ*rsJ/(norm(rsJ)^3) 
% euJ = -uJ*(rEJ)/(norm(rEJ)^3)
% a3B = -uJ*(rsJ/(norm(rsJ)^3) - (rEJ)/(norm(rEJ)^3)) % km/s^2

% ------------------------------------------------------------------------
%%% State Transition Matrix
% ------------------------------------------------------------------------
%%% Reshape (n^2,1) stm to (n,n)
stm = reshape(X_n(7:end),6,6);

%%% Build A matrix and evaluate at current state
A = zeros(6,6);
A(1:3,4:6) = eye(3,3);
% Asym(4,1)
A(4,1) = ((3*uEur*xsE*abs(xsE)*sign(xsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) - uJ*(1/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(3/2) - (3*abs(xsE - xsJ)*sign(xsE - xsJ)*(xsE - xsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2)) - uEur/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(3/2));
% Asym(4,2)
A(4,2) = ((3*uEur*xsE*abs(ysE)*sign(ysE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(ysE - ysJ)*sign(ysE - ysJ)*(xsE - xsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2));
% Asym(4,3)
A(4,3) = ((3*uEur*xsE*abs(zsE)*sign(zsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(zsE - zsJ)*sign(zsE - zsJ)*(xsE - xsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2));
% Asym(5,1)
A(5,1) = ((3*uEur*ysE*abs(xsE)*sign(xsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(xsE - xsJ)*sign(xsE - xsJ)*(ysE - ysJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2));
% Asym(5,2)
A(5,2) = ((3*uEur*ysE*abs(ysE)*sign(ysE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) - uJ*(1/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(3/2) - (3*abs(ysE - ysJ)*sign(ysE - ysJ)*(ysE - ysJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2)) - uEur/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(3/2));
% Asym(5,3)
A(5,3) = ((3*uEur*ysE*abs(zsE)*sign(zsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(zsE - zsJ)*sign(zsE - zsJ)*(ysE - ysJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2));
% Asym(6,1)
A(6,1) = ((3*uEur*zsE*abs(xsE)*sign(xsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(xsE - xsJ)*sign(xsE - xsJ)*(zsE - zsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2));
% Asym(6,2)
A(6,2) = ((3*uEur*zsE*abs(ysE)*sign(ysE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) + (3*uJ*abs(ysE - ysJ)*sign(ysE - ysJ)*(zsE - zsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2));
% Asym(6,3)
A(6,3) = ((3*uEur*zsE*abs(zsE)*sign(zsE))/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(5/2) - uJ*(1/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(3/2) - (3*abs(zsE - zsJ)*sign(zsE - zsJ)*(zsE - zsJ))/(abs(xsE - xsJ)^2 + abs(ysE - ysJ)^2 + abs(zsE - zsJ)^2)^(5/2)) - uEur/(abs(xsE)^2 + abs(ysE)^2 + abs(zsE)^2)^(3/2));

%%% Calculate new STM
stm_dot = A*stm;

%%% Output the derivative of the state
dX_n = zeros(6+6^2,1);
dX_n(1:3) = vP_n; % km/s
dX_n(4:6) = ddxP_n; % km/s^2
dX_n(7:end) = reshape(stm_dot,36,1);
end

