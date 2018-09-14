function [dX] = integrator_FPdata(t,X,rad2,uJ,uEur,aEur,wEur)
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
dX = zeros(6,1);
dX(1:3) = vsE; % km/s
dX(4:6) = ddrH; % km/s^2
end

