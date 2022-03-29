function [dX] = Int_CR3Bn(t,X,prms)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u,n)

%%% Preallocate state output
dX = zeros(6,1);

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((X(1)+prms.u)^2 + X(2)^2 + X(3)^2);
r2 = sqrt((X(1)+prms.u-1)^2 + X(2)^2 + X(3)^2);

%%% Equations of Motion
ddx = 2*prms.n*X(5) + (prms.n^2)*X(1) - (1-prms.u)*(X(1)+prms.u)/(r1^3) - prms.u*(X(1)+prms.u-1)/(r2^3);
ddy = -2*prms.n*X(4) + (prms.n^2)*X(2) -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(2);
ddz = -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(3);

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2
end

