function [dX] = Int_CR3BnSTM(t,X,prms)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [42x1]
%          prms - (u)

%%% Preallocate state output
dX = zeros(6+36,1);

%%% Distances to primary (1) and secondary (2) bodies
% r1 = sqrt(((X(1)+(1-prms.u))+prms.u)^2 + X(2)^2 + X(3)^2);
% r2 = sqrt(((X(1)+(1-prms.u))+prms.u-1)^2 + X(2)^2 + X(3)^2);
% r1 = sqrt((X(1)+1)^2 + X(2)^2 + X(3)^2);
% r2 = sqrt(X(1)^2 + X(2)^2 + X(3)^2);

%%% Equations of Motion
% ddx = 2*prms.n*X(5) + (prms.n^2)*(X(1)+(1-prms.u)) - (1-prms.u)*((X(1)+(1-prms.u))+prms.u)/(r1^3) - prms.u*((X(1)+(1-prms.u))+prms.u-1)/(r2^3);
% ddy = -2*prms.n*X(4) + (prms.n^2)*X(2) -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(2);
% ddz = -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(3);
ddx = (prms.n^2*(2*X(1) - 2*prms.u + 2))/2 + 2*prms.n*X(5) + ((2*X(1) + 2)*(prms.u - 1))/(2*((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(3/2)) - (prms.u*X(1))/(X(1)^2 + X(2)^2 + X(3)^2)^(3/2);
ddy = prms.n^2*X(2) - 2*prms.n*X(4) + (X(2)*(prms.u - 1))/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(3/2) - (prms.u*X(2))/(X(1)^2 + X(2)^2 + X(3)^2)^(3/2);
ddz = (X(3)*(prms.u - 1))/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(3/2) - (prms.u*X(3))/(X(1)^2 + X(2)^2 + X(3)^2)^(3/2);

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

%%% Reshaping STM
stm = reshape(X(7:end),6,6);

%%% Evaluate A matrix ... from A = jacobian(EOM_CR3BP,state)
A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4,5) = 2*prms.n;
A(5,4) = -2*prms.n;
% A(4,1) = (prms.u - 1)/((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(3/2) + prms.n^2 - prms.u/((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(3/2) + (3*prms.u*(2*prms.u + 2*(X(1)+(1-prms.u)) - 2)^2)/(4*((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2)) - (3*(2*prms.u + 2*(X(1)+(1-prms.u)))^2*(prms.u - 1))/(4*((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2));
% A(4,2) = (3*prms.u*X(2)*(2*prms.u + 2*(X(1)+(1-prms.u)) - 2))/(2*((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2)) - (3*X(2)*(2*prms.u + 2*(X(1)+(1-prms.u)))*(prms.u - 1))/(2*((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2));
% A(4,3) = (3*prms.u*X(3)*(2*prms.u + 2*(X(1)+(1-prms.u)) - 2))/(2*((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2)) - (3*X(3)*(2*prms.u + 2*(X(1)+(1-prms.u)))*(prms.u - 1))/(2*((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2));
% A(5,1) = (3*prms.u*X(2)*(2*prms.u + 2*(X(1)+(1-prms.u)) - 2))/(2*((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2)) - (3*X(2)*(2*prms.u + 2*(X(1)+(1-prms.u)))*(prms.u - 1))/(2*((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2));
% A(5,2) = (prms.u - 1)/((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(3/2) + prms.n^2 - prms.u/((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(3/2) - (3*X(2)^2*(prms.u - 1))/((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2) + (3*prms.u*X(2)^2)/((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2);
% A(5,3) = (3*prms.u*X(2)*X(3))/((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(2)*X(3)*(prms.u - 1))/((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2);
% A(6,1) = (3*prms.u*X(3)*(2*prms.u + 2*(X(1)+(1-prms.u)) - 2))/(2*((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2)) - (3*X(3)*(2*prms.u + 2*(X(1)+(1-prms.u)))*(prms.u - 1))/(2*((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2));
% A(6,2) = (3*prms.u*X(2)*X(3))/((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(2)*X(3)*(prms.u - 1))/((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2);
% A(6,3) = (prms.u - 1)/((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(3/2) - prms.u/((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(3/2) - (3*X(3)^2*(prms.u - 1))/((prms.u + (X(1)+(1-prms.u)))^2 + X(2)^2 + X(3)^2)^(5/2) + (3*prms.u*X(3)^2)/((prms.u + (X(1)+(1-prms.u)) - 1)^2 + X(2)^2 + X(3)^2)^(5/2);
A(4,1) = prms.n^2 - prms.u/(X(1)^2 + X(2)^2 + X(3)^2)^(3/2) + (prms.u - 1)/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(3/2) - (3*(2*X(1) + 2)^2*(prms.u - 1))/(4*((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2)) + (3*prms.u*X(1)^2)/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2);
A(4,2) = (3*prms.u*X(1)*X(2))/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(2)*(2*X(1) + 2)*(prms.u - 1))/(2*((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2));
A(4,3) = (3*prms.u*X(1)*X(3))/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(3)*(2*X(1) + 2)*(prms.u - 1))/(2*((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2));
A(5,1) = (3*prms.u*X(1)*X(2))/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(2)*(2*X(1) + 2)*(prms.u - 1))/(2*((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2));
A(5,2) = prms.n^2 - prms.u/(X(1)^2 + X(2)^2 + X(3)^2)^(3/2) + (prms.u - 1)/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(3/2) - (3*X(2)^2*(prms.u - 1))/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2) + (3*prms.u*X(2)^2)/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2);
A(5,3) = (3*prms.u*X(2)*X(3))/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(2)*X(3)*(prms.u - 1))/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2);
A(6,1) = (3*prms.u*X(1)*X(3))/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(3)*(2*X(1) + 2)*(prms.u - 1))/(2*((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2));
A(6,2) = (3*prms.u*X(2)*X(3))/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2) - (3*X(2)*X(3)*(prms.u - 1))/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2);
A(6,3) = (prms.u - 1)/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(3/2) - prms.u/(X(1)^2 + X(2)^2 + X(3)^2)^(3/2) - (3*X(3)^2*(prms.u - 1))/((X(1) + 1)^2 + X(2)^2 + X(3)^2)^(5/2) + (3*prms.u*X(3)^2)/(X(1)^2 + X(2)^2 + X(3)^2)^(5/2);



%%% Calculate stmDot and output result
stmDot = A*stm;
dX(7:end) = reshape(stmDot,36,1);

end

