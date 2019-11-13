function [dX] = Int_2BI_constPert(t,X,u,aPert)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity and a constant perturbation vector
%%% Inputs:
%          t - time vector
%          X - [6x1] initial state
%          u - gravitational parameter of primary body
%          aPert - [3x1] perturbing acceleration

%%% Preallocate state output
dX = zeros(6,1);

%%% Distances from primary body to spacecraft
r = [X(1); X(2); X(3)]; % km
r_mag = norm(r);

%%% 2B accelerations
a_2B = r * (-u / (r_mag^3));

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)];    % km/s
dX(4:6) = [a_2B(1) + aPert(1); a_2B(2) + aPert(2); a_2B(3) + aPert(3)]; % km/s^2

end
