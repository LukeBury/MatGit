function [dX] = Int_2BI(t,X,u)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body

%%% Preallocate state output
dX = zeros(6,1);

%%% Distances from primary body to spacecraft
r = [X(1); X(2); X(3)]; % km
r_mag = norm(r);

%%% 2B accelerations
a_2B = r * (-u / (r_mag^3));

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)];    % km/s
dX(4:6) = [a_2B(1); a_2B(2); a_2B(3)]; % km/s^2

end
