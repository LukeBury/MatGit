function [A] = get_Amat_CR3BP(mu, r_n, n)
%%% Description
%       Returns the A-matrix for the classical CR3BP
%       
% --------------------------------------------------------------
%%% Inputs
%       mu  - [1x1] Mass ratio of CR3BP system
%       r_n - [3x1] Normalized position vector to evaluate dynamics at
%       n   - [1x1] Normalized mean motion of system
% --------------------------------------------------------------
%%% Outputs
%       A - [6x6] A matrix for the classical CR3BP (Jacobian of state wrt dynamics)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Unpack position
% ------------------------------------
%%% Enforcing column vector
r_n = r_n(:);

%%% Unpacking
x_n = r_n(1);
y_n = r_n(2);
z_n = r_n(3);

% ------------------------------------
%%% Calculate A Matrix
% ------------------------------------
% %%% Full A matrix
% Amat = @(u, x, y, z) ...
%     [zeros(3,3),     eye(3);...
%      A_LL(u,x,y,z), [0,2,0;-2,0,0;0,0,0]];

A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4,5) = 2*n;
A(5,4) = -2*n;
A(4,1) = (mu - 1)/((mu + x_n)^2 + y_n^2 + z_n^2)^(3/2) - mu/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(3/2) + (3*mu*(2*mu + 2*x_n - 2)*(mu + x_n - 1))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2)) - (3*(2*mu + 2*x_n)*(mu + x_n)*(mu - 1))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) + n^2;
A(4,2) = (3*mu*y_n*(mu + x_n - 1))/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*y_n*(mu + x_n)*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2);
A(4,3) = (3*mu*z_n*(mu + x_n - 1))/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*z_n*(mu + x_n)*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2);
A(5,1) = -y_n*((3*(2*mu + 2*x_n)*(mu - 1))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) - (3*mu*(2*mu + 2*x_n - 2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2)));
A(5,2) = (mu - 1)/((mu + x_n)^2 + y_n^2 + z_n^2)^(3/2) + y_n*((3*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*y_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) - mu/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(3/2) + n^2;
A(5,3) = y_n*((3*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*z_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2));
A(6,1) = -z_n*((3*(2*mu + 2*x_n)*(mu - 1))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) - (3*mu*(2*mu + 2*x_n - 2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2)));
A(6,2) = z_n*((3*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*y_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2));
A(6,3) = (mu - 1)/((mu + x_n)^2 + y_n^2 + z_n^2)^(3/2) + z_n*((3*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*z_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) - mu/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(3/2);

end

