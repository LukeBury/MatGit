function [A] = get_Amat_CR3BP_J2(mu, r_n, J2p, J2s, Rp, Rs)
%%% Description
%       Returns the A-matrix for the modified CR3BP with J2 for both the
%       primary and secondary bodies
%       
% --------------------------------------------------------------
%%% Inputs
%       mu  - [1x1] Mass ratio of CR3BP system
%       r_n - [3x1] Normalized position vector to evaluate dynamics at
%       J2p - [1x1] J2 value of primary body
%       J2s - [1x1] J2 value of secondary body
%       Rp  - [1x1] Normalized radius of primary body
%       Rs  - [1x1] Normalized radius of secondary body
% --------------------------------------------------------------
%%% Outputs
%       A - [6x6] A matrix for the modified CR3BP (Jacobian of state wrt dynamics)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Unpack position
% ------------------------------------
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
A(4,5) = 2;
A(5,4) = -2;
A(4,1) = (mu + x_n - 1)*((3*mu*(2*mu + 2*x_n - 2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2)) - (3*J2s*Rs^2*mu*(2*mu + 2*x_n - 2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2)) + (21*J2s*Rs^2*mu*(2*mu + 2*x_n - 2)*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(4*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2))) + (mu - 1)/((mu + x_n)^2 + y_n^2 + z_n^2)^(3/2) - (mu + x_n)*((3*(2*mu + 2*x_n)*(mu - 1))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) - (J2p*Rp^2*(2*mu + 2*x_n)*(3*mu - 3))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2)) + (7*J2p*Rp^2*(2*mu + 2*x_n)*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(4*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2))) - mu/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(3/2) - (3*J2s*Rs^2*mu*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2)) + (J2p*Rp^2*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2)) + 1;
A(4,2) = (mu + x_n - 1)*((3*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*J2s*Rs^2*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2) + (21*J2s*Rs^2*mu*y_n*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2))) - (mu + x_n)*((3*y_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2) - (J2p*Rp^2*y_n*(3*mu - 3))/((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2) + (7*J2p*Rp^2*y_n*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)));
A(4,3) = (mu + x_n - 1)*((3*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) + (12*J2s*Rs^2*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2) + (21*J2s*Rs^2*mu*z_n*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2))) - (mu + x_n)*((3*z_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2) + (4*J2p*Rp^2*z_n*(3*mu - 3))/((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2) + (7*J2p*Rp^2*z_n*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)));
A(5,1) = -y_n*((3*(2*mu + 2*x_n)*(mu - 1))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) - (3*mu*(2*mu + 2*x_n - 2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2)) - (J2p*Rp^2*(2*mu + 2*x_n)*(3*mu - 3))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2)) + (3*J2s*Rs^2*mu*(2*mu + 2*x_n - 2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2)) - (21*J2s*Rs^2*mu*(2*mu + 2*x_n - 2)*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(4*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2)) + (7*J2p*Rp^2*(2*mu + 2*x_n)*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(4*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)));
A(5,2) = (mu - 1)/((mu + x_n)^2 + y_n^2 + z_n^2)^(3/2) + y_n*((3*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*y_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2) + (J2p*Rp^2*y_n*(3*mu - 3))/((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2) - (3*J2s*Rs^2*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2) - (7*J2p*Rp^2*y_n*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)) + (21*J2s*Rs^2*mu*y_n*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2))) - mu/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(3/2) - (3*J2s*Rs^2*mu*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2)) + (J2p*Rp^2*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2)) + 1;
A(5,3) = y_n*((3*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*z_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2) - (4*J2p*Rp^2*z_n*(3*mu - 3))/((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2) + (12*J2s*Rs^2*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2) - (7*J2p*Rp^2*z_n*(3*mu - 3)*((mu + x_n)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)) + (21*J2s*Rs^2*mu*z_n*((mu + x_n - 1)^2 + y_n^2 - 4*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2)));
A(6,1) = -z_n*((3*(2*mu + 2*x_n)*(mu - 1))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2)) - (3*mu*(2*mu + 2*x_n - 2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2)) - (J2p*Rp^2*(6*mu + 6*x_n)*(3*mu - 3))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2)) + (3*J2s*Rs^2*mu*(6*mu + 6*x_n - 6))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2)) + (7*J2p*Rp^2*(2*mu + 2*x_n)*(3*mu - 3)*(3*(mu + x_n)^2 + 3*y_n^2 - 2*z_n^2))/(4*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)) - (21*J2s*Rs^2*mu*(2*mu + 2*x_n - 2)*(3*(mu + x_n - 1)^2 + 3*y_n^2 - 2*z_n^2))/(4*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2)));
A(6,2) = z_n*((3*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*y_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2) + (3*J2p*Rp^2*y_n*(3*mu - 3))/((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2) - (9*J2s*Rs^2*mu*y_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2) + (21*J2s*Rs^2*mu*y_n*(3*(mu + x_n - 1)^2 + 3*y_n^2 - 2*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2)) - (7*J2p*Rp^2*y_n*(3*mu - 3)*(3*(mu + x_n)^2 + 3*y_n^2 - 2*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)));
A(6,3) = z_n*((3*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(5/2) - (3*z_n*(mu - 1))/((mu + x_n)^2 + y_n^2 + z_n^2)^(5/2) - (2*J2p*Rp^2*z_n*(3*mu - 3))/((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2) + (6*J2s*Rs^2*mu*z_n)/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2) - (7*J2p*Rp^2*z_n*(3*mu - 3)*(3*(mu + x_n)^2 + 3*y_n^2 - 2*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(9/2)) + (21*J2s*Rs^2*mu*z_n*(3*(mu + x_n - 1)^2 + 3*y_n^2 - 2*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(9/2))) + (mu - 1)/((mu + x_n)^2 + y_n^2 + z_n^2)^(3/2) - mu/((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(3/2) + (J2p*Rp^2*(3*mu - 3)*(3*(mu + x_n)^2 + 3*y_n^2 - 2*z_n^2))/(2*((mu + x_n)^2 + y_n^2 + z_n^2)^(7/2)) - (3*J2s*Rs^2*mu*(3*(mu + x_n - 1)^2 + 3*y_n^2 - 2*z_n^2))/(2*((mu + x_n - 1)^2 + y_n^2 + z_n^2)^(7/2));

end

