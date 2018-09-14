clear
clc
% Given
xbar0 = [1;1;1];
Pbar0 = [4 0 0; 0 2 0; 0 0 1];
Htilde1 = [1 0 0];
y1 = 2;
phi_t_t0 = [1 1 .5; 0 1 1; 0 0 1];
R = 1;

% Batch
H1 = Htilde1*phi_t_t0;
xhatBATCH = inv(inv(Pbar0) + H1'*inv(R)*H1)*(H1'*inv(R)*y1 + inv(Pbar0)*xbar0)

% CKF
xbar1 = phi_t_t0*xbar0;
Pbar1 = phi_t_t0*Pbar0*phi_t_t0';
K1 = Pbar1*Htilde1'*inv(Htilde1*Pbar1*Htilde1' + R);
xhat1 = xbar1 + K1*(y1 - Htilde1*xbar1);
xhatCKF = inv(phi_t_t0)*xhat1