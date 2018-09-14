clear
clc

% Given
y = [1;2;1];
R = [.5 0 0; 0 1 0; 0 0 1];
H = [1;1;1];
xbar = 2;
Pbar = .5;
W = inv(R);

% Batch
xhatBATCH = inv(inv(Pbar) + H'*W*H)*(H'*W*y + inv(Pbar)*xbar)

% Covariance (from batch information method)
P = inv(inv(Pbar) + H'*W*H)




