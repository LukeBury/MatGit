clear
clc
close all
A = [-0.0188 11.5959 0 -32.2;...
    -0.0007 -0.5357 1 0;...
    0.000048 -0.4944 -0.4935 0;...
    0 0 1 0];
B = [0; 0; -0.5632; 0];
% ------------------------------------------------------------------------
%%% 4b
% ------------------------------------------------------------------------
fprintf('======================================================================\n')
fprintf('================================= 4b =================================')
eigs = eig(A);
[frequency,d,poles] = damp(A);
poles
frequency % rad/s

fprintf('Norm of each pole of A: [%1.4f, %1.4f, %1.4f, %1.4f]''\n'...
    ,norm(eigs(1)),norm(eigs(2)),norm(eigs(3)),norm(eigs(4)))
fprintf('The system is stable since the magnitude of all poles are < 1\n\n')
fprintf('The system''s phugoid pole is -0.522 +- 0.703i\n')
fprintf('The system''s short period pole is -0.00193 +- 0.125i\n')
% ------------------------------------------------------------------------
%%% 4c
% ------------------------------------------------------------------------
fprintf('======================================================================\n')
fprintf('================================= 4c =================================\n')
fprintf('x(k+1) = F*x(k) + G*u(k)')
dt = 0.01;  % s

C = eye(4,4);
D = zeros(4,1);

sys = ss(A,B,C,D);
sysd = c2d(sys,dt,'ZOH');

F = sysd.A
G = sysd.B
H = sysd.C;
M = sysd.D;

fprintf('======================================================================\n')
fprintf('================================= 4d =================================')

[frequency,d,poles] = damp(F);
poles
frequency % rad/s

fprintf('The system''s phugoid pole is 0.995 +- 0.00699i because\n')
fprintf('   it has the lowest natural frequency.\n')
fprintf('The system''s short period pole is 0.9999 +- 0.0012i\n')
fprintf('   because it has the highest natural frequency')



