clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================


% ========================================================================
%%% Work
% ========================================================================
A = [0, 1; 0, 0];
S = eye(2);
Q = zeros(2,2);
R = 1;
x0 = [1; 1];

B = [0; 1];

tf = 10;


M = [A, -B*inv(R)*B';...
    -Q, -A']

syms lam21 lam22 real
expm(M*tf)*[x0; -lam21; -lam22]
clear lam21 lam22

vec = @(lam21, lam22) expm(M*tf)*[x0; -lam21; -lam22]


lam21 = -213/3533;
lam22 = -1483/3533;


y = vec(lam21,lam22)

xf = y(1:2)


%%% verifying solution
lam = [lam21; lam22];

-lam'*A*x0 - .5*lam'*B*R*B'*lam - xf'*S*A*xf + .5*xf'*S*B*R*B'*S*xf

expm(M*tf)*[x0; -lam] - [xf;xf]





















