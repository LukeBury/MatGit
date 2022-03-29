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
hessian = @(x1,x2,x3,l1,l2) ...
    [2 + 2*l1 + 2*l2, 0, 0, 2*x1, 2*x1;...
    0, 2 + 2*l1 + 2*l2, 0, 2*x2, 2*x2;...
    0, 0, 2 + 2*l1, 2*x3, 0;...
    2*x1, 2*x2, 2*x3, 0, 0;...
    2*x1, 2*x2, 0, 0, 0];

x1 = 0; x2 = 1; x3 = 1; l1 = -1; l2 = 0;
h1 = hessian(x1, x2, x3, l1, l2);

x1 = 0; x2 = -1; x3 = 1; l1 = -1; l2 = 0;
h2 = hessian(x1, x2, x3, l1, l2);

x1 = 1; x2 = 0; x3 = 1; l1 = -1; l2 = 2;
h3 = hessian(x1, x2, x3, l1, l2);

x1 = -1; x2 = 0; x3 = 1; l1 = -1; l2 = -2;
h4 = hessian(x1, x2, x3, l1, l2);


[v1, e1] = eig(h1);
[v2, e2] = eig(h2);
[v3, e3] = eig(h3);
[v4, e4] = eig(h4);

diag(e1)'
diag(e2)'
diag(e3)'
diag(e4)'




% 
% x1 = 0; x2 = 1; x3 = -1; l1 = -1; l2 = 0;
% h1 = hessian(x1, x2, x3, l1, l2);
% 
% x1 = 0; x2 = -1; x3 = -1; l1 = -1; l2 = 0;
% h2 = hessian(x1, x2, x3, l1, l2);
% 
% x1 = 1; x2 = 0; x3 = -1; l1 = -1; l2 = 2;
% h3 = hessian(x1, x2, x3, l1, l2);
% 
% x1 = -1; x2 = 0; x3 = -1; l1 = -1; l2 = -2;
% h4 = hessian(x1, x2, x3, l1, l2);
% 
% 
% [v1, e1] = eig(h1);
% [v2, e2] = eig(h2);
% [v3, e3] = eig(h3);
% [v4, e4] = eig(h4);
% 
% diag(e1)'
% diag(e2)'
% diag(e3)'
% diag(e4)'









