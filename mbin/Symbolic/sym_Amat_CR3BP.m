% ========================================================================
%%% Description
% ========================================================================
% Symbolically create the A matrix for the classical CR3BP

% Created: 11/19/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% 
% ========================================================================
clear
clc
close all

%%% Create symbollic variables and tell Matlab they're all real numbers
syms x y z xd yd zd u R1 R2 r1 r2 n real

%%% Build position vectors 
r1 = ((x+u)^2 + y^2 + z^2)^(1/2);
r2 = ((x-1+u)^2 + y^2 + z^2)^(1/2);

%%% Write acceleration terms
xdd = 2*n*yd + (n^2)*x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ydd = -2*n*xd + (n^2)*y -((1-u)/(r1^3) + u/(r2^3))*y;
zdd = -((1-u)/(r1^3) + u/(r2^3))*z;

%%% Create full equatio of motion vector (\dot{X})
EOM_CR3BP = [xd; yd; zd; xdd; ydd; zdd];

%%% Create full state vector (X)
state = [x;y;z;xd;yd;zd];

%%% Create A matrix by taking Jacobian
Amat_dFdX = simplify(jacobian(EOM_CR3BP, state))





