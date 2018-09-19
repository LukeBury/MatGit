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
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Running
% ========================================================================

syms a e real

% eqnOG = sqrt(2*a2/(a1*(a1+a2))) + sqrt((2*a1)/(a2*(a1+a2))) - 2/sqrt(a2) + sqrt(2/(a1*(1+a1)) + sqrt(2/(a2*(a+a2)));

eqnFull = -sqrt((2+2*e)/(2+3*a*e+e)) + sqrt((2+2*a*e)/(2+a*e+3*e)) - 2/sqrt(1+e) + sqrt(2/(2 + 3*a*e)) + sqrt(2/(2+3*e));
pretty(eqnFull)

e = 0;
a = 3;
vpa(subs(eqnFull))

e = 0.5;
a = 3;
vpa(subs(eqnFull))

e = 0;
a = 0;
vpa(subs(eqnFull))




