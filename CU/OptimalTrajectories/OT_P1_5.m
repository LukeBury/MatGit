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
%%% l>r case (can be used for l<r case)
% ========================================================================
syms a v r u real
eq1 = sqrt(2*a/(1+a));
deq1 = diff(eq1,a);
pretty(simplify(deq1))
a = 1;
subs(deq1)

fprintf('==================\n')

syms a v r u real
eq2 = sqrt(2/a + v*v*r/u);
deq2 = diff(eq2,a);
pretty(simplify(deq2))
a = 1;
pretty(subs(deq2))

fprintf('==================\n')

syms a v r u real
eq3 = sqrt(2/(a*(1+a)));
deq3 = diff(eq3,a);
pretty(simplify(deq3))
a = 1;
subs(deq3)

















