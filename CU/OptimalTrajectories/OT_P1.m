clear
clc
close all
addpath(genpath('/Users/CU_Google_Drive/lukebury/Documents/MATLAB/mbin'))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Running
% ========================================================================

% syms u r1 r2 r3 real
% assume(r2>r3)
% assume(r3>r2)
% 
% J_H = sqrt(2*u*r2/((r1+r2)*r1)) - sqrt(u/r1) + sqrt(u/r2) - sqrt(2*u*r1/((r1+r2)*r2));
% 
% J_Bi = sqrt(2*u*r3/((r1+r3)*r1)) - sqrt(u/r1) + sqrt(2*u*r2/((r2+r3)*r3)) - sqrt(2*u*r1/((r1+r3)*r3)) + sqrt(2*u*r3/((r2+r3)*r2)) - sqrt(u/r2);
% 
% normalizer = sqrt(u/r1);
% 
% Jt_H = J_H / normalizer;
% 
% Jt_Bi = J_Bi / normalizer;
% 
% r = r2/r1;
% l = r3/r1;


% --------------------------------
% notes = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1) - (sqrt(2)-1)*(1 + 1/sqrt(r));
% 
% 
% me = @(r) sqrt(2/(1+r)) + (1/r)*(2-sqrt(2)) - sqrt(2)/(r*sqrt(1+r)) - sqrt(2/r);
% 
% 
% % Target solution
% sol = @(r) r^3 - (7 + 4*sqrt(2))*r^2 + (3 + 4*sqrt(2))*r - 1;
% 
% 
% % these forms are equal for hohman
% Jt_H_HW    = @(r) sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)));
% Jt_H_Notes = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1);
% 
% Jt_P = @(r) (sqrt(2)-1)*(1 + 1/sqrt(r))

% --------------------------------

syms r real
assume(r>1)

inequality = (sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)))) - ((sqrt(2)-1)*(1 + 1/sqrt(r)));
inequality_f = @(r) (sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)))) - ((sqrt(2)-1)*(1 + 1/sqrt(r)));

x = 0:.001:20;
inequalitySolution = x.^3 - (7 + 4*sqrt(2)).*x.^2 + (3 + 4*sqrt(2)).*x - 1;

hold all
inSo = plot(x,inequalitySolution,'r','linewidth',1.5);


inequalityFromEquations = (sqrt(2.*x./(1+x)) - 1 + 1./sqrt(x) - sqrt(2./(x.*(1+x)))) - ((sqrt(2)-1).*(1 + 1./sqrt(x)));

% figure
inEq = plot(x,inequalityFromEquations,'b','linewidth',1.5);

plot([min(x) max(x)],[0 0],'k')

legend([inSo inEq], 'Given Solution','From Equations')
PlotBoi2('r','Cost',14,'LaTex')
ylim([-10 10])










