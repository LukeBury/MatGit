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

r = 0:.001:14;
eq1 = @(r) r.^3 - (7 + 4*sqrt(2)).*r.^2 + (3 + 4*sqrt(2)).*r - 1;
eq2 = @(r) r.^3 - (3 + 4*sqrt(2)).*r.^2 + (7 + 4*sqrt(2)).*r - 1;


figure; hold all
p1 = plot(r, eq1(r),'b','linewidth',1.5);
p2 = plot(r, eq2(r),'r','linewidth',1.5);
plot([min(r) max(r)],[0 0],'k')

PlotBoi2('$r$','Cost',14,'LaTex')
ylim([-10 10])
legend([p1 p2],'part 1','part 2')

eqn1_roots = roots([1, -(7 + 4*sqrt(2)),(3 + 4*sqrt(2)),-1])
eqn2_roots = roots([1, -(3 + 4*sqrt(2)),(7 + 4*sqrt(2)),-1])
% fzero(inequalitySolution,0)
% fzero(inequalitySolution,0.6)
% fzero(inequalitySolution,12)





