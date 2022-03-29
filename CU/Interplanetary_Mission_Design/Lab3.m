clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Plotting
% ========================================================================

transferDuration = [150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250];

C3_values = [20.49, 18.49, 17.15, 16.31, 15.89, 15.86, 16.27, 17.29, 19.32, 23.45, 32.97];

vInf_values = [5.30, 4.66, 4.13, 3.68, 3.32, 3.05, 2.87, 2.80, 2.86, 3.12, 3.73];


figure; hold all
plot(transferDuration, C3_values,'--o','linewidth',2,'markersize',8,'markerfacecolor',colors.std.blue,'markeredgecolor',colors.std.black)
PlotBoi2('Transfer Duration, $Days$','C$_3$ at Earth Departure, $km^2/s^2$',18,'LaTex')

figure; hold all
plot(transferDuration, vInf_values,'--o','linewidth',2,'markersize',8,'markerfacecolor',colors.std.blue,'markeredgecolor',colors.std.black)
PlotBoi2('Transfer Duration, $Days$','$v_\infty$ at Mars arrival, $km/s$',18,'LaTex')



