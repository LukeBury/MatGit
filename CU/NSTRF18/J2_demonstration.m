clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin')
addpath('/Users/lukebury/Documents/MATLAB/CU/bin/LowEnergy')
addpath('/Users/lukebury/Documents/MATLAB/CU/bin')
colors = get_color_palettes();

syms rsb2 rb2b1 Rsb1 ub1 ub2 Rb1 J2

% rsb2 = [xsb2; ysb2; zsb2];     
% rb2b1 = [xb2b1; yb2b1; zb2b1]; 
% Rsb1 = rsb2 + rb2b1;
rsb2 = [xsb2; ysb2; zsb2];     
rb2b1 = [xb2b1; yb2b1; zb2b1]; 
Rsb1 = rsb2 + rb2b1;

U_2B = ub2/rsb2;

U_3B = ub1/rb2b1 - ub1/rsb1;

% U_J2 = -(3*u*RE*RE*J2*z*z)/(2*(r^5)) + (u*RE*RE*J2)/(2*(r^3))
U_J2_noZ = (ub1*Rb1*Rb1*J2)/(2*(r^3))




% U(x,y) = -(1-u)/r1 - u/r2
% where r1 = dist to primary and r2 = dist to secondary






Utot_uJ2J3 = -u/r + U_J2J3; % Total potential

EQM = [dx; dy; dz; diff(Utot_uJ2J3, x); diff(Utot_uJ2J3, y); diff(Utot_uJ2J3,z)];