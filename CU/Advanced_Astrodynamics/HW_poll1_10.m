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
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% 
% ------------------------------------
% syms p1 p2 q1 q2 P1 P2 Q1 Q2 real
syms p1P1 p2P1 q1P1 q2P1 real
syms p1P2 p2P2 q1P2 q2P2 real
syms p1Q1 p2Q1 q1Q1 q2Q1 real
syms p1Q2 p2Q2 q1Q2 q2Q2 real

mat = [p1P1 p2P1 q1P1 q2P1;...
       p1P2 p2P2 q1P2 q2P2;...
       p1Q1 p2Q1 q1Q1 q2Q1;...
       p1Q2 p2Q2 q1Q2 q2Q2]'


% M = [cos(Q2), -(Q1^-1)*sin(Q2), P2*(Q1^-2)*sin(Q2), -P1*sin(Q2)-P2*(Q1^-1)*cos(Q2);...
%      sin(Q2), (Q1^-1)*cos(Q2), -P2*(Q1^-2)*cos(Q2), P1*cos(Q2)-P2*(Q1^-1)*sin(Q2);...
%      0, 0, cos(Q2), -Q1*sin(Q2);...
%      0, 0, sin(Q2), Q1*cos(Q2)];

J = [zeros(2,2), eye(2,2);...
     -eye(2,2), zeros(2,2)];
result = mat' * J * mat;

result(1:2,1:2)
result(1:2,3:4)
result(3:4,1:2)
result(3:4,3:4)






































