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
%%% Problem 2
% ========================================================================
% ------------------------------------
%%% Setup
% ------------------------------------
%%% General parameters
AU_km = 1.495978707e8; % km

%%% Transfer orbit parameters
rPer_trans = 2.17 * AU_km; % km
rApo_trans = 2.57 * AU_km; % km
a_trans = (rPer_trans + rApo_trans) / 2; % km
e_trans = (rApo_trans - rPer_trans) / (rApo_trans + rPer_trans);
Tp_trans = 2*pi*sqrt((a_trans^3)/bodies.sun.u);
meanMot_trans = 2*pi/Tp_trans;

%%% Ceres orbit parameters
Tp_Ceres = 1682 * 86400;
e_Ceres = 0.0758;
a_Ceres = (bodies.sun.u * (Tp_Ceres/(2*pi))^2)^(1/3); % km
rPer_Ceres = a_Ceres * (1 - e_Ceres); % km

%%% Need true anomaly at Ceres arrival
ta_CeresArrival = acos(((a_trans*(1-e_trans^2))/rPer_Ceres - 1) / e_trans) % rad

%%% Need to turn this true anomaly into a delta-T
E_CeresArrival = acos((e_trans + cos(ta_CeresArrival))/(1 + e_trans*cos(ta_CeresArrival))); % rad
M_CeresArrival = E_CeresArrival - e_trans*sin(E_CeresArrival) % rad

transferTime_sec  = M_CeresArrival / meanMot_trans; % sec
transferTime_days = transferTime_sec/86400 % days




% %%% Date 1 (Vesta Departure)
% year_1   = 2012;
% month_1  = 9;
% day_1    = 5;
% hour_1   = 0;
% minute_1 = 0;
% second_1 = 0;
% [ JD_1 ] = calendar2julian(year_1, month_1, day_1, hour_1, minute_1, second_1);


% ------------------------------------
%%% 
% ------------------------------------

















