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
%%% Earth Section
% ------------------------------------
Earth = bodies.earth;
Sun = bodies.sun;
Saturn = bodies.saturn;

altitude = 300; % km
RSc0 = Earth.R + altitude; % km

VSc0 = sqrt(Earth.u / RSc0);

% ------------------------------------
%%% Transfer section
% ------------------------------------
rpTrans = Earth.a; % km
raTrans = Saturn.a; % km
aTrans = (rpTrans + raTrans) / 2; % km

vSc_rpTrans = visviva_v(rpTrans, aTrans, Sun.u);
vSc_raTrans = visviva_v(raTrans, aTrans, Sun.u);

TpTrans = pi*sqrt((aTrans^3)/Sun.u); % sec

% ------------------------------------
%%% Saturn section
% ------------------------------------
RSc_SOI = 20000 + Saturn.R;

% V2_TVI = sqrt(VInf_Launch^2 + 2*Earth.u/(Earth.R+300));
VInf_SOI = sqrt(Sun.u / Saturn.a) - vSc_raTrans;

VSc_SOI_m = sqrt(VInf_SOI^2 + 2*Saturn.u/RSc_SOI);

days_per_year = 365.242189;
Tp_SaturnOrbit = days_per_year*0.5*86400;
a_SaturnOrbit = (((Tp_SaturnOrbit/(2*pi))^2)*Saturn.u)^(1/3);
rp_SaturnOrbit = RSc_SOI; % From Cassini

VSc_SOI_p = visviva_v(rp_SaturnOrbit, a_SaturnOrbit, Saturn.u);
% ------------------------------------
%%% DVs
% ------------------------------------
VInf_TVI = vSc_rpTrans - sqrt(Sun.u / Earth.a);
VSc_TVI = sqrt(VInf_TVI^2 + 2*Earth.u/RSc0);

DV_TVI = VSc_TVI - VSc0
DV_SOI = VSc_SOI_m - VSc_SOI_p





[DV1, DV2, DVTot, TpTrans_sec] = interplanetaryHohmannTransfer_DVs(Sun.u, Earth.u, Saturn.u, Earth.a, Saturn.a, RSc0, rp_SaturnOrbit, VSc0, VSc_SOI_p)























