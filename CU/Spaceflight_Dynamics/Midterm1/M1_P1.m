clear
clc
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
alt1 = 30; % km
alt2 = 250; % km
% Moon radius
Rm = 1740.0; % km
% Moon gravitational parameter
um = 0.01213 * 398600.4415; % km^3/s^2
% ------------------------------------------------------------------------
%%% Problem 1b
% ------------------------------------------------------------------------
%%% Calculating orbital info
% Transfer orbit info
rpt = Rm + alt1; % km
rat = Rm + alt2; % km
at = (rpt + rat)/2; % km

% Velocity at burnout point (periapsis of transfer)
vpt = sqrt(2*um/rpt - um/at); % km/s

%%% Printing Answer
fprintf('------------------ 1b ------------------\n')
fprintf('Velocity at burnout point = %f km/s\n\n',vpt)

% ------------------------------------------------------------------------
%%% Problem 1c
% ------------------------------------------------------------------------
%%% Calculating velocities for dV
% Velocity at transfer apoapsis
vat = sqrt(2*um/rat - um/at); % km/s

% Velocity of circular rendezvous orbit
vc = sqrt(um/rat); % km/s

% Delta V required
dV = vc - vat; % km/s

%%% Printing Answer
fprintf('------------------ 1c ------------------\n')
fprintf('Required dV = %f km/s\n\n',dV)

% ------------------------------------------------------------------------
%%% Problem 1d
% ------------------------------------------------------------------------
%%% Coast time = 1/2 Period
dt = pi*sqrt((at^3)/um); % sec

%%% Printing Answer
fprintf('------------------ 1d ------------------\n')
fprintf('Coast time = %f seconds\n\n',dt)

% ------------------------------------------------------------------------
%%% Problem 1e
% ------------------------------------------------------------------------
% CM mean motion
nCM = sqrt(um/(rat^3)); % rad/s
M = nCM*dt; % rads
e = 0; % eccentricity
E = M2E(M,e); % rads
ta = E2T(E,e); % rads
phase = 180 - ta*180/pi; % deg

%%% Printing Answer
fprintf('------------------ 1e ------------------\n')
fprintf('Phase angle = %f degrees\n\n',phase)

















