function [DV1, DV2, DVTot, TpTrans_sec] = interplanetaryHohmannTransfer_DVs(centralBody_mu, body1_mu, body2_mu, body1_R, body2_R, RSc_body1, RSc_body2, VSc_body1, VSc_body2)
%%% Description
%       Computes parameters for interplanetary hohmann transfer
% --------------------------------------------------------------
%%% Inputs
%       centralBody_mu - [1x1] gravitational parameter of central body (km^2/s^3) 
%       body1_mu       - [1x1] gravitational parameter of first body (starting location) (km^2/s^3)
%       body2_mu       - [1x1] gravitational parameter of second body (final location) (km^2/s^3)
%       body1_R        - [1x1] semi-major axis of body1 about central body (km)
%       body2_R        - [1x1] semi-major axis of body2 about central body (km)
%       RSc_body1      - [1x1] Distance of spacecraft from body1 in initial orbit (km)
%       RSc_body2      - [1x1] Distance of spacecraft from body2 in final orbit (km)
%       VSc_body1      - [1x1] Velocity of spacecraft with respect to body1 in initial orbit (km/s)
%       VSc_body2      - [1x1] Velocity of spacecraft with respect to body2 in final orbit (km/s)
% --------------------------------------------------------------
%%% Outputs
%       DV1         - [1x1] Magnitude of first delta-V which is required to initiate the transfer (km/s)
%       DV2         - [1x1] Magnitude of second delta-V which completes the transfer (km/s)
%       DVTot       - [1x1] Total delta-V (km/s)
%       TpTrans_sec - [1x1] Time required for hohmann transfer (s)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Heliocentric velocities of each body
% ------------------------------------
V_body1 = sqrt(centralBody_mu / body1_R); % km/s
V_body2 = sqrt(centralBody_mu / body2_R); % km/s

% ------------------------------------
%%% Use VisViva to look at heliocentric spacecraft velocities at rp and ra
%%% of transfer orbit
% ------------------------------------
%%% rp, ra, and semimajor axis of transfer orbit
if     body1_R > body2_R
    rpTrans = body2_R;                % km
    raTrans = body1_R;                % km
elseif body1_R < body2_R
    rpTrans = body1_R;                % km
    raTrans = body2_R;                % km
end

aTrans = (rpTrans + raTrans) / 2; % km

%%% heliocentric Sc velocities at rp and ra of transfer orbit
VSc_rpTrans = visviva_v(rpTrans, aTrans, centralBody_mu); % km/s
VSc_raTrans = visviva_v(raTrans, aTrans, centralBody_mu); % km/s

%%% Transfer time
TpTrans_sec = pi*sqrt((aTrans^3)/centralBody_mu); % sec
% ------------------------------------
%%% V-infinities of Sc at both ends of the transfer
% ------------------------------------
if     body1_R > body2_R
    VInf_body1 = V_body1 - vSc_raTrans; % km/s
    VInf_body2 = VSc_rpTrans - V_body2; % km/s
elseif body1_R < body2_R
    VInf_body1 = VSc_rpTrans - V_body1; % km/s
    VInf_body2 = V_body2 - VSc_raTrans; % km/s
end

% ------------------------------------
%%% Sc velocities on either side of each burn (Plus & Minus)
% ------------------------------------
%%% Burn 1
VSc_DV1_m = VSc_body1;                                  % km/s
VSc_DV1_p =  sqrt(VInf_body1^2 + 2*body1_mu/RSc_body1); % km/s

%%% Burn 2
VSc_DV2_m = sqrt(VInf_body2^2 + 2*body2_mu/RSc_body2); % km/s
VSc_DV2_p = VSc_body2;                                 % km/s

% ------------------------------------
%%% Calculate DV values
% ------------------------------------
DV1 = abs(VSc_DV1_p - VSc_DV1_m); % km/s
DV2 = abs(VSc_DV2_m - VSc_DV2_p); % km/s
DVTot = DV1 + DV2;           % km/s
TpTrans_sec
warning('Could use a quick looking over')
end