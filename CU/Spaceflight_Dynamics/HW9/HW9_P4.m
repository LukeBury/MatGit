clear
clc

% ------------------------------------------------------------------------
%%% Setup
% ------------------------------------------------------------------------
%%% Mass
m_Sun = 1.9891e30; % kg, Sun
m_J = 1.8988e27; % kg, Jupiter
m_S = 5.685e26; % kg, Saturn

%%% Semimajor Axis of Circular Orbit
a_J = 778298361; % km, Jupiter
a_S = 1429394133; % km, Saturn

%%% Satellite distance
a_JSat = 29540000; % km, Jupiter satellite
a_SSat = 24504879; % km, Saturn satellite

% ------------------------------------------------------------------------
%%% Jupiter SOI 
% ------------------------------------------------------------------------
rSOI_J = ((m_J/m_Sun)^(2/5))*a_J % km, Jupiter SOI
SOI_67p = .67*rSOI_J % km, 62% of SOI
a_JSat/SOI_67p
fprintf('Stable, but barely\n')

% ------------------------------------------------------------------------
%%% Saturn SOI 
% ------------------------------------------------------------------------
rSOI_S = ((m_S/m_Sun)^(2/5))*a_S % km, Saturn SOI
SOI_67p = .67*rSOI_S % km, 62% of SOI
a_SSat/SOI_67p
fprintf('Stable, but barely, but less barely\n')

%Question, in regard to the last sentence of #4, does he mean "is it
%possible more satellites could be discovered beyond this one?"
