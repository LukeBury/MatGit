function [r_heliocentric, v_heliocentric] = JulianDate_2_HeliocentricState(JD, body)
%%% Description
%       Provide a julian date and body name, receive heliocentric state
%       
% --------------------------------------------------------------
%%% Inputs
%       JD - Julian date (days) [1x1]
%       body - string of capitalized body name... ex: 'Earth'
% --------------------------------------------------------------
%%% Outputs
%       r_heliocentric - heliocentric position of body (km) [3x1]
%       v_heliocentric - heliocentric velocity of body (km/s) [3x1]
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Using JD to get orbital elements
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD, body, 'radians');

%%% Gravitational parameter of sun
mu_sun = 1.32712440018e11; % km^3/s^2

%%% Acquire heliocentric states
[r_heliocentric, v_heliocentric] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, mu_sun);

end