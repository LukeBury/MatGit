function [ JD ] = calendar2julian(yr, mo, d, h, min, s)
%%% Description
%       Converts Calendar date to Julian date 
% --------------------------------------------------------------
%%% Inputs
%       yr  - year 
%       mo  - month
%       d   - day
%       h   - hour
%       min - minute
%       s   - seconds
% --------------------------------------------------------------
%%% Outputs
%       JD - Julian Date
% ===============================================================
%%% Correcting month/year for January and February
if mo == 1 || mo == 2
    yr = yr - 1;
    mo = mo + 12;
end

%%% Calculating JD
B = 2 - floor(yr/100) + floor(floor(yr/100)/4);

C = (((s/60) + min)/60 + h) / 24;

JD = floor(365.25*(yr + 4716)) + floor(30.6001*(mo+1)) + d + B - 1524.5 + C;

%%% This one also works
% JD = 367*yr - floor((7*(yr + floor((mo + 9)/12)))/4) + floor(275*mo/9) + d + 1721013.5 + ((((s/60) + min)/60 + h) / 24);

warning('Not equipped for leap seconds')
end









