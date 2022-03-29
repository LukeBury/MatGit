function [Year, Mon, Day, hr, min, s] = julian2calendar(JD)
%%% Description
%       Converts Julian date to Calendar date (algorithm from Vallado)
% --------------------------------------------------------------
%%% Inputs
%       JD - Julian Date
% --------------------------------------------------------------
%%% Outputs
%       yr  - year 
%       mo  - month
%       d   - day
%       h   - hour
%       min - minute
%       s   - seconds
% ===============================================================
T_1900 = (JD - 2415019.5)/365.25;
Year = 1900 + floor(T_1900);
LeapYrs = floor((Year - 1900 - 1)*0.25);
Days = (JD - 2415019.5) - ((Year - 1900)*(365) + LeapYrs);

if Days < 1.0
    Year = Year - 1;
    LeapYrs = floor((Year - 1900 - 1)*0.25);
    Days = (JD - 2415019.5) - ((Year - 1900)*365 + LeapYrs);
end

LMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
if mod(Year,4) == 0
    LMonth(2) = 29;
end

DayOfYr = floor(Days);

daysInMonthSummation = 0;
Mon = 0;
while (daysInMonthSummation + 1) <= DayOfYr
    Mon = Mon + 1;
    daysInMonthSummation = daysInMonthSummation + LMonth(Mon);
end

Day = DayOfYr - daysInMonthSummation + LMonth(Mon);

Tau = (Days - DayOfYr)*24;
hr = floor(Tau);
min = floor((Tau - hr)*60);
s = (Tau - hr - min/60)*3600;


end
	   
