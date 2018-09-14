%%% Accepts year, month, day, hour, and minute (w/ fraction) of time
%%% Converts into Julian date
function [Jd]=JD(yr, mo, day, hr, min)

Jd=367*yr - floor((7*(yr+floor((mo+9)/12)))/4) + floor(275*mo/9) + day +...
    1721013.5 + (hr+min/60)/24;
end