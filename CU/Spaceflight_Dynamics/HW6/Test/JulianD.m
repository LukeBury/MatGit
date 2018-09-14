function [Jd]=JulianD(mo, day, yr, hr, min)

Jd=367*yr-floor((7*((yr)+floor((mo+9)/12)))/4)+floor(275*mo/9)+day+1721013.5+(hr+min/60)/24;