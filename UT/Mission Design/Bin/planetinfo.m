function [mup,rpl,rap,dp,wp]=planetinfo(pn,JD,n)

%---inputs--------
% pn(n):    planet number, increasing w/ heliocentric radius (Mercury=1...
%           Neptune=8,ect) 
%           Note: Pluto=9, Sun=10
% JD(n):    Julian date
% n(1):     number of planets in vector

%---outputs-------
% MUp(n):      gravitational parameter
% rpl(n):       equatorial radius of planet  
% rap(n):      right ascension of planet's rotational axis
% dp(n):       declination of planet's rotational axis
% wp(n):       location of prime meridian of date along the planet's equator
%           with respect to acending node of the planet's equator on 
%           Earth's equatorial plane of J2000

%Note: Parameters retrieved from "The Astronomical Almanac" 2012 and
% "The Explanatory Supplement to the Astronomical Almanac" (1992)
Rj=71492;   %km
muj = 1.26712764e+008;      %km3/s2 

% all MU: km3/s2
mu(1) = 5959.916;  %Io
mu(2) = 3202.739;  %Europa
mu(3) = 9887.834;  %Ganymede
mu(4) = 7197.289;  %Callisto    
      

r(1)=1821.6; %Io all r: km
r(2)=1560.8; %Europa
r(3)=2631.2; %Ganymede
r(4)=2410.3; %Callisto

mup(1:n)=0; %initialize vector
rpl(1:n)=0; %initialize vector

for i=1:n
    if (pn(i)>=1)&&(pn(i)<=10)
mup(i)=mu(pn(i));
rpl(i)=r(pn(i));
    elseif pn(i)==0
mup(i)=muj;
rpl(i)=Rj;        
    end
end

rap(1:n)=0; %initialize vectors
dp(1:n)=90; 
wp(1:n)=0;

