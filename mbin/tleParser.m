% i ----- inclination (rad)
% raan -- right ascension of node (rad)
% e ----- eccentricity
% w ----- argument of periapse (rad)
% M ----- mean anomaly (rad)
% n ----- mean motion (rad/sec)
% a ----- semimajor axis (km)
% year -- last two digits of epoch year
% day --- epoch day of year plus fraction
function [a, e, i, raan, w, M, n, year, day] = tleParser(tlePath)
tleFile = fopen(tlePath, 'r');
line1 = fgetl(tleFile);
line2 = fgetl(tleFile);

u = 398600.4418; % km^3 / s^2
deg2rad = pi/180;

i = str2num(line2(9:16))* deg2rad;
raan = str2num(line2(18:25))* deg2rad;
e = str2num(line2(27:33)) * 10^-7;
w = str2num(line2(35:42)) * deg2rad;
M = str2num(line2(44:51)) * deg2rad;
n = str2num(line2(53:63)) * 2*pi/(24*60*60); %rev/day to rad/sec
a = (u/(n^2))^(1/3);
year = str2num(line1(19:20));
day = str2num(line1(21:32));
    
end

