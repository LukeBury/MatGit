function [Rtopo] = ecef2topo(Recef,lat,lon,alt)
Re=6378.137; %radius of the Earth
r=Re+alt; %magnitude position of the topocentric frame with respect
%to the center of the earth
rsite=[r*cos(lat)*cos(lon);r*cos(lat)*sin(lon);r*sin(lat)];
%vector form of the postion of the topocentric frame with resepct ot the
%center of the Earth
pecef=Recef-rsite; %range from satellite to topocentric frame center
range=norm(pecef);
R2nlat=[cos(pi/2-lat),0,-sin(pi/2-lat);0,1,0;sin(pi/2-lat),0,cos(pi/2-lat)];
%1st rotation of the axis about the ECEF 
R3lon=[cos(lon),sin(lon),0;-sin(lon),cos(lon),0;0,0,1];
%2nd rotation of the axis about the new transfer frame
Rtopo=R2nlat*R3lon*pecef; %final product of adding both rotations together
%to get topocentric frame

