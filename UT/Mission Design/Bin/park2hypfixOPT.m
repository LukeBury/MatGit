function [dvd]= park2hypfixOPT(vinfdhat,vinfdmag,ra,dec,coepd,mud)
% warning off
% Determines the hyperbolic departure trajectory from a parking orbit of
% orbital elements a, e, and i

% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory 

%------INPUTS-----
% vinfdhat(3): vinf dir vector
% vinfdmag: vinf magnitude
% ra: right ascension of north pole of rotation
% dec: declination of north pole of rotation
% coepd(6): classical orbital elements of parking departure orbit 
% mud:  gravitational parameter of departure planet
 
%------OUTPUTS-----
% dvd: Delta-V required for hyperbolic depature trajectory planet frame

% zero=1e-12;
ap = coepd(1);
ep = coepd(2);
ip = abs(coepd(3));
op = coepd(4);
wp = coepd(5);
tap = coepd(6);

ip=mod(ip,360);
if ip>180 %adjust retrograde inclination to b/w 90-180 deg
    ip=360-ip;
end

% Transform v-inf from ecliptic to planetocentric equatorial
eps= 0; % Obliquity of ecliptic at J2000
npp=[cosd(ra)*cosd(dec);
    cosd(eps)*sind(ra)*cosd(dec)+sind(eps)*sind(dec); 
    -sind(eps)*sind(ra)*cosd(dec)+cosd(eps)*sind(dec)];
npe=[0; sind(eps); cosd(eps)];
nhat = [npe(2)*npp(3) - npe(3)*npp(2);...%cross product 
        npe(3)*npp(1) - npe(1)*npp(3);...%line of node of planet equator
        npe(1)*npp(2) - npe(2)*npp(1)];
nhat=nhat/sqrt(nhat'*nhat);
yp = [npp(2)*nhat(3) - npp(3)*nhat(2);...%cross product 
      npp(3)*nhat(1) - npp(1)*nhat(3);...%y-axis of planet equatorial frame in ecliptic frame
      npp(1)*nhat(2) - npp(2)*nhat(1)];
rot=[nhat yp npp]'; %transformation from ecliptic to planetary frame
vinfdhat=rot*vinfdhat;%transform vector from ecliptic to planetary frames 

    [rd, vd]=orbel2rv(mud,ap,ep,ip,op,wp,tap); 
    %parking orbit state at departure in planet frame
    rdmag = sqrt(rd'*rd);
    rdhat = rd/rdmag;
hhat3 = rd(1)*vd(2) - rd(2)*vd(1); %z-component
    
% Determine the departure velocity (same location) vector
D=sqrt(mud/(rdmag*(1+rdhat'*vinfdhat))+vinfdmag^2/4);
vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rdhat;
h3=rd(1)*vphd(2)-rd(2)*vphd(1);
if (sign(hhat3)~=sign(h3))
    D=-D;
    vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rdhat;
end
dvd=vphd-vd;
dvd=sqrt(dvd'*dvd);
end