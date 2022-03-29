function [dva]= hyp2parkfixOPT(vinfahat,vinfamag,ra,dec,coepa,mua)
% warning off
% Determines the hyperbolic arrival trajectory to a parking orbit of
% orbital elements a, e, and i

% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory

%------INPUTS-----
% vinfahat(3): vinf dir vector
% vinfamag: vinf magnitude
% ra: right ascension of north pole of rotation
% dec: declination of north pole of rotation
% coepa(6): classical orbital elements of parking departure orbit 
% mua:  gravitational parameter of departure planet
 
%------OUTPUTS-----
% dva: Delta-V required for hyperbolic depature trajectory planet frame

% zero=1e-12;
ap = coepa(1);
ep = coepa(2);
ip = abs(coepa(3));
op = coepa(4);
wp = coepa(5);
tap = coepa(6);

ip=mod(ip,360);
if ip>180 %adjust retrograde inclination to b/w 90-180 deg
    ip=360-ip;
end

% Transform v-inf from ecliptic to planetocentric equatorial
eps= 0; % (deg) Obliquity of ecliptic at J2000
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
vinfahat=rot*vinfahat;

[ra, va]=orbel2rv(mua,ap,ep,ip,op,wp,tap);
    %parking orbit state at arrival in planet frame
    ramag = sqrt(ra'*ra);
    rahat = ra/ramag;
hhat3 = ra(1)*va(2) - ra(2)*va(1); %z-component

% Determine the departure velocity (same location) vector
D=sqrt(mua/(ramag*(1-rahat'*vinfahat))+vinfamag^2/4);
vpha=(D+0.5*vinfamag)*vinfahat-(D-0.5*vinfamag)*rahat;
h3=ra(1)*vpha(2)-ra(2)*vpha(1); %z-component
if (sign(hhat3)~=sign(h3))
    D=-D;
    vpha=(D+0.5*vinfamag)*vinfahat-(D-0.5*vinfamag)*rahat;
end
dva = va-vpha; %delta-V in planet frame
dva=sqrt(dva'*dva);
end


