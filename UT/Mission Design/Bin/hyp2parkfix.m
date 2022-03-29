function [dva,dvae,ra,va,rppae,vppae,vinfpl,vPplhat,vaplhat,vpha,vphae,betaa,deltaa,Tpa,tsoia]...
    = hyp2parkfix(vPa,va,vinfahat,vinfamag,ra,dec,coepa,mua,rsoia)
% warning off
% Determines the hyperbolic arrival trajectory to a parking orbit of
% orbital elements a, e, and i
% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory
%------INPUTS-----
% rPa(3): Position of arrival planet
% vPa(3): Velocity of arrival planet
% va(3): Arrival velocity required in heliocentric frame
% vinfa(3): vinf vector
% vinfahat(3): vinf dir vector
% vinfamag: vinf magnitude
% rpl: equatorial radius of planet
% ra: right ascension of north pole of rotation
% dec: declination of north pole of rotation
% coepa(6): classical orbital elements of parking departure orbit 
% mua:  gravitational parameter of departure planet
% rsoi: radius of sphere of influence of planet
 
%------OUTPUTS-----
% dva(3): Delta-V required for hyperbolic depature trajectory planet frame
% dvae(3): Delta-V required for hyperbolic depature trajectory helio frame
% rppa(3): location at periapsis of parking orbit
% vppa(3): Velocity at periapsis of parking orbit 
% vinfpl(3): Hyperbolic depature V infinity
% vPplhat(3): velocity of planet in planetary frame
% vpha: Velocity of Hyperbolic trajectory at periapsis
% betaa: Angle between periapsis location and vinf asymptote of hyperbolic
%       orbit
% deltaa: Aiming radius of hyperbolic trajectory
% Tpa: Period of parking orbit [sec]
% tsoia: time on hyperbolic orbit between periapsis and rsoi intersection

% zero=1e-12;
ap = coepa(1);
ep = coepa(2);
ip = abs(coepa(3));
op = coepa(4);
wp = coepa(5);
tap = coepa(6);

vinfamag2=vinfamag*vinfamag;
ah=-mua/vinfamag2;
ip=mod(ip,360);
if ip>180 %adjust retrograde inclination to b/w 90-180 deg
    ip=360-ip;
end

% Transform v-inf from ecliptic to planetocentric equatorial
eps= 23.439291111; %(deg) Obliquity of ecliptic at J2000

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
vPahat=vPa/sqrt(vPa'*vPa);
% vPplhat=[(vPahat'*nhat); (vPahat'*yp); (vPahat'*npp)];
vPplhat=rot*vPahat;
vahat=va/sqrt(va'*va);
% vaplhat=[(vahat'*nhat); (vahat'*yp); (vahat'*npp)];
vaplhat=rot*vahat;
% vinfahat=[(vinfahat'*nhat); (vinfahat'*yp); (vinfahat'*npp)]; 
vinfahat=rot*vinfahat;
vinfpl=vinfahat*vinfamag;

    [ra, va]=orbel2rv(mua,ap,ep,ip,op,wp,tap);
    %parking orbit state at arrival in planet frame
    ramag = sqrt(ra'*ra);
    rahat = ra/ramag;
    
hhat3 = ra(1)*va(2) - ra(2)*va(1); %z-component

% Determine the departure velocity (same location) vector
D=sqrt(mua/(ramag*(1-rahat'*vinfahat))+vinfamag2/4);
vpha=(D+0.5*vinfamag)*vinfahat-(D-0.5*vinfamag)*rahat;
h3=ra(1)*vpha(2)-ra(2)*vpha(1); %z-component
if (sign(hhat3)~=sign(h3))
    D=-D;
    vpha=(D+0.5*vinfamag)*vinfahat-(D-0.5*vinfamag)*rahat;
end

% Determine hyperbolic trajectory data
    ev=vpha'*vpha*ra/mua-ra'*vpha*vpha/mua-rahat;
    eh = sqrt(ev'*ev);
    rph=ah*(1-eh);
    deltaa = rph*sqrt(1-2*ah/rph);
    betaa = acos(1/(1-rph/ah));
    Tpa=2*pi*(sqrt(ap^3/mua));

dva = va-vpha; %delta-V in planet frame
rot=rot'; %transformation from planetary frame to ecliptic frame
dvae=rot*dva; %delta-V in ecliptic frame
rppae=rot*ra; %hyperbolic periapse in ecliptic frame (planet centered)
vppae=rot*va; %parking orbit periapse velocity in ecliptic frame 
vphae=rot*vpha; %hyperbolic periapse velocity in ecliptic frame 


% Determine the time from periapsis of hyperbolic trajectory to rsoi
hhmag = sqrt(mua*ah*(1-eh*eh));
F = acosh((1-rsoia/ah)/eh);
Mh = eh*sinh(F)-F;
tsoia = hhmag^3*Mh/(mua*mua*(eh*eh-1)^(3/2));
if ep>1.0
Tpa=tsoia;    
end

end
