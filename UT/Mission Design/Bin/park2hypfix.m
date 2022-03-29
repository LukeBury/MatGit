function [dvd,dvde,rd,vd,rppde,vppde,vinfpl,vPplhat,vdplhat,vphd,vphde,betad,deltad,Tpd,tsoid]...
    = park2hypfix(vPd,vd,vinfdhat,vinfdmag,ra,dec,coepd,mud,rsoid)
% warning off
% Determines the hyperbolic departure trajectory from a parking orbit of
% orbital elements a, e, and i

% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory 

%------INPUTS-----
% vPd(3): Velocity of arrival planet
% vd(3): Arrival velocity required in heliocentric frame
% vinfdhat(3): vinf dir vector
% vinfdmag: vinf magnitude
% ra: right ascension of north pole of rotation
% dec: declination of north pole of rotation
% coepd(6): classical orbital elements of parking departure orbit 
% mud:  gravitational parameter of departure planet
% rsoi: radius of sphere of influence of planet
 
%------OUTPUTS-----
% dvd(3): Delta-V required for hyperbolic depature trajectory planet frame
% dvde(3): Delta-V required for hyperbolic depature trajectory helio frame
% rppd(3): location at periapsis of parking orbit
% vppd(3): Velocity at periapsis of parking orbit 
% vinfpl(3): Hyperbolic depature V infinity
% vPplhat(3): velocity of planet in planetary frame
% vphd: Velocity of Hyperbolic trajectory at periapsis
% betad: Angle between periapsis location and vinf asymptote of hyperbolic
%       orbit
% deltad: Aiming radius of hyperbolic trajectory
% Tpd: Period of parking orbit [sec]
% tsoid: time on hyperbolic orbit between periapsis and rsoi intersection


% zero=1e-12;
ap = coepd(1);
ep = coepd(2);
ip = abs(coepd(3));
op = coepd(4);
wp = coepd(5);
tap = coepd(6);

vinfdmag2=vinfdmag*vinfdmag;
ah=-mud/vinfdmag2;
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
vPdhat=vPd/sqrt(vPd'*vPd);
% vPplhat=[(vPdhat'*nhat); (vPdhat'*yp); (vPdhat'*npp)]; %dot
vPplhat=rot*vPdhat;%transform vector from ecliptic to planetary frame
vdhat=vd/sqrt(vd'*vd);
% vdplhat=[(vdhat'*nhat); (vdhat'*yp); (vdhat'*npp)];%dot
vdplhat=rot*vdhat;%transform vector from ecliptic to planetary frame
% vinfdhat=[(vinfdhat'*nhat); (vinfdhat'*yp); (vinfdhat'*npp)];%dot 
vinfdhat=rot*vinfdhat;%transform vector from ecliptic to planetary frames 
vinfpl=vinfdhat*vinfdmag; %v-inf in planet frame
  
    [rd, vd]=orbel2rv(mud,ap,ep,ip,op,wp,tap); 
    %parking orbit state at departure in planet frame
    rdmag = sqrt(rd'*rd);
    rdhat = rd/rdmag;

hhat3 = rd(1)*vd(2) - rd(2)*vd(1); %z-component
    
% Determine the departure velocity (same location) vector
D=sqrt(mud/(rdmag*(1+rdhat'*vinfdhat))+vinfdmag2/4);
vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rdhat;
h3=rd(1)*vphd(2)-rd(2)*vphd(1); %z-component
if (sign(hhat3)~=sign(h3))
    D=-D;
    vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rdhat;
end

% Determine hyperbolic trajectory data
    ev=vphd'*vphd*rd/mud-rd'*vphd*vphd/mud-rdhat;
    eh = sqrt(ev'*ev);
    rph=ah*(1-eh);
    deltad = rph*sqrt(1-2*ah/rph);
    betad = acos(1/(1-rph/ah));
    Tpd=2*pi*(sqrt(ap^3/mud));

dvd = vphd-vd; %delta-V in planet frame   
rot=rot'; %transformation from planetary frame to ecliptic frame
dvde=rot*dvd; %delta-V in ecliptic frame
rppde=rot*rd; %hyperbolic periapse in ecliptic frame (planet centered)
vppde=rot*vd; %parking orbit periapse velocity in ecliptic frame 
vphde=rot*vphd; %hyperbolic periapse velocity in ecliptic frame 

% Determine the time from periapsis of hyperbolic trajectory to rsoi
hhmag = sqrt(mud*ah*(1-eh*eh));
F = acosh((1-rsoid/ah)/eh); 
Mh = eh*sinh(F)-F;
tsoid = hhmag^3*Mh/(mud*mud*(eh*eh-1)^(3/2));
if ep>1.0
Tpd=tsoid;    
end

end

