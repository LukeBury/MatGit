function [dvd,dvde,rppd,vppd,rppde,vppde,vinfpl,vPplhat,vdplhat,vphd,vphde,betad,deltad,Tpd,tsoid,coepd]...
    = park2hypopt(vPd,vd,vinfdhat,vinfdmag,ra,dec,coepd,mud,rsoid)
% warning off
% Determines the hyperbolic departure trajectory from a parking orbit of
% orbital elements a, e, and i

% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory 

%------INPUTS-----
% rPd(3): Position of arrival planet
% vPd(3): Velocity of arrival planet
% vd(3): Arrival velocity required in heliocentric frame
% vinfd(3): vinf vector
% vinfdhat(3): vinf dir vector
% vinfdmag: vinf magnitude
% rpl: equatorial radius of planet
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
% coepd(6): Classical orbital elements of parking orbit (Determined for minimum
%       arrival delta-V

% zero=1e-12;
ap = coepd(1);
ep = coepd(2);
ip = abs(coepd(3));
% op = coepd(4);
% wp = coepd(5);
% tap = 0; %true anomaly; set to zero b/c departure at periapsis
deg2rad=pi/180;
vinfdmag2=vinfdmag*vinfdmag;
rppdmag = ap*(1-ep); %periapsis radius (km)
ip=mod(ip,360);
if ip>180 %adjust retrograde inclination to b/w 90-180 deg
    ip=360-ip;
end
ip=ip*deg2rad;
ah=-mud/vinfdmag2;
deltad = rppdmag*sqrt(1+2*mud/(rppdmag*vinfdmag2));
betad = acos(1/(1+rppdmag*vinfdmag2/mud));
Tpd=2*pi*(sqrt(ap*ap*ap/mud));

% Transform v-inf from ecliptic to planetocentric equatorial 

ecl= 0; % Obliquity of ecliptic at J2000

npp=[cosd(ra)*cosd(dec);
    cosd(ecl)*sind(ra)*cosd(dec)+sind(ecl)*sind(dec); 
    -sind(ecl)*sind(ra)*cosd(dec)+cosd(ecl)*sind(dec)];
npe=[0; sind(ecl); cosd(ecl)];
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
vcmag=sqrt(mud/rppdmag); %velocity of circular orbit with same periapsis 

dvmag=abs(sqrt(2*vcmag*vcmag+vinfdmag2)-vcmag*sqrt(1+ep)); %scalar delta-V 
decinf=asin(vinfdhat(3)); %declination of v-inf vector planet frame
rainf=atan2(vinfdhat(2),vinfdhat(1)); %right ascension of v-inf vector planet frame
if (abs(decinf)<=ip)&&(ip<=pi-abs(decinf))  % Check if capable of coplanar periapse transfer
   if ip==0
        alpha=pi+rainf; % or  pi-rainf-B
        op=0; %RAAN undefined
        wp =rainf-acos(-1/(1-rppdmag/ah)); % argument of periapse
    elseif ip==pi    
        alpha=pi+rainf; % or  pi-rainf-B
        op=0; %RAAN undefined
        wp =-rainf-acos(-1/(1-rppdmag/ah)); % argument of periapse
    else
        alpha=pi+rainf-acos(tan(decinf)*cot(ip)); % or pi-rainf-B  
        op = alpha+pi/2; %detailed eqn below
%         op=pi+rainf+asin(cot(pi/2-decinf)/tan(ip))
        wp=acos(sin(decinf)/sin(ip))-asin(1/(1-rppdmag/ah));
%       RAAN=360+rainf*180/pi-asind(cot(B)/tand(ip)) %Other solution
%       W=-acosd(cos(B)/sind(ip))+asind(1/(1+rppdmag*vinfdmag2/mud)) %Other solution
    end
else
    alpha=pi+rainf; 
    ip=decinf;
    op = alpha+pi/2; %detailed eqn below
%       op=pi+rainf+asin(cot(decinfp)/tan(ip))
    wp=-asin(1/(1-rppdmag/ah));
end

hhat=[cos(alpha)*sin(ip); sin(alpha)*sin(ip); cos(ip)]; 
%hyperbolic trajectory angular momentum unit vector
    rppdhat= [vinfdhat(2)*hhat(3) - vinfdhat(3)*hhat(2);...%cross product 
            vinfdhat(3)*hhat(1) - vinfdhat(1)*hhat(3);...
            vinfdhat(1)*hhat(2) - vinfdhat(2)*hhat(1)];
    rppdhat=cos(pi-betad)*vinfdhat+sin(pi-betad)*rppdhat; 
%location unit vector at periapsis of parking orbit
rppd=rppdmag*rppdhat; %location vector at periapsis of parking orbit
dvd = dvmag*[hhat(2)*rppdhat(3) - hhat(3)*rppdhat(2);...%cross product 
            hhat(3)*rppdhat(1) - hhat(1)*rppdhat(3);...
            hhat(1)*rppdhat(2) - hhat(2)*rppdhat(1)];
D=sqrt(mud/(rppdmag*(1+rppdhat'*vinfdhat))+vinfdmag2/4);
vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rppdhat; 
%hyperbolic velocity at periapsis
vppd=vphd-dvd; %parking orbit velocity at periapsis

% [ap,ep,ip,op,wp,~,~] = rv2orbel(rppd,vppd,mud); %find orbital elements
%Find remaining orbital elements (op and wp)
% evec=rppdhat*ep;
% hh=rppdmag*sqrt(vppd'*vppd)*hhat;
% n=[-hh(2);hh(1);0];
% nmag=sqrt(n'*n);
% if ep>zero
%     if (abs(ip)<zero)||(180-abs(ip)<zero) 
%         % If Elliptical Equatorial:
%         op=0;         % set o = 0 and w = whatTrue 
%         wp = acosd(evec(1)/ep);
%         if ((evec(2)<0)&&(abs(ip)<zero))||((-evec(2)<0)&&(180-abs(ip)<zero))
%             wp = 360 -wp;
%         end
%     else            % Normal cases
%         op = acosd(n(1)/nmag);
%         if n(2)<0
%             op = 360-op;
%         end
%         wp = acosd((n'*evec)/(nmag*ep));
%         if evec(3)<0
%             wp = 360-wp;
%         end
%     end 
% else
%     if ((abs(ip)>zero))&&(180-abs(ip)>zero)    % If Circular Inclined: 
%         %replaced u/ta w/ wp to determine orbit periapse location
%         if (n'*rppd)/(nmag*rppdmag)>1.0
%             wp=0; %u
%         else
%             wp = acosd((n'*rppd)/(nmag*rppdmag)); %u
%         end
%         if rppd(3)<0
%             wp = 360-wp;
%         end
%         op = acosd(n(1)/nmag);
%         if n(2)<0
%             op = 360-op;
%         end
%     else
%         % If Circular Equatorial: 
%         op=0;
%         wp = acosd(rppdhat(1)); %lamdaTrue
% 
%         if (180-abs(ip)<zero) 
%             wp = -wp;  %lamdaTrue
%         end
%         
%         if (rppd(2)<0)
%             wp = 360-wp;  %lamdaTrue 
%         end
%     end
% end
coepd = [ap ep ip/deg2rad op/deg2rad wp/deg2rad 0]; %orbital elements in planet frame

rot=rot'; %transformation from planetary frame to ecliptic frame
dvde=rot*dvd; %delta-V in ecliptic frame
rppde=rot*rppd; %hyperbolic periapse in ecliptic frame (planet centered)
vppde=rot*vppd; %parking orbit periapse velocity in ecliptic frame 
vphde=rot*vphd; %hyperbolic periapse velocity in ecliptic frame 

% Determine the time from periapsis of hyperbolic trajectory to rsoi
eh = 1-rppdmag/ah;
hhmag = rppdmag*sqrt(vinfdmag2+2*mud/rppdmag);
F = acosh((1-rsoid/ah)/eh);
Mh = eh*sinh(F)-F;
tsoid = hhmag*hhmag*hhmag*Mh/(mud*mud*(eh*eh-1)^(3/2));
if ep>1.0
Tpd=tsoid;    
end

end