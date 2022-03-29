function [dvd]= park2hypdOPT(vinfdhat,vinfdmag,ra,dec,coepd,mud)
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
% deg2rad=pi/180;
ap = coepd(1);
ep = coepd(2);
ip = abs(coepd(3));
% op = coepd(4);
% wp = coepd(5);
% tap = 0; %true anomaly; set to zero b/c departure at periapsis
ip=mod(ip,360);
if ip>180 %adjust retrograde inclination to b/w 90-180 deg
    ip=360-ip;
end
ip=ip*pi/180;
rppdmag = ap*(1-ep);
vppdmag=sqrt(mud*(2/rppdmag-1/ap));
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
decinf=asin(vinfdhat(3)); %declination of v-inf vector planet frame
if (abs(decinf)<=ip)&&(ip<=pi-abs(decinf))  % Check if capable of coplanar periapse transfer
%     rppdmag = ap*(1-ep); %periapsis radius (km)
%     rapdmag = ap*(1+ep); %apoapsis radius (km)
    vcmag=(mud/rppdmag); %velocity of circular orbit with same periapsis 
%     dvd=sqrt(2*vcmag+vinfdmag2)-sqrt(2*vcmag*rapdmag/(rppdmag+rapdmag)); %scalar delta-V 
    dvd=abs(sqrt(2*vcmag+vinfdmag*vinfdmag)-vppdmag); %scalar delta-V 

else
    decinfp=pi/2-decinf; %ajusted declination from pole
    hhat=-sign(decinf)*sin(ip)/sin(decinfp)*vinfdhat+[0; 0; sin(ip+decinfp)/sin(decinfp)];
    if ip<=pi/2
        node=[-hhat(2); hhat(1); 0]/sin(ip);
    else
        node=[hhat(2); -hhat(1); 0]/sin(ip);
    end

if ip==0 % If Prograde Equatorial: 
        % set w = 0, RAAN = 0, and TA (or w) = lamdaTrue 
    op=0;
    wpi = rainf-acos(-1/(1-rppdmag/ah))*cos(decinf); %initial guess (rad)
%      wp =rainf-acos(-1/(1-rppdmag/ah))*cos(decinfp) %initial guess (rad)
    [~,dvd]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,rppdmag,vppdmag,ip,op),wpi,...
        optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
elseif ip==pi % If retrograde Equatorial:
             % set RAAN = 0 and w = whatTrue 
    op=0;
    wpi = -rainf-acos(-1/(1-rppdmag/ah))*cos(decinf); %initial guess (rad)
%     wp =-rainf-acos(-1/(1-rppdmag/ah))*cos(decinfp) %initial guess (rad)
    [~,dvd]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,rppdmag,vppdmag,ip,op),wpi,...
        optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
else
    op=atan2(node(2),node(1));
%     ip+decinfp,asin((1/(1-rppdmag/ah)))
%     wpi=pi-asin((1/(1-rppdmag/ah))/sin(ip+decinfp))
%     wpi=pi-asin((1/(1-rppdmag/ah))*cos(ip+decinfp))%initial guess (rad)
    wpi=-sign(decinf)*pi/2+pi/4;
    [~,dvd]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,rppdmag,vppdmag,ip,op),wpi,...
        optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));

end
end
    
function f = dvmin(x,vinfdhat,vinfdmag,mud,rppmag,vppmag,ip,op)

%     [rpp, vpp] = orbel2rvRAD(mud,ap,ep,ip,op,x,0);
%---------------------------------------------------------------------
% rmag = a*(1-e^2)/(1+e*cos(TA));
% vmag = sqrt(mu*(2/rmag-1/a));
% FPA = atan2(e*sin(TA),1+e*cos(TA)); %b/c TA & FPA = 0
% taw=x; %TA+w;
coso=cos(op);
sino=sin(op);
cosi=cos(ip);
sini=sin(ip);
costaw=cos(x);
sintaw=sin(x);
    rpp = rppmag*[(costaw*coso - sintaw*cosi*sino);
        (costaw*sino + sintaw*cosi*coso);
        (sintaw*sini)];
%     taw=taw; %taw-FPA=0
%     costaw=cos(taw);
%     sintaw=sin(taw);
    vpp = vppmag*[(-sintaw*coso - costaw*cosi*sino);
        (-sintaw*sino + costaw*cosi*coso);
        (costaw*sini)];
%---------------------------------------------------------------------    rppmag = sqrt(rpp'*rpp);
    rpphat = rpp/rppmag;
    
% Determine the departure velocity (same location) vector
D=sqrt(mud/(rppmag*(1+rpphat'*vinfdhat))+vinfdmag*vinfdmag/4);
vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rpphat;
    
f=sqrt((vphd-vpp)'*(vphd-vpp));
end
end

