function [dva]= hypa2parkOPT(vinfahat,vinfamag,ra,dec,coepa,mua)
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
% op = coepa(4);
% wp = coepa(5);
% tap = 0; %true anomaly; set to zero b/c arrival at periapsis
deg2rad=pi/180;
ip=mod(ip,360);
if ip>180 %adjust retrograde inclination to b/w 90-180 deg
    ip=360-ip;
end
ip=ip*deg2rad;
rppamag = ap*(1-ep);
vppamag=sqrt(mua*(2/rppamag-1/ap));

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
decinf=asin(vinfahat(3)); %declination of v-inf vector planet frame

if (abs(decinf)<=ip)&&(ip<=pi-abs(decinf)) % Check if capable of coplanar periapse transfer
    vcmag=(mua/rppamag); %velocity of circular orbit with same periapsis 
%     dva=sqrt(2*vcmag+vinfamag*vinfamag)-sqrt(2*vcmag*rapamag/(rppamag+rapamag)); %scalar delta-V 
    dva=abs(sqrt(2*vcmag+vinfamag*vinfamag)-vppamag); %scalar delta-V 

else
    ah=-mua/vinfamag*vinfamag;
    rainf=atan2(vinfahat(2),vinfahat(1)); %right ascension of v-inf vector planet frame
    decinfp=pi/2-decinf; %ajusted declination from pole
    hhat=-sign(decinf)*sin(ip)/sin(decinfp)*vinfahat+[0; 0; sin(ip+decinfp)/sin(decinfp)];
    if ip<=pi/2
        node=[-hhat(2); hhat(1); 0]/sin(ip);
    else
        node=[hhat(2); -hhat(1); 0]/sin(ip);
    end
    
    if ip==0 % If Prograde Equatorial: 
            % set w = 0, RAAN = 0, and TA (or w) = lamdaTrue 
        wpi = (pi+rainf)+acos(-1/(1-rppamag/ah))*cos(decinf); %initial guess (rad)
%         wpi =rainf-acos(-1/(1+ap*(1-ep)*vinfamag*vinfamag/mua))*cos(decinfp); %initial guess (rad)
        [~,dva]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,rppamag,vppamag,ip,0),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));

    elseif ip==pi % If retrograde Equatorial:
                 % set RAAN = 0 and w = whatTrue 
        wpi = -(pi+rainf)+acos(-1/(1-rppamag/ah))*cos(decinf); %initial guess (rad)         
%         wpi =-rainf-acos(-1/(1+ap*(1-ep)*vinfamag*vinfamag/mua))*cos(decinfp); %initial guess (rad)
        [~,dva]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,rppamag,vppamag,ip,0),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));

    else
        op=atan2(node(2),node(1));
        wpi=sign(decinf)*pi/2-pi/4;
%         wpi=-asin(-(1/(1+ap*(1-ep)*vinfamag*vinfamag/mua))/sin(ip+decinfp)); %initial guess (rad)
        [~,dva]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,rppamag,vppamag,ip,op),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    end
end
end
function f = dvmin(x,vinfahat,vinfamag,mua,rppmag,vppmag,ip,op)

%     [rpp, vpp] = orbel2rvRAD(mua,ap,ep,ip,op,x,0);
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
%---------------------------------------------------------------------

%     rppmag = sqrt(rpp'*rpp);
    rpphat = rpp/rppmag;
    
% Determine the departure velocity (same location) vector
D=sqrt(mua/(rppmag*(1-rpphat'*vinfahat))+vinfamag*vinfamag/4);
vpha=(D+0.5*vinfamag)*vinfahat-(D-0.5*vinfamag)*rpphat;   
    
f=sqrt((vpha-vpp)'*(vpha-vpp));
end