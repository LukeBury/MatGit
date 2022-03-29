function [dva]= hypa2vinfOPT(vinfahat,vinfamag,rpl,ra,dec,vinfavmag,ip,altp,mua,dvopt)
% warning off
% Determines the hyperbolic arrival trajectory to a parking orbit of
% orbital elements a, e, and i

%------INPUTS-----
% vPa(3): Velocity of arrival planet
% va(3): Arrival velocity required in heliocentric frame
% vinfahat(3): vinf dir vector
% vinfamag: vinf magnitude
% rpl: equatorial radius of planet
% ra: right ascension of north pole of rotation
% dec: declination of north pole of rotation
% vinfavmag: vinf available from launch
% ip: launch inclination
% altp: altitude of launch periapsis
% mud:  gravitational parameter of departure planet
% rsoi: radius of sphere of influence of planet
% dvopt: optimization of dv at launch periapsis
 
%------OUTPUTS-----
% dva: Delta-V required for hyperbolic depature trajectory planet frame

% zero=1e-12;
deg2rad=pi/180;
ip=abs(ip)*deg2rad;
rppamag = rpl+altp;
vppamag=sqrt(2*mua/rppamag+vinfavmag*vinfavmag);
vinfamag2=vinfamag*vinfamag;
% tap = 0; %true anomaly; set to zero b/c arrival at periapsis

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
vinfahat=rot*vinfahat;
decinf=asin(vinfahat(3)); %declination of v-inf vector planet frame
rainf=atan2(vinfahat(2),vinfahat(1)); %right ascension of v-inf vector planet frame
if (abs(decinf)<=ip)&&(ip<=pi-abs(decinf))  % Check if capable of coplanar periapse transfer
   vphamag = sqrt(2*mua/rppamag+vinfamag*vinfamag); 
    if dvopt == 1
            if vphamag<vppamag
                dva=0; 
            else
                dva=abs(vphamag-vppamag); %delta-V in planet frame  
            end
    else
        dva=abs(vphamag-vppamag); %delta-V in planet frame 
    end
else
    ah= -mua/vinfamag2;
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
%         wpi =rainf-acos(-1/(1+rppamag*vinfamag*vinfamag/mua))*cos(decinfp); %initial guess (rad)
        [~,dva]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,rppamag,vppamag,ip,0,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    
    elseif ip==pi % If retrograde Equatorial:
                 % set RAAN = 0 and w = whatTrue 
        wpi = -(pi+rainf)+acos(-1/(1-rppamag/ah))*cos(decinf); %initial guess (rad)                          
%         wpi =-rainf-acos(-1/(1+rppamag*vinfamag*vinfamag/mua))*cos(decinfp); %initial guess (rad)
        [~,dva]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,rppamag,vppamag,ip,0,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
   
    else
        op=atan2(node(2),node(1));
%         wpi=-asin(-(1/(1+rppamag*vinfamag*vinfamag/mua))/sin(ip+decinfp)); %initial guess (rad)
        wpi=sign(decinf)*pi/2-pi/4;
        [~,dva]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,rppamag,vppamag,ip,op,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
  
    end

end

end
function f = dvmin(x,vinfahat,vinfamag,mua,rppmag,vppmag,ip,op,dvopt)

%     [rpp, vpp] = orbel2rvRAD(mua,ap,ep,ip,op,x,0);
%---------------------------------------------------------------------
% rmag = a*(1-e^2)/(1+e*cos(TA));
% vmag = sqrt(mu*(2/rmag-1/a));
% FPA = atan2(e*sin(TA),1+e*cos(TA)); %b/c TA & FPA = 0
% taw=TA+w;
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
vphamag=sqrt(vpha'*vpha);
if dvopt == 1
    if vphamag>vppmag
        f=sqrt((vpha-vpp)'*(vpha-vpp));  %basic DV to match entry velocity
    else 
        costheta=(vpha'*vpp/(vphamag*vppmag));
        vphmag=vphamag/costheta; %can potentially be undefined if cos(theta)=0
        if vphmag<vppmag %check if perpendicular DV aligns entry but < vppmag
            f=sqrt(abs(vphmag*vphmag-vphamag*vphamag)); %perpendicular DV
        else
            f=sqrt((vpha-vpp)'*(vpha-vpp));  %basic DV to match entry velocity
        end
    end
else
    f=sqrt((vpha-vpp)'*(vpha-vpp));  
end    
end