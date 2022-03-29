function [dvd]= vinf2hypdOPT(vinfdhat,vinfdmag,rpl,ra,dec,vinfavmag,ip,altp,mud,dvopt)
% warning off
% Determines the hyperbolic departure trajectory from a parking orbit of
% orbital elements a, e, and i

% NOTE: The periapsis location of parking orbit is configured as periapsis 
% location of hyperbolic departure trajectory 

%------INPUTS-----
% vinfdhat(3): vinf dir vector
% vinfdmag: vinf magnitude
% rpl: equatorial radius of planet
% ra: right ascension of north pole of rotation
% dec: declination of north pole of rotation
% vinfavmag: vinf available from launch
% ip: launch inclination
% altp: altitude of launch periapsis 
% mud: gravitational parameter 
% dvopt: optimization of dv at launch periapsis
 
%------OUTPUTS-----
% dvd: Delta-V required for hyperbolic depature trajectory planet frame


% zero=1e-12;
%Determine initial launch trajectory parameters
% deg2rad=pi/180;
ip=abs(ip)*pi/180;
rppdmag = rpl+altp;
vppdmag=sqrt(2*mud/rppdmag+vinfavmag*vinfavmag);
vinfdmag2=vinfdmag*vinfdmag;

% tap = 0; %true anomaly; set to zero b/c departure at periapsis

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
vinfdhat=rot*vinfdhat;%transform vector from ecliptic to planetary frames
decinf=asin(vinfdhat(3)); %declination of v-inf vector planet frame
rainf=atan2(vinfdhat(2),vinfdhat(1)); %right ascension of v-inf vector planet frame
if (abs(decinf)<=ip)&&(ip<=pi-abs(decinf))  % Check if capable of coplanar periapse transfer
    vphdmag=sqrt(2*mud/rppdmag+vinfdmag2);
    if dvopt == 1
            if vphdmag<vppdmag
                dvd=0; 
            else
                dvd=abs(vppdmag-vphdmag); %delta-V in planet frame  
            end
    else
        dvd=abs(vppdmag-vphdmag); %delta-V in planet frame 
    end
else
    ah= -mud/vinfdmag2;
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
        [~,dvd]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,rppdmag,vppdmag,ip,op,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    elseif ip==pi % If retrograde Equatorial:
                 % set RAAN = 0 and w = whatTrue 
        op=0;
        wpi = -rainf-acos(-1/(1-rppdmag/ah))*cos(decinf); %initial guess (rad)
    %     wp =-rainf-acos(-1/(1-rppdmag/ah))*cos(decinfp) %initial guess (rad)
        [~,dvd]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,rppdmag,vppdmag,ip,op,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    else
        op=atan2(node(2),node(1));
    %     ip+decinfp,asin((1/(1-rppdmag/ah)))
    %     wpi=pi-asin((1/(1-rppdmag/ah))/sin(ip+decinfp))
    %     wpi=pi-asin((1/(1-rppdmag/ah))*cos(ip+decinfp))%initial guess (rad)
        wpi=-sign(decinf)*pi/2+pi/4;
        [~,dvd]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,rppdmag,vppdmag,ip,op,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    %    tic
    %    [wp,fval]=fminbnd(@(x)dvmin(x,vinfdhat,vinfdmag,mud,ap,ep,ip,op),0,2*pi,...
    %        optimset('Display','off','TolFun',1e-8,'TolX',1e-8))
    %    toc
    end

end
end
function f = dvmin(x,vinfdhat,vinfdmag,mud,rppmag,vppmag,ip,op,dvopt)

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
    rpphat = rpp/rppmag;
    
% Determine the departure velocity (same location) vector
D=sqrt(mud/(rppmag*(1+rpphat'*vinfdhat))+vinfdmag*vinfdmag/4);
vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rpphat;
vphdmag=sqrt(vphd'*vphd);
if dvopt == 1
    if vphdmag>vppmag
        f=sqrt((vpp-vphd)'*(vpp-vphd));  %basic DV to match entry velocity
    else 
        costheta=(vphd'*vpp/(vphdmag*vppmag));
        vphmag=vphdmag/costheta; %can potentially be undefined if cos(theta)=0
        if vphmag<vppmag %check if perpendicular DV aligns launch but < vppmag
            f=sqrt(abs(vphmag*vphmag-vphdmag*vphdmag)); %perpendicular DV
        else
            f=sqrt((vpp-vphd)'*(vpp-vphd));  %basic DV to match entry velocity
        end
    end
else
    f=sqrt((vpp-vphd)'*(vpp-vphd));  
end    
end