function [dvd,dvde,rppd,vppd,rppde,vppde,vinfpl,vPplhat,vdplhat,vphd,vphde,betad,deltad,Tpd,tsoid,coepd]...
    = vinf2hypd(vPd,vd,vinfdhat,vinfdmag,rpl,ra,dec,vinfavmag,ip,altp,mud,rsoid,dvopt)
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
% rpl: equatorial radius of planet
% ra: right ascension of north pole of rotation
% dec: declination of north pole of rotation
% vinfavmag: vinf available from launch
% ip: launch inclination
% altp: altitude of launch periapsis 
% rsoi:  radius of sphere of influence of planet
% mud: gravitational parameter 
% dvopt: optimization of dv at launch periapsis
 
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
%Determine initial launch trajectory parameters
deg2rad=pi/180;
ip=abs(ip)*deg2rad;
rppdmag = rpl+altp;
vinfavmag2=vinfavmag*vinfavmag;
vinfdmag2=vinfdmag*vinfdmag;
ap = -mud/vinfavmag2;
ah= -mud/vinfdmag2;
ep = 1+rppdmag*vinfavmag2/mud;
% tap = 0; %true anomaly; set to zero b/c departure at periapsis
vppdmag=sqrt(2*mud/rppdmag+vinfavmag2);


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
    
    betad = acos(1/(1+rppdmag*vinfdmag2/mud));
    hhat=[cos(alpha)*sin(ip); sin(alpha)*sin(ip); cos(ip)]; 
    %hyperbolic trajectory angular momentum unit vector
    rppdhat= [vinfdhat(2)*hhat(3) - vinfdhat(3)*hhat(2);...%cross product 
            vinfdhat(3)*hhat(1) - vinfdhat(1)*hhat(3);...
            vinfdhat(1)*hhat(2) - vinfdhat(2)*hhat(1)];
    rppdhat=cos(pi-betad)*vinfdhat+sin(pi-betad)*rppdhat; 
    %location unit vector at periapsis of parking orbit
    rppd=rppdmag*rppdhat; %location vector at periapsis of parking orbit
    vppdhat = [hhat(2)*rppdhat(3) - hhat(3)*rppdhat(2);...%cross product 
            hhat(3)*rppdhat(1) - hhat(1)*rppdhat(3);...
            hhat(1)*rppdhat(2) - hhat(2)*rppdhat(1)];
    vppd=sqrt(2*mud/rppdmag+vinfavmag2)*vppdhat;
    vphdmag=sqrt(2*mud/rppdmag+vinfdmag2);
    vphd = vppdhat*vphdmag;
    if dvopt == 1
            if vphdmag<vppdmag
                dvd=[0; 0; 0]; 
            else
                dvd=(vppdmag-vphdmag)*vppdhat; %delta-V in planet frame  
            end
    else
        dvd=(vppdmag-vphdmag)*vppdhat; %delta-V in planet frame 
    end
    
%     [~,~,~,op,wp,~,~]=rv2orbel(rppd, vppd,mud); %don't need to overwrite a,e,i
%Find remaining orbital elements (op and wp)
%     evec=rppdhat*ep;
%     hh=rppdmag*sqrt(vppd'*vppd)*hhat;
%     n=[-hh(2);hh(1);0];
%     nmag=sqrt(n'*n);
%     if ep>zero
%         if (abs(ip)<zero)||(180-abs(ip)<zero) 
%             % If Elliptical Equatorial:
%             op=0;         % set o = 0 and w = whatTrue 
%             wp = acosd(evec(1)/ep);
%             if ((evec(2)<0)&&(abs(ip)<zero))||((-evec(2)<0)&&(180-abs(ip)<zero))
%                 wp = 360 -wp;
%             end
%         else            % Normal cases
%             op = acosd(n(1)/nmag);
%             if n(2)<0
%                 op = 360-op;
%             end
%             wp = acosd((n'*evec)/(nmag*ep));
%             if evec(3)<0
%                 wp = 360-wp;
%             end
%         end 
%     else
%         if ((abs(ip)>zero))&&(180-abs(ip)>zero)    % If Circular Inclined: 
%             %replaced u/ta w/ wp to determine orbit periapse location
%             if (n'*rppd)/(nmag*rppdmag)>1.0
%                 wp=0; %u
%             else
%                 wp = acosd((n'*rppd)/(nmag*rppdmag)); %u
%             end
%             if rppd(3)<0
%                 wp = 360-wp;
%             end
%             op = acosd(n(1)/nmag);
%             if n(2)<0
%                 op = 360-op;
%             end
%         else
%             % If Circular Equatorial: 
%             op=0;
%             wp = acosd(rppdhat(1)); %lamdaTrue
% 
%             if (180-abs(ip)<zero) 
%                 wp = -wp;  %lamdaTrue
%             end
% 
%             if (rppd(2)<0)
%                 wp = 360-wp;  %lamdaTrue 
%             end
%         end
%     end
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
        [wp,fval]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,ap,ep,ip,op,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    elseif ip==pi % If retrograde Equatorial:
                 % set RAAN = 0 and w = whatTrue 
        op=0;
        wpi = -rainf-acos(-1/(1-rppdmag/ah))*cos(decinf); %initial guess (rad)
    %     wp =-rainf-acos(-1/(1-rppdmag/ah))*cos(decinfp) %initial guess (rad)
        [wp,fval]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,ap,ep,ip,op,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    else                           
        op=atan2(node(2),node(1));
    %     ip+decinfp,asin((1/(1-rppdmag/ah)))
    %     wpi=pi-asin((1/(1-rppdmag/ah))/sin(ip+decinfp))
    %     wpi=pi-asin((1/(1-rppdmag/ah))*cos(ip+decinfp))%initial guess (rad)
        wpi=-sign(decinf)*pi/2+pi/4;
        [wp,fval]=fminunc(@(x)dvmin(x,vinfdhat,vinfdmag,mud,ap,ep,ip,op,dvopt),wpi,...
            optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));
    %    tic
    %    [wp,fval]=fminbnd(@(x)dvmin(x,vinfdhat,vinfdmag,mud,ap,ep,ip,op),0,2*pi,...
    %        optimset('Display','off','TolFun',1e-8,'TolX',1e-8))
    %    toc
    end

    [rppd, vppd]=orbel2rvRAD(mud,ap,ep,ip,op,wp,0);
    %parking orbit state at departure in planet frame
    %     rppdmag = sqrt(rppd'*rppd);
    rppdhat = rppd/rppdmag; 
 
    % Determine the departure velocity (same location) vector
D=sqrt(mud/(rppdmag*(1+rppdhat'*vinfdhat))+vinfdmag2/4);
vphd=(D+0.5*vinfdmag)*vinfdhat+(D-0.5*vinfdmag)*rppdhat;
vphdmag=sqrt(vphd'*vphd);

    if dvopt == 1
        if vphdmag>vppdmag
            dvd=vppd-vphd;  %basic DV to match entry velocity
        else 
            costheta=(vphd'*vppd/(vphdmag*vppdmag));
            vphmag=vphdmag/costheta; %can potentially be undefined if cos(theta)=0
            if vphmag<vppdmag %check if perpendicular DV aligns entry but < vppmag
                %perpendicular DV
                dvd=(vppd/vppdmag)*vphmag-vphd;
            else
                dvd=vppd-vphd;  %basic DV to match entry velocity
            end
        end
    else
        dvd=vppd-vphd;  %basic DV to match entry velocity
    end
end

    % Determine hyperbolic trajectory data
    ev = vphd'*vphd*rppd/mud-rppd'*vphd*vphd/mud-rppdhat;
    eh = sqrt(ev'*ev);
    rph = ah*(1-eh);
    deltad = rph*sqrt(1-2*ah/rph);
    betad = acos(1/(1-rph/ah));
    Tpd = 2*pi*(sqrt(rph*rph*rph/mud));
    
coepd = [rppdmag 0 ip/deg2rad op/deg2rad wp/deg2rad 0]; %orbital elements of theoretical parking orbit
    
rot=rot'; %transformation from planetary frame to ecliptic frame
dvde=rot*dvd; %delta-V in ecliptic frame
rppde=rot*rppd; %hyperbolic periapse in ecliptic frame (planet centered)
vppde=rot*vppd; %parking orbit periapse velocity in ecliptic frame 
vphde=rot*vphd; %hyperbolic periapse velocity in ecliptic frame 

% Determine the time from periapsis of hyperbolic trajectory to rsoi
hhmag = sqrt(mud*ah*(1-eh*eh));
F = acosh((1-rsoid/ah)/eh);
Mh = eh*sinh(F)-F;
tsoid = hhmag*hhmag*hhmag*Mh/(mud*mud*(eh*eh-1)^(3/2));

end
function f = dvmin(x,vinfdhat,vinfdmag,mud,ap,ep,ip,op,dvopt)

    [rpp, vpp] = orbel2rvRAD(mud,ap,ep,ip,op,x,0);
    rppmag = sqrt(rpp'*rpp);
    vppmag=sqrt(vpp'*vpp);
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