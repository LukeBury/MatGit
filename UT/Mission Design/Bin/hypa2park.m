function [dva,dvae,rppa,vppa,rppae,vppae,vinfpl,vPplhat,vaplhat,vpha,vphae,betaa,deltaa,Tpa,tsoia,coepa]...
    = hypa2park(vPa,va,vinfahat,vinfamag,ra,dec,coepa,mua,rsoia)
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
% coepa(6): Classical orbital elements of parking orbit (Determined for minimum
%       arrival delta-V

% zero=1e-12;
ap = coepa(1);
ep = coepa(2);
ip = abs(coepa(3));
% op = coepa(4);
% wp = coepa(5);
% tap = 0; %true anomaly; set to zero b/c arrival at periapsis
deg2rad=pi/180;
vinfamag2=vinfamag*vinfamag;
ah=-mua/vinfamag2;
rppamag = ap*(1-ep);
ip=mod(ip,360);
if ip>180 %adjust retrograde inclination to b/w 90-180 deg
    ip=360-ip;
end
ip=ip*deg2rad;

% Transform v-inf from ecliptic to planetocentric equatorial
eps= 0; %(deg) Obliquity of ecliptic at J2000

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
decinf=asin(vinfahat(3)); %declination of v-inf vector planet frame
rainf=atan2(vinfahat(2),vinfahat(1)); %right ascension of v-inf vector planet frame
if (abs(decinf)<=ip)&&(ip<=pi-abs(decinf))  % Check if capable of coplanar periapse transfer
    if ip==0
        alpha=pi+rainf; % or  pi-rainf-B
        op=0; %RAAN undefined
        wp =alpha+acos(-1/(1-rppamag/ah)); % argument of periapse
    elseif ip==pi    
        alpha=pi+rainf; % or  pi-rainf-B
        op=0; %RAAN undefined
        wp =-alpha+acos(-1/(1-rppamag/ah)); % argument of periapse
    else
        alpha=pi+rainf-acos(tan(decinf)*cot(ip)); % or pi-rainf-B  
        op = alpha+pi/2; %detailed eqn below
%         op=pi+rainf+asin(cot(pi/2-decinf)/tan(ip))
        wp=acos(sin(decinf)/sin(ip))+asin(1/(1-rppamag/ah));
%       RAAN=360+rainf*180/pi-asind(cot(B)/tand(ip)) %Other solution
%       W=-acosd(cos(B)/sind(ip))+asind(1/(1+rppamag*vinfamag2/mua)) %Other solution
    end

    betaa = acos(1/(1-rppamag/ah));
    hhat=[cos(alpha)*sin(ip); sin(alpha)*sin(ip); cos(ip)]; 
    %hyperbolic trajectory angular momentum unit vector
    rppahat= [vinfahat(2)*hhat(3) - vinfahat(3)*hhat(2);...%cross product 
            vinfahat(3)*hhat(1) - vinfahat(1)*hhat(3);...
            vinfahat(1)*hhat(2) - vinfahat(2)*hhat(1)];
    rppahat=-cos(pi-betaa)*vinfahat+sin(pi-betaa)*rppahat; %negative sign b/c arrival
%     rotuv(vinfahat,hhat,-betaa)
    %location unit vector at periapsis of parking orbit
    rppa=rppamag*rppahat; %location vector at periapsis of parking orbit
    vppa = [hhat(2)*rppahat(3) - hhat(3)*rppahat(2);...%cross product 
            hhat(3)*rppahat(1) - hhat(1)*rppahat(3);...
            hhat(1)*rppahat(2) - hhat(2)*rppahat(1)];
    vppa=sqrt(2*mua/rppamag-mua/ap)*vppa;

%     [~,~,~,op,wp,~,~]=rv2orbel(rppa,vppa,mua); %don't need to overwrite a,e,i
%     evec=rppahat*ep;
%     hh=rppamag*sqrt(vppa'*vppa)*hhat;
%     n=[-hh(2);hh(1);0];
%     nmag=sqrt(n'*n);
%     if ep>zero
%         if ((ip)<zero)||(180-(ip)<zero) 
%             % If Elliptical Equatorial:
%             op=0;         % set o = 0 and w = whatTrue 
%             wp = acosd(evec(1)/ep);
%             if ((evec(2)<0)&&((ip)<zero))||((-evec(2)<0)&&(180-(ip)<zero))
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
%         if (((ip)>zero))&&(180-(ip)>zero)    % If Circular Inclined: 
%             %replaced u/ta w/ wp to determine orbit periapse location
%             if (n'*rppa)/(nmag*rppamag)>1.0
%                 wp=0; %u
%             else
%                 wp = acosd((n'*rppa)/(nmag*rppamag)); %u
%             end
%             if rppa(3)<0
%                 wp = 360-wp;
%             end
%             op = acosd(n(1)/nmag);
%             if n(2)<0
%                 op = 360-op;
%             end
%         else
%             % If Circular Equatorial: 
%             op=0;
%             wp = acosd(rppahat(1)); %lamdaTrue
% 
%             if (180-(ip)<zero) 
%                 wp = -wp;  %lamdaTrue
%             end
% 
%             if (rppa(2)<0)
%                 wp = 360-wp;  %lamdaTrue 
%             end
%         end
%     end

else
    decinfp=pi/2-decinf; %ajusted declination from pole
    hhat=-sign(decinf)*sin(ip)/sin(decinfp)*vinfahat+[0; 0; sin(ip+decinfp)/sin(decinfp)];

    if ip<=pi/2
        node=[-hhat(2); hhat(1); 0]/sin(ip);
    else
        node=[hhat(2); -hhat(1); 0]/sin(ip);
    end

if ip==0 % If Prograde Equatorial: 
        % set w = 0, RAAN = 0, and TA (or w) = lamdaTrue 
    op=0;
    wpi = (pi+rainf)+acos(-1/(1-rppamag/ah))*cos(decinf); %initial guess (rad)
%      wp =rainf-acos(-1/(1-rppamag/ah))*cos(decinfp) %initial guess (rad)
    [wp,fval]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,ap,ep,ip,op),wpi,...
        optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));

elseif ip==pi % If retrograde Equatorial:
             % set RAAN = 0 and w = whatTrue 
    op=0;
    wpi = -(pi+rainf)+acos(-1/(1-rppamag/ah))*cos(decinf); %initial guess (rad)
%     wp =-rainf-acos(-1/(1-rppamag/ah))*cos(decinfp) %initial guess (rad)
    [wp,fval]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,ap,ep,ip,op),wpi,...
        optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));

else
    op=atan2(node(2),node(1));
%     ip+decinfp,asin((1/(1-rppamag/ah)))
%     wpi=-asin(-(1/(1-rppamag/ah))/sin(ip+decinfp))%initial guess (rad)
%     wp=2*pi-asin((1/(1-rppamag/ah))/sin(ip+decinfp))%initial guess (rad)
%     wpi=(pi/2-asin(1/(1-rppamag/ah)))*cos(ip+decinfp)
%     phi=asin(cos(ip+decinfp)/cos(asin((1/(1-rppamag/ah)))))
%     wpi=cos(phi)
    wpi=sign(decinf)*pi/2-pi/4;
%     tic
    [wp,fval]=fminunc(@(x)dvmin(x,vinfahat,vinfamag,mua,ap,ep,ip,op),wpi,...
        optimoptions('fminunc','Display','off','TolFun',1e-8,'TolX',1e-8,'Algorithm','quasi-newton'));

    %    toc
%     tic
%     if decinf>0
%         xmin=0;
%         xmax=pi/2;
%     else
%         xmin=-pi/2;
%         xmax=0;
%     end
%         
%    [wp,fval]=fminbnd(@(x)dvmin(x,vinfahat,vinfamag,mua,ap,ep,ip,op),xmin,xmax,...
%        optimset('Display','off','TolX',1e-12))
%    toc
end
    [rppa, vppa]=orbel2rvRAD(mua,ap,ep,ip,op,wp,0);
    %parking orbit state at arrival in planet frame
%     rppamag = sqrt(rppa'*rppa);
    rppahat = rppa/rppamag;
% wpi=asin(sin(acos(-(vppa/norm(vppa))'*vinfahat))/sin(ip+decinfp))

end

% Determine the departure velocity (same location) vector
D=sqrt(mua/(rppamag*(1-rppahat'*vinfahat))+vinfamag2/4);
vpha=(D+0.5*vinfamag)*vinfahat-(D-0.5*vinfamag)*rppahat;

% Determine hyperbolic trajectory data
    ev=vpha'*vpha*rppa/mua-rppa'*vpha*vpha/mua-rppahat;
    eh = sqrt(ev'*ev);
    rph=ah*(1-eh);
    deltaa = rph*sqrt(1-2*ah/rph);
    betaa = acos(1/(1-rph/ah));
    Tpa=2*pi*(sqrt(ap*ap*ap/mua));

coepa = [ap ep ip/deg2rad op/deg2rad wp/deg2rad 0]; %orbital elements in planet frame
dva = vppa-vpha; %delta-V in planet frame
rot=rot'; %transformation from planetary frame to ecliptic frame
dvae=rot*dva; %delta-V in ecliptic frame
rppae=rot*rppa; %hyperbolic periapse in ecliptic frame (planet centered)
vppae=rot*vppa; %parking orbit periapse velocity in ecliptic frame 
vphae=rot*vpha; %hyperbolic periapse velocity in ecliptic frame 


% Determine the time from periapsis of hyperbolic trajectory to rsoi
hhmag = sqrt(mua*ah*(1-eh*eh));
F = acosh((1-rsoia/ah)/eh);
Mh = eh*sinh(F)-F;
tsoia = hhmag*hhmag*hhmag*Mh/(mua*mua*(eh*eh-1)^(3/2));
if ep>1.0
Tpa=tsoia;    
end

function f = dvmin(x,vinfahat,vinfamag,mua,ap,ep,ip,op)

    [rpp, vpp] = orbel2rvRAD(mua,ap,ep,ip,op,x,0);
    rppmag = sqrt(rpp'*rpp);
    rpphat = rpp/rppmag;
    
% Determine the departure velocity (same location) vector
D=sqrt(mua/(rppmag*(1-rpphat'*vinfahat))+vinfamag*vinfamag/4);
vpha=(D+0.5*vinfamag)*vinfahat-(D-0.5*vinfamag)*rpphat;   
    
f=sqrt((vpha-vpp)'*(vpha-vpp));
end
end
