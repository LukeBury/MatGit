function [dvd,rfb,rpfb,vin,vout,rpin,rpout,vpin,vpout,betain,betaout,deltain,...
     deltaout,tsoiin,tsoiout]=flybyunp(rpmag,rP,vinfinhat,vinfinmag,vinfouthat,vinfoutmag,...
    mu,rsoi,optr)
% warning off
% Determines the powered gravity assist trajectory to a parking orbit of

%------INPUTS-----
% rpmag: periapsis radius
% rP(3): Position of arrival planet
% vinfinhat(3): vinf incoming dir vector
% vinfinmag: vinf incoming magnitude
% vinfouthat(3): vinf outgoing dir vector
% vinfoutmag: vinf outgoing magnitude
% mu: gravitational parameter of planet
% rsoi: radius of sphere of influence of planet
% optsoi: optimization of dv
% optr: optimization of rp
%------OUTPUTS-----
% dvd(3): Delta-V required for hyperbolic depature trajectory
% rfb(3): location of dv
% rpfb(3): location of periapsis
% vin(3): velocity incoming at dv
% vout(3): velocity outgoing at dv
% rpin(3): location of incoming periapsis
% rpout(3): location of outgoing periapsis
% vpin(3): Velocity at periapsis incoming
% vpout(3): Velocity at periapsis outgoing
% betain: Angle between periapsis location and vinf asymptote of incoming 
%       hyperbolic orbit
% betaout: Angle between periapsis location and vinf asymptote of outgoing 
%       hyperbolic orbit
% deltain: Aiming radius of incoming hyperbolic trajectory
% deltaout: Aiming radius of outgoing hyperbolic trajectory
% tsoiin: time on incoming hyperbolic orbit between periapsis and rsoi intersection
% tsoiin: time on outgoing hyperbolic orbit between periapsis and rsoi intersection
zero=1e-12;
vinfinmag2=vinfinmag*vinfinmag;
vinfoutmag2=vinfoutmag*vinfoutmag;
ain = -mu/vinfinmag2;    %semi-transverse axis of hyperbolic trajectory
aout= -mu/vinfoutmag2;
ein=1-rpmag/ain;
eout=1-rpmag/aout;
phig=asin(1/eout)+asin(1/ein); %maximum natural gravitational turn angle 
phir=acos(vinfouthat'*vinfinhat); %dot prod
rsoib=rsoi; %set upper bound of transfer radius

% dv=vinfinmag*sqrt(vinfoutmag2/vinfinmag2+1+2*vinfoutmag/vinfinmag*cos(phir))
% dv=vinfinmag*sqrt(vinfoutmag2/vinfinmag2+1-2*vinfoutmag/vinfinmag*cos(phir))

hhat = [vinfinhat(2)*vinfouthat(3) - vinfinhat(3)*vinfouthat(2);...%cross product 
        vinfinhat(3)*vinfouthat(1) - vinfinhat(1)*vinfouthat(3);...
        vinfinhat(1)*vinfouthat(2) - vinfinhat(2)*vinfouthat(1)];
hhat = hhat/sqrt(hhat'*hhat);
    if (sqrt(hhat'*hhat)<zero)||(phir<zero)||(phir>pi-zero)
% hhat = cross(vinfinhat,rP)/norm(cross(vinfinhat,rP));
hhat = [vinfinhat(2)*rP(3) - vinfinhat(3)*rP(2);...%cross product 
        vinfinhat(3)*rP(1) - vinfinhat(1)*rP(3);...
        vinfinhat(1)*rP(2) - vinfinhat(2)*rP(1)];
hhat = hhat/sqrt(hhat'*hhat);    
    end
    
%find optimal periapse radius, if selected
if optr
   rpig = mu/(vinfinmag*vinfoutmag)*(sqrt(2/(1-cos(phir)))-1); %initial guess
   cphi0=cos(phir);
   f=1; %initialize f
   while abs(f)>1e-8
   rpig0=rpig;
   rpig2=rpig*rpig;
   term=sqrt(rpig2*vinfinmag2*vinfoutmag2*(rpig*vinfinmag2+2*mu)*(rpig*vinfoutmag2+2*mu)/mu^4);
   cphi=mu*mu*(term-1)/((rpig*vinfinmag2+mu)*(rpig*vinfoutmag2+mu));
   f=cphi-cphi0;
   df=(rpig2*rpig*vinfinmag^6*vinfoutmag2+rpig2*rpig*vinfinmag2*vinfoutmag^6+3*rpig2*vinfinmag2*vinfinmag2*vinfoutmag2*mu+...
    3*rpig2*vinfinmag2*vinfoutmag2*vinfoutmag2*mu+term*(2*mu*mu*rpig*vinfinmag2*vinfoutmag2+mu^3*vinfinmag2+mu^3*vinfoutmag2)+...
    4*mu^2*rpig*vinfinmag2*vinfoutmag2)/(term*(rpig*vinfinmag2+mu)^2*(rpig*vinfoutmag2+mu)^2);
   rpig=abs(rpig0-f/df);
   end
   if rpig >= rpmag
        if rpig<rsoi
        rpmag = rpig;
        ein=1-rpmag/ain;
        eout=1-rpmag/aout;
       else
        rpmag=rsoi; %adjust periapse to rsoi
        ein=1-rpmag/ain;
        eout=1-rpmag/aout;
        rsoib=2*rsoi;
       end
   end
end
    
    
    
if phir>phig
    if vinfoutmag>=vinfinmag %Constrained RF Case
%      1   
hhmag=sqrt(ain*(1-ein*ein)*mu);
betain = acos(1/ein);
deltain = hhmag/vinfinmag;
% rpinhat = rotuv(vinfinhat,hhat,-betain);
rpinhat= [vinfinhat(2)*hhat(3) - vinfinhat(3)*hhat(2);...%cross product 
        vinfinhat(3)*hhat(1) - vinfinhat(1)*hhat(3);...
        vinfinhat(1)*hhat(2) - vinfinhat(2)*hhat(1)];
rpinhat=-cos(pi-betain)*vinfinhat+sin(pi-betain)*rpinhat;
rpin = rpmag*rpinhat;
rpfb=rpin;
vpinmag=sqrt(2*mu/rpmag+vinfinmag2);
% vpinhat=rotuv(vinfinhat,hhat,pi/2-betain);
vpinhat=[hhat(2)*rpinhat(3) - hhat(3)*rpinhat(2);...%cross product 
          hhat(3)*rpinhat(1) - hhat(1)*rpinhat(3);...
          hhat(1)*rpinhat(2) - hhat(2)*rpinhat(1)];
vpin=vpinhat*vpinmag;

rfb = rpfb;  %r in ijk frame km
vin = vpin;   %v in ijk frame km/sec
costhetar=(rpinhat'*vinfouthat);
D=sqrt(mu/(rpmag*(1+costhetar))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat;  
h3=rfb(1)*vout(2)-rfb(2)*vout(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat; 
end

F = acosh((-rsoib/ain+1)/ein);
Mh = ein*sinh(F)-F;
tsoiin = hhmag*hhmag*hhmag*Mh/(mu*mu*(ein*ein-1)^(3/2));

evec=((vout'*vout)*rfb-(rfb'*vout)*vout)/mu-rpinhat;
eout=sqrt(evec'*evec);
rpouthat=evec/eout;
hhmag=sqrt(aout*(1-eout*eout)*mu);
betaout = acos(1/eout);
deltaout = hhmag/vinfoutmag;
rpoutmag=aout*(1-eout);
rpout=rpouthat*rpoutmag;
vpoutmag=sqrt(2*mu/rpoutmag+vinfoutmag2);
vpouthat=[hhat(2)*rpouthat(3) - hhat(3)*rpouthat(2);...%cross product 
          hhat(3)*rpouthat(1) - hhat(1)*rpouthat(3);...
          hhat(1)*rpouthat(2) - hhat(2)*rpouthat(1)];
vpout=vpouthat*vpoutmag;
F = acosh((-rsoib/aout+1)/eout);
Mh = eout*sinh(F)-F;
% Determine the time from rsoi to hyperbolic (flyby) periapsis 
tsoiout = hhmag*hhmag*hhmag*Mh/(mu*mu*(eout*eout-1)^(3/2));


    elseif vinfoutmag<vinfinmag %Constrained RF Case
%     2        

hhmag=sqrt(aout*(1-eout*eout)*mu);
betaout = acos(1/eout);
deltaout = hhmag/vinfoutmag;
% rpouthat = rotuv(-vinfouthat,hhat,betaout);
rpouthat= [vinfouthat(2)*hhat(3) - vinfouthat(3)*hhat(2);...%cross product 
        vinfouthat(3)*hhat(1) - vinfouthat(1)*hhat(3);...
        vinfouthat(1)*hhat(2) - vinfouthat(2)*hhat(1)];
rpouthat=cos(pi-betaout)*vinfouthat+sin(pi-betaout)*rpouthat;
rpout = rpmag*rpouthat;
rpfb=rpout;
vpoutmag=sqrt(2*mu/rpmag+vinfoutmag2);
% vpouthat=rotuv(vinfouthat,hhat,betaout-pi/2);
vpouthat=[hhat(2)*rpouthat(3) - hhat(3)*rpouthat(2);...%cross product 
          hhat(3)*rpouthat(1) - hhat(1)*rpouthat(3);...
          hhat(1)*rpouthat(2) - hhat(2)*rpouthat(1)];
vpout=vpouthat*vpoutmag;

rfb = rpout;  %r in ijk frame km
vout = vpout;   %v in ijk frame km/sec
costhetar=(-rpouthat'*vinfinhat);
D=sqrt(mu/(rpmag*(1+costhetar))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
h3=rfb(1)*vin(2)-rfb(2)*vin(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
end

F = acosh((-rsoib/aout+1)/eout);
Mh = eout*sinh(F)-F;
tsoiout = hhmag*hhmag*hhmag*Mh/(mu*mu*(eout*eout-1)^(3/2));

evec=((vin'*vin)*rfb-(rfb'*vin)*vin)/mu-rpouthat;
ein=sqrt(evec'*evec);
rpinhat=evec/ein;
hhmag=sqrt(ain*(1-ein*ein)*mu);
betain = acos(1/ein);
deltain = hhmag/vinfinmag;
rpinmag=ain*(1-ein);
rpin=rpinhat*rpinmag;
vpinmag=sqrt(2*mu/rpinmag+vinfinmag2);
vpinhat=[hhat(2)*rpinhat(3) - hhat(3)*rpinhat(2);...%cross product 
          hhat(3)*rpinhat(1) - hhat(1)*rpinhat(3);...
          hhat(1)*rpinhat(2) - hhat(2)*rpinhat(1)];
vpin=vpinhat*vpinmag;
F = acosh((-rsoib/ain+1)/ein);
Mh = ein*sinh(F)-F;
% Determine the time from rsoi to hyperbolic (flyby) periapsis 
tsoiin = hhmag*hhmag*hhmag*Mh/(mu*mu*(ein*ein-1)^(3/2));

    
    
    end
elseif phir<=phig
    if vinfoutmag>=vinfinmag %F case or Constrained RF Case
%     3     
hhmag=sqrt(aout*(1-eout*eout)*mu);
betaout = acos(1/eout);
deltaout = hhmag/vinfoutmag;
% rpouthat = rotuv(-vinfouthat,hhat,betaout);
rpouthat= [vinfouthat(2)*hhat(3) - vinfouthat(3)*hhat(2);...%cross product 
        vinfouthat(3)*hhat(1) - vinfouthat(1)*hhat(3);...
        vinfouthat(1)*hhat(2) - vinfouthat(2)*hhat(1)];
rpouthat=cos(pi-betaout)*vinfouthat+sin(pi-betaout)*rpouthat;
rpout = rpmag*rpouthat;
rpfb=rpout;
vpoutmag=sqrt(2*mu/rpmag+vinfoutmag2);
% vpouthat=rotuv(vinfouthat,hhat,betaout-pi/2);
vpouthat=[hhat(2)*rpouthat(3) - hhat(3)*rpouthat(2);...%cross product 
          hhat(3)*rpouthat(1) - hhat(1)*rpouthat(3);...
          hhat(1)*rpouthat(2) - hhat(2)*rpouthat(1)];
vpout=vpouthat*vpoutmag;

rfb = rpout;  %r in ijk frame km
vout = vpout;   %v in ijk frame km/sec
costhetar=(-rpouthat'*vinfinhat);
D=sqrt(mu/(rpmag*(1+costhetar))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
h3=rfb(1)*vin(2)-rfb(2)*vin(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
end

F = acosh((-rsoib/aout+1)/eout);
Mh = eout*sinh(F)-F;
tsoiout = hhmag*hhmag*hhmag*Mh/(mu*mu*(eout*eout-1)^(3/2));

evec=((vin'*vin)*rfb-(rfb'*vin)*vin)/mu-rpouthat;
ein=sqrt(evec'*evec);
rpinhat=evec/ein;
hhmag=sqrt(ain*(1-ein*ein)*mu);
betain = acos(1/ein);
deltain = hhmag/vinfinmag;
rpinmag=ain*(1-ein);
rpin=rpinhat*rpinmag;
vpinmag=sqrt(2*mu/rpinmag+vinfinmag2);
vpinhat=[hhat(2)*rpinhat(3) - hhat(3)*rpinhat(2);...%cross product 
          hhat(3)*rpinhat(1) - hhat(1)*rpinhat(3);...
          hhat(1)*rpinhat(2) - hhat(2)*rpinhat(1)];
vpin=vpinhat*vpinmag;
F = acosh((-rsoib/ain+1)/ein);
Mh = ein*sinh(F)-F;
% Determine the time from rsoi to hyperbolic (flyby) periapsis 
tsoiin = hhmag*hhmag*hhmag*Mh/(mu*mu*(ein*ein-1)^(3/2));

    elseif vinfoutmag<vinfinmag %F case or Constrained RF Case
%     4
hhmag=sqrt(ain*(1-ein*ein)*mu);
betain = acos(1/ein);
deltain = hhmag/vinfinmag;
% rpinhat = rotuv(vinfinhat,hhat,-betain)
rpinhat= [vinfinhat(2)*hhat(3) - vinfinhat(3)*hhat(2);...%cross product 
        vinfinhat(3)*hhat(1) - vinfinhat(1)*hhat(3);...
        vinfinhat(1)*hhat(2) - vinfinhat(2)*hhat(1)];
rpinhat=-cos(pi-betain)*vinfinhat+sin(pi-betain)*rpinhat;
rpin = rpmag*rpinhat;
rpfb=rpin;
vpinmag=sqrt(2*mu/rpmag+vinfinmag2);
% vpinhat=rotuv(vinfinhat,hhat,pi/2-betain);
vpinhat=[hhat(2)*rpinhat(3) - hhat(3)*rpinhat(2);...%cross product 
          hhat(3)*rpinhat(1) - hhat(1)*rpinhat(3);...
          hhat(1)*rpinhat(2) - hhat(2)*rpinhat(1)];
vpin=vpinhat*vpinmag;

rfb = rpfb;  %r in ijk frame km
vin = vpin;   %v in ijk frame km/sec
costhetar=(rpinhat'*vinfouthat);
D=sqrt(mu/(rpmag*(1+costhetar))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat;  
h3=rfb(1)*vout(2)-rfb(2)*vout(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat; 
end

F = acosh((-rsoib/ain+1)/ein);
Mh = ein*sinh(F)-F;
tsoiin = hhmag*hhmag*hhmag*Mh/(mu*mu*(ein*ein-1)^(3/2));

evec=((vout'*vout)*rfb-(rfb'*vout)*vout)/mu-rpinhat;
eout=sqrt(evec'*evec);
rpouthat=evec/eout;
hhmag=sqrt(aout*(1-eout*eout)*mu);
betaout = acos(1/eout);
deltaout = hhmag/vinfoutmag;
rpoutmag=aout*(1-eout);
rpout=rpouthat*rpoutmag;
vpoutmag=sqrt(2*mu/rpoutmag+vinfoutmag2);
vpouthat=[hhat(2)*rpouthat(3) - hhat(3)*rpouthat(2);...%cross product 
          hhat(3)*rpouthat(1) - hhat(1)*rpouthat(3);...
          hhat(1)*rpouthat(2) - hhat(2)*rpouthat(1)];
vpout=vpouthat*vpoutmag;
F = acosh((-rsoib/aout+1)/eout);
Mh = eout*sinh(F)-F;
% Determine the time from rsoi to hyperbolic (flyby) periapsis 
tsoiout = hhmag*hhmag*hhmag*Mh/(mu*mu*(eout*eout-1)^(3/2));

    end
end
dvd = vout - vin;
end