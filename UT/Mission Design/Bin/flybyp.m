function [dvd,rfb,rpfb,vin,vout,rpin,rpout,vpin,vpout,betain,betaout,deltain,...
     deltaout,tsoiin,tsoiout]=flybyp(rpmag,rP,vinfinhat,vinfinmag,vinfouthat,vinfoutmag,...
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
ein=1+rpmag*vinfinmag2/mu;
eout=1+rpmag*vinfoutmag2/mu;
ain = -mu/vinfinmag2;    %semi-transverse axis of hyperbolic trajectory
aout= -mu/vinfoutmag2;
phig=asin(1/eout)+asin(1/ein); %maximum natural gravitational turn angle 
phir=acos(vinfouthat'*vinfinhat); %dot prod
rsoib=rsoi; %set upper bound of transfer radius
 

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
    
if phir>phig
    if vinfoutmag>=vinfinmag %Constrained RF Case
%      1 
hhmag=sqrt(ain*(1-ein*ein)*mu);
hh=hhat*hhmag;
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

%Find Classic Orbital Elements for hyperbolic trajectory  
n= [-hh(2);hh(1);0];% cross product (zhat,h)
nmag = sqrt(n'*n);
pin = ain*(1-ein*ein);
iin= acos(hh(3)/hhmag);
evec = rpinhat*ein;
 if (abs(iin)<zero)||(pi-abs(iin)<zero) 
    Oin=0;         % set o = 0 and w = whatTrue 
    win = acos(evec(1)/ein);
    if ((evec(2)<0)&&(abs(iin)<zero))||((-evec(2)<0)&&(pi-abs(iin)<zero))
        win = 2*pi -win;
    end
 else
      Oin = acos(n(1)/nmag);
    if n(2)<0
        Oin = 2*pi-Oin;
    end
    win = acos((n'*evec)/(nmag*ein));
    if evec(3)<0
        win = 2*pi-win;
    end
 end
 % Rotation 3x2 matrix for converting between pqw and ijk frame
 cOin=cos(Oin);cwin=cos(win);ciin=cos(iin);
 sOin=sin(Oin);swin=sin(win);siin=sin(iin);
Rijkpqw = [(cOin*cwin-sOin*swin*ciin) (-cOin*swin-sOin*cwin*ciin);
           (sOin*cwin+cOin*swin*ciin) (-sOin*swin+cOin*cwin*ciin) ;
           (swin*siin)  (cwin*siin)];

taub=acos((ain*(1-ein*ein)/rsoib-1)/ein); %True Anomaly upper bound at rsoi

f=Inf; %initialize f root zero-function
tata=phir+pi-acos(-1/ein)-acos(-1/eout);
H=aout/ain*(1-eout^2)/(1-ein^2)/eout;
ta0=tata/2; %initialize true anomaly

while abs(f)>1e-8 
    K=H*(1+ein*cos(ta0))-1/eout;
    f=ta0-tata+acos(K);
    df=1+H*ein*sin(ta0)/sqrt(1-K^2);
    ddf=(-K*(df-1)^2+H*ein*cos(ta0))/sqrt(1-K^2);
%        tamin=ta0-f/df;
    talb=ta0-2*f*df/(2*df^2-f*ddf);
    ta0=talb;
end


[tafb,fval]=fminbnd(@(ta)dvin2out(ta,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag),...
talb,taub,optimset('Display','off','TolFun',1e-8,'TolX',1e-8));

costa=cos(tafb);
sinta=sin(tafb);
rpqw = [pin*costa/(1+ein*costa);
        pin*sinta/(1+ein*costa)];
vpqw = [-sqrt(mu/pin)*sinta;
        sqrt(mu/pin)*(ein+costa)];
rfb = (Rijkpqw*rpqw);  %r in ijk frame km
vin = (Rijkpqw*vpqw);   %v in ijk frame km/sec
rfbmag=sqrt(rfb'*rfb);
rfbhat= rfb/rfbmag;
costhetar=(rfbhat'*vinfouthat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;  


h3=rfb(1)*vout(2)-rfb(2)*vout(1);
if sign(hhat(3))~=sign(h3)&&((tafb+acos(-1/ein)-pi)<phir)
D=-D; 
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat; 
end

F = acosh((-rsoib/ain+1)/ein);
Mh = ein*sinh(F)-F;
tsoiin = hhmag*hhmag*hhmag*Mh/(mu*mu*(ein*ein-1)^(3/2));

evec=((vout'*vout)*rfb-(rfb'*vout)*vout)/mu-rfbhat;
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
hh=hhat*hhmag;
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
%Find Classic Orbital Elements for hyperbolic trajectory  
n= [-hh(2);hh(1);0];% cross product (zhat,h)
nmag = sqrt(n'*n);
pout = aout*(1-eout*eout);
iout= acos(hh(3)/hhmag);
evec = rpouthat*eout;
 if (abs(iout)<zero)||(pi-abs(iout)<zero) 
    Oout=0;         % set o = 0 and w = whatTrue 
    wout = acos(evec(1)/eout);
    if ((evec(2)<0)&&(abs(iout)<zero))||((-evec(2)<0)&&(pi-abs(iout)<zero))
        wout = 2*pi -wout;
    end
 else
      Oout = acos(n(1)/nmag);
    if n(2)<0
        Oout = 2*pi-Oout;
    end
    wout = acos((n'*evec)/(nmag*eout));
    if evec(3)<0
        wout = 2*pi-wout;
    end
 end
 % Rotation 3x2 matrix for converting between pqw and ijk frame
 cOout=cos(Oout);cwout=cos(wout);ciout=cos(iout);
 sOout=sin(Oout);swout=sin(wout);siout=sin(iout);
Rijkpqw = [(cOout*cwout-sOout*swout*ciout) (-cOout*swout-sOout*cwout*ciout);
           (sOout*cwout+cOout*swout*ciout) (-sOout*swout+cOout*cwout*ciout) ;
           (swout*siout)  (cwout*siout)];

talb=-acos((aout*(1-eout*eout)/rsoib-1)/eout); %True Anomaly upper bound at rsoi


f=Inf; %initialize f root zero-function
tata=phir+pi-acos(-1/ein)-acos(-1/eout);
H=ain/aout*(1-ein^2)/(1-eout^2)/ein;
ta0=tata/2; %initialize true anomaly

while abs(f)>1e-8 
    K=H*(1+eout*cos(ta0))-1/ein;
    f=ta0-tata+acos(K);
    df=1+H*eout*sin(ta0)/sqrt(1-K^2);
    ddf=(-K*(df-1)^2+H*eout*cos(ta0))/sqrt(1-K^2);
%        tamin=ta0-f/df;
    tamin=ta0-2*f*df/(2*df^2-f*ddf);
    ta0=tamin;
end

taub=-tamin;


[tafb,fval]=fminbnd(@(ta)dvout2in(ta,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag),...
talb,taub,optimset('Display','off','TolFun',1e-8,'TolX',1e-8));
costa=cos(tafb);
sinta=sin(tafb);
rpqw = [pout*costa/(1+eout*costa);
        pout*sinta/(1+eout*costa)];
vpqw = [-sqrt(mu/pout)*sinta;
        sqrt(mu/pout)*(eout+costa)];
rfb = (Rijkpqw*rpqw);  %r in ijk frame km
vout = (Rijkpqw*vpqw);   %v in ijk frame km/sec
rfbmag=sqrt(rfb'*rfb);
rfbhat= rfb/rfbmag;
costhetar=(-rfbhat'*vinfinhat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    

h3=rfb(1)*vin(2)-rfb(2)*vin(1);
if sign(hhat(3))~=sign(h3)&&((abs(tafb-acos(-1/eout))-pi)<phir)
D=-D;
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    
end

F = acosh((-rsoib/aout+1)/eout);
Mh = eout*sinh(F)-F;
tsoiout = hhmag*hhmag*hhmag*Mh/(mu*mu*(eout*eout-1)^(3/2));

evec=((vin'*vin)*rfb-(rfb'*vin)*vin)/mu-rfbhat;
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
sigma=(pi-phir)/4;
delta=atan((vinfinmag-vinfoutmag)/(vinfinmag+vinfoutmag)*tan(sigma));
F=(sigma-delta);
% G=(sigma+delta);
L=sqrt(2*mu/rpmag);
rfbmag=rpmag*L*L*sin(F)^2/(vinfoutmag2*cos(delta)^2*(2*cos(sigma)^2-cos(delta)^2));
% DV=(vinfinmag+vinfoutmag)*sin(delta)
rfbhat=rotuv(vinfinhat,hhat,-2*F);
costhetar=(rfbhat'*vinfouthat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;
rfb=rfbhat*rfbmag;
evec=((vout'*vout)*rfb-(rfb'*vout)*vout)/mu-rfbhat;
rpfbmag=aout*(1-sqrt(evec'*evec));


if rpfbmag<rpmag || rpfbmag>rsoi || ~optr %check if RF Case or fixed R option
    
    if (rpfbmag>rsoi)&&(optr)
       rpmag=rsoi; %adjust periapse just below rsoi
       eout=1+rpmag*vinfoutmag2/mu;
       rsoib=2*rsoi; %adjust upper bound transfer radius
    end
    %Optimized constrained RF Case
    hhmag=sqrt(aout*(1-eout*eout)*mu);
    hh=hhat*hhmag;
    betaout = acos(1/eout);
    deltaout = hhmag/vinfoutmag;
%     rpouthat = rotuv(-vinfouthat,hhat,betaout);
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
%Find Classic Orbital Elements for hyperbolic trajectory  
    n= [-hh(2);hh(1);0];% cross product (zhat,h)
    nmag = sqrt(n'*n);
    pout = aout*(1-eout*eout);
    iout= acos(hh(3)/hhmag);
    evec = rpouthat*eout;
     if (abs(iout)<zero)||(pi-abs(iout)<zero) 
        Oout=0;         % set o = 0 and w = whatTrue 
        wout = acos(evec(1)/eout);
        if ((evec(2)<0)&&(abs(iout)<zero))||((-evec(2)<0)&&(pi-abs(iout)<zero))
            wout = 2*pi -wout;
        end
     else
          Oout = acos(n(1)/nmag);
        if n(2)<0
            Oout = 2*pi-Oout;
        end
        wout = acos((n'*evec)/(nmag*eout));
        if evec(3)<0
            wout = 2*pi-wout;
        end
     end
     % Rotation 3x2 matrix for converting between pqw and ijk frame
 cOout=cos(Oout);cwout=cos(wout);ciout=cos(iout);
 sOout=sin(Oout);swout=sin(wout);siout=sin(iout);
Rijkpqw = [(cOout*cwout-sOout*swout*ciout) (-cOout*swout-sOout*cwout*ciout);
           (sOout*cwout+cOout*swout*ciout) (-sOout*swout+cOout*cwout*ciout) ;
           (swout*siout)  (cwout*siout)];
        
    talb=-acos((aout*(1-eout*eout)/rsoib-1)/eout); %True Anomaly upper bound at rsoi
    taub=0; %True Anomaly lower bound at periapse

    [tafb,fval]=fminbnd(@(ta)dvout2in(ta,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag),...
    talb,taub,optimset('Display','off','TolFun',1e-8,'TolX',1e-8));
    costa=cos(tafb);
    sinta=sin(tafb);
    rpqw = [pout*costa/(1+eout*costa);
            pout*sinta/(1+eout*costa)];
    vpqw = [-sqrt(mu/pout)*sinta;
            sqrt(mu/pout)*(eout+costa)];
    rfb = (Rijkpqw*rpqw);  %r in ijk frame km
    vout = (Rijkpqw*vpqw);   %v in ijk frame km/sec
    rfbmag=sqrt(rfb'*rfb);
    rfbhat= rfb/rfbmag;
    costhetar=(-rfbhat'*vinfinhat);
    D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag2/4);
    vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;  
    
else
    eout = sqrt(evec'*evec); %reassign eout for F Case
    hhmag=sqrt(aout*(1-eout*eout)*mu);
    betaout = acos(1/eout);
    deltaout = hhmag/vinfoutmag;
    rpouthat=evec/eout;
    rpfb=rpfbmag*rpouthat;
    rpout=rpfb;
    vpoutmag=sqrt(2*mu/rpfbmag+vinfoutmag2);
    % vpouthat=rotuv(vinfouthat,hhat,betaout-pi/2);
    vpouthat=[hhat(2)*rpouthat(3) - hhat(3)*rpouthat(2);...%cross product 
              hhat(3)*rpouthat(1) - hhat(1)*rpouthat(3);...
              hhat(1)*rpouthat(2) - hhat(2)*rpouthat(1)];
    vpout=vpouthat*vpoutmag;
    costhetar=(-rfbhat'*vinfinhat);
    D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag2/4);
    vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;  
    
end
F = acosh((-rsoib/aout+1)/eout);
Mh = eout*sinh(F)-F;
tsoiout = hhmag*hhmag*hhmag*Mh/(mu*mu*(eout*eout-1)^(3/2));

evec=((vin'*vin)*rfb-(rfb'*vin)*vin)/mu-rfbhat;
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
sigma=(pi-phir)/4;
delta=atan((vinfinmag-vinfoutmag)/(vinfinmag+vinfoutmag)*tan(sigma));
F=(sigma-delta);
% G=(sigma+delta);
L=sqrt(2*mu/rpmag);
rfbmag=rpmag*L*L*sin(F)^2/(vinfoutmag2*cos(delta)^2*(2*cos(sigma)^2-cos(delta)^2));
% DV=(vinfinmag+vinfoutmag)*sin(delta)
rfbhat=rotuv(vinfinhat,hhat,-2*F);
costhetar=(-rfbhat'*vinfinhat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;
rfb=rfbhat*rfbmag;
evec=((vin'*vin)*rfb-(rfb'*vin)*vin)/mu-rfbhat;
rpfbmag=ain*(1-sqrt(evec'*evec));

if rpfbmag<rpmag  || rpfbmag>rsoi || ~optr %check if RF Case or fixed R option
    
    if (rpfbmag>rsoi)&&(optr)
       rpmag=rsoi; %adjust periapse just below rsoi
       ein=1+rpmag*vinfinmag2/mu;
       rsoib=2*rsoi; %adjust upper bound transfer radius
    end
    %Optimized constrained RF Case
    hhmag=sqrt(ain*(1-ein*ein)*mu);
    hh=hhat*hhmag;
    betain = acos(1/ein);
    deltain = hhmag/vinfinmag;
%     rpinhat = rotuv(vinfinhat,hhat,-betain);
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
%Find Classic Orbital Elements for hyperbolic trajectory  
    n= [-hh(2);hh(1);0];% cross product (zhat,h)
    nmag = sqrt(n'*n);
    pin = ain*(1-ein*ein);
    iin= acos(hh(3)/hhmag);
    evec = rpinhat*ein;
     if (abs(iin)<zero)||(pi-abs(iin)<zero) 
        Oin=0;         % set o = 0 and w = whatTrue 
        win = acos(evec(1)/ein);
        if ((evec(2)<0)&&(abs(iin)<zero))||((-evec(2)<0)&&(pi-abs(iin)<zero))
            win = 2*pi -win;
        end
     else
          Oin = acos(n(1)/nmag);
        if n(2)<0
            Oin = 2*pi-Oin;
        end
        win = acos((n'*evec)/(nmag*ein));
        if evec(3)<0
            win = 2*pi-win;
        end
     end
     % Rotation 3x2 matrix for converting between pqw and ijk frame
 cOin=cos(Oin);cwin=cos(win);ciin=cos(iin);
 sOin=sin(Oin);swin=sin(win);siin=sin(iin);
Rijkpqw = [(cOin*cwin-sOin*swin*ciin) (-cOin*swin-sOin*cwin*ciin);
           (sOin*cwin+cOin*swin*ciin) (-sOin*swin+cOin*cwin*ciin) ;
           (swin*siin)  (cwin*siin)];
        
    taub=acos((ain*(1-ein*ein)/rsoib-1)/ein); %True Anomaly upper bound at rsoi
    talb=0; %True Anomaly lower bound at periapse

    [tafb,fval]=fminbnd(@(ta)dvin2out(ta,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag),...
    talb,taub,optimset('Display','off','TolFun',1e-8,'TolX',1e-8));
    costa=cos(tafb);
    sinta=sin(tafb);
    rpqw = [pin*costa/(1+ein*costa);
            pin*sinta/(1+ein*costa)];
    vpqw = [-sqrt(mu/pin)*sinta;
            sqrt(mu/pin)*(ein+costa)];
    rfb = (Rijkpqw*rpqw);  %r in ijk frame km
    vin = (Rijkpqw*vpqw);   %v in ijk frame km/sec
    rfbmag=sqrt(rfb'*rfb);
    rfbhat= rfb/rfbmag;
    costhetar=(rfbhat'*vinfouthat);
    D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag2/4);
    vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;  
    
else
    ein = sqrt(evec'*evec); %reassign eout for F Case
    hhmag=sqrt(ain*(1-ein*ein)*mu);
    betain = acos(1/ein);
    deltain = hhmag/vinfinmag;
    rpinhat =evec/ein;
    rpfb=rpfbmag*rpinhat;
    rpin=rpfb;
    vpinmag=sqrt(2*mu/rpfbmag+vinfinmag2);
    % vpinhat=rotuv(vinfinhat,hhat,pi/2-betain);
    vpinhat=[hhat(2)*rpinhat(3) - hhat(3)*rpinhat(2);...%cross product 
              hhat(3)*rpinhat(1) - hhat(1)*rpinhat(3);...
              hhat(1)*rpinhat(2) - hhat(2)*rpinhat(1)];
    vpin=vpinhat*vpinmag;
    costhetar=(rfbhat'*vinfouthat);
    D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag2/4);
    vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;  
    
end
F = acosh((-rsoib/ain+1)/ein);
Mh = ein*sinh(F)-F;
tsoiin = hhmag*hhmag*hhmag*Mh/(mu*mu*(ein*ein-1)^(3/2));

evec=((vout'*vout)*rfb-(rfb'*vout)*vout)/mu-rfbhat;
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

function f = dvin2out(ta,Rijkpqw,vinfoutmag,vinfouthat,hhat3,p,e,mu,phir,rpmag)
    
    
costa=cos(ta);
sinta=sin(ta);
rpqw = p/(1+e*costa)*[costa; sinta];
vpqw =  sqrt(mu/p)*[-sinta; (e+costa)];
rfb = (Rijkpqw*rpqw);  %r in ijk frame km
vfb = (Rijkpqw*vpqw);   %v in ijk frame km/sec
rfbmag=sqrt(rfb'*rfb);
rfbhat= rfb/rfbmag;
costhetar=(rfbhat'*vinfouthat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag*vinfoutmag/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;  
h3=rfb(1)*vout(2)-rfb(2)*vout(1);
if (sign(hhat3)~=sign(h3))&&((ta+acos(-1/e)-pi)<phir)
D=-D;
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat; 
end

        f=sqrt((vfb-vout)'*(vfb-vout));

end

function f = dvout2in(ta,Rijkpqw,vinfinmag,vinfinhat,hhat3,p,e,mu,phir,rpmag)

    
costa=cos(ta);
sinta=sin(ta);
rpqw = [p*costa/(1+e*costa);
        p*sinta/(1+e*costa)];
vpqw = [-sqrt(mu/p)*sinta;
        sqrt(mu/p)*(e+costa)];
rfb = (Rijkpqw*rpqw);  %r in ijk frame km
vfb = (Rijkpqw*vpqw);   %v in ijk frame km/sec
rfbmag=sqrt(rfb'*rfb);
rfbhat= rfb/rfbmag;
costhetar=(-rfbhat'*vinfinhat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag*vinfinmag/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    
h3=rfb(1)*vin(2)-rfb(2)*vin(1);
if (sign(hhat3)~=sign(h3))&&((abs(ta-acos(-1/e))-pi)<phir)
D=-D;
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    
end

        f=sqrt((vfb-vin)'*(vfb-vin));

end