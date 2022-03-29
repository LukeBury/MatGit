function [dvd,rfb,rpfb,vin,vout,rpin,rpout,vpin,vpout,betain,betaout,deltain,...
     deltaout,tsoiin,tsoiout]=flybysoi(rpmag,rP,vinfinhat,vinfinmag,vinfouthat,vinfoutmag,...
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
% phig=asin(1/eout)+asin(1/ein); %maximum natural gravitational turn angle 
phir=acos(vinfouthat'*vinfinhat); %dot prod
phimax=2*asin(1/(1+rpmag*min([vinfinmag vinfoutmag])^2/mu)); %maximum available turning angle
phisoi=2*asin(1/(1+rsoi*min([vinfinmag vinfoutmag])^2/mu)); %max available turning angle at rsoi
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
    

    if vinfoutmag>=vinfinmag %Constrained RF Case
%      1   
if (phimax<phir)||~optr
    dvd=vinfoutmag*vinfouthat-vinfinmag*rotuv(vinfinhat,hhat,phimax);
        %use min periapse radius
else
    if phir>phisoi
    dvd=(vinfoutmag-vinfinmag)*vinfouthat;
    rpmag= mu/vinfinmag2*(1/sin(phir/2)-1); %determine periapse for phir
    ein=1-rpmag/ain;
    else
    dvd=vinfoutmag*vinfouthat-vinfinmag*rotuv(vinfinhat,hhat,phisoi);
    rpmag=rsoi; %set to max periapse radius to rsoi  
    ein=1-rpmag/ain;
    rsoib=2*rsoi; %adjust max radius bound
    end
end     
     
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

%Determine state [r v] at rsoi from f and g functions
tasoi=acos((ain*(1-ein*ein)/rsoi-1)/ein); %True Anomaly upper bound at rsoi
vr0=vpin'*rpinhat;
sdta=sin(tasoi);
cdta=cos(tasoi);
f=1-mu*rsoi*(1-cdta)/hhmag^2;
g=rsoi*rpmag*sdta/hhmag;
fdot=mu/hhmag*(vr0/hhmag*(1-cdta)-sdta/rpmag);
gdot=1-mu*rpmag/hhmag^2*(1-cdta);
rfb=f*rpin+g*vpin;
rfbhat=rfb/rsoi;
vin=fdot*rpin+gdot*vpin;

costhetar=(rfbhat'*vinfouthat);
D=sqrt(mu/(rsoi*(1+costhetar))+vinfoutmag2/4);
vout1=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;  
D=-D;
vout2=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat; 
if norm(vout1-vin-dvd) < norm(vout2-vin-dvd)
    vout=vout1;
else
    vout=vout2;
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
F = acosh((-rsoib*1.1/aout+1)/eout); %hyperbolic anomaly just above rsoi
Mh = eout*sinh(F)-F;
% Determine the time from rsoi to hyperbolic (flyby) periapsis 
tsoiout = hhmag*hhmag*hhmag*Mh/(mu*mu*(eout*eout-1)^(3/2));

    elseif vinfoutmag<vinfinmag %Constrained RF Case
%     2        
if (phimax<phir)||~optr
    dvd=vinfoutmag*rotuv(vinfouthat,hhat,-phimax)-vinfinmag*vinfinhat;
           %use min periapse radius
else
    if phir>phisoi
    dvd=(vinfoutmag-vinfinmag)*vinfinhat;
    rpmag= mu/vinfoutmag2*(1/sin(phir/2)-1); %determine periapse for phir
    eout=1-rpmag/aout;
    else
    dvd=vinfoutmag*rotuv(vinfouthat,hhat,-phisoi)-vinfinmag*vinfinhat;
    rpmag=rsoi; %set to max periapse radius  
    eout=1-rpmag/aout;
    rsoib=2*rsoi; %adjust max radius bound
    end       
end

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

%Determine state [r v] at rsoi from f and g functions
tasoi=-acos((aout*(1-eout*eout)/rsoi-1)/eout); %True Anomaly upper bound at rsoi
vr0=vpout'*rpouthat;
sdta=sin(tasoi);
cdta=cos(tasoi);
f=1-mu*rsoi*(1-cdta)/hhmag^2;
g=rsoi*rpmag*sdta/hhmag;
fdot=mu/hhmag*(vr0/hhmag*(1-cdta)-sdta/rpmag);
gdot=1-mu*rpmag/hhmag^2*(1-cdta);
rfb=f*rpout+g*vpout;
rfbhat=rfb/rsoi;
vout=fdot*rpout+gdot*vpout;

costhetar=(-rfbhat'*vinfinhat);
D=sqrt(mu/(rsoi*(1+costhetar))+vinfinmag2/4);
vin1=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    
D=-D;
vin2=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    
if norm(vout-vin1-dvd) < norm(vout-vin2-dvd)
    vin=vin1;
else
    vin=vin2;
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
F = acosh((-rsoib*1.1/ain+1)/ein); %hyperbolic anomaly just above rsoi
Mh = ein*sinh(F)-F;
% Determine the time from rsoi to hyperbolic (flyby) periapsis 
tsoiin = hhmag*hhmag*hhmag*Mh/(mu*mu*(ein*ein-1)^(3/2));
    end
end