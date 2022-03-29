function [dvd,tafb,rfbmag]=flybypOPT(rpmag,rP,vinfinhat,vinfinmag,vinfouthat,vinfoutmag,mu,rsoi,optr)
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
% dvd: Delta-V required for hyperbolic depature trajectory

zero=1e-12;
tol=1e-6;
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
% rpinhat = rotuv(vinfinhat,hhat,-acos(1/ein));
nuinf=acos(-1/ein);
rpinhat= [vinfinhat(2)*hhat(3) - vinfinhat(3)*hhat(2);...%cross product 
        vinfinhat(3)*hhat(1) - vinfinhat(1)*hhat(3);...
        vinfinhat(1)*hhat(2) - vinfinhat(2)*hhat(1)];
rpinhat=-cos(nuinf)*vinfinhat+sin(nuinf)*rpinhat;

%Find Classic Orbital Elements for hyperbolic trajectory  
n= [-hh(2);hh(1);0]; % cross product (zhat,h)
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

while abs(f)>1e-7 
    K=H*(1+ein*cos(ta0))-1/eout;
    f=ta0-tata+acos(K);
    df=1+H*ein*sin(ta0)/sqrt(1-K^2);
    ddf=(-K*(df-1)^2+H*ein*cos(ta0))/sqrt(1-K^2);
%        tamin=ta0-f/df;
    talb=ta0-2*f*df/(2*df^2-f*ddf);
    if ~isreal(talb)
       talb=ta0;
       break
    end
    ta0=talb;
end


[tafb,dvd]=fminbnd(@(ta)dvin2out(ta,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag),...
talb,taub,optimset('Display','off','TolFun',tol,'TolX',tol));

%check boundary conditions for edge solutions
% dvperi = dvin2out(talb,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag);
dvsoi = dvin2out(taub,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag);
% [dvd,idv]=min([dvperi dvd dvsoi]);
% tavec=[talb tafb taub];
[dvd,idv]=min([dvd dvsoi]);
tavec=[tafb taub];
tafb=tavec(idv); %select min dv solution


    elseif vinfoutmag<vinfinmag %Constrained RF Case
%     2        

hhmag=sqrt(aout*(1-eout*eout)*mu);
hh=hhat*hhmag;
% rpouthat = rotuv(-vinfouthat,hhat,acos(1/eout));
nuinf=acos(-1/eout);
rpouthat= [vinfouthat(2)*hhat(3) - vinfouthat(3)*hhat(2);...%cross product 
        vinfouthat(3)*hhat(1) - vinfouthat(1)*hhat(3);...
        vinfouthat(1)*hhat(2) - vinfouthat(2)*hhat(1)];
rpouthat=cos(nuinf)*vinfouthat+sin(nuinf)*rpouthat;

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

while abs(f)>1e-7 
    K=H*(1+eout*cos(ta0))-1/ein;
    f=ta0-tata+acos(K);
    df=1+H*eout*sin(ta0)/sqrt(1-K^2);
    ddf=(-K*(df-1)^2+H*eout*cos(ta0))/sqrt(1-K^2);
%        tamin=ta0-f/df;
    tamin=ta0-2*f*df/(2*df^2-f*ddf);
    if ~isreal(tamin) %check for imaginary case
      tamin=ta0;
      break
    end
    ta0=tamin;
end

taub=-tamin;

[tafb,dvd]=fminbnd(@(ta)dvout2in(ta,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag),...
talb,taub,optimset('Display','off','TolFun',tol,'TolX',tol));

%check boundary conditions for edge solutions
% dvperi = dvout2in(taub,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag);
dvsoi = dvout2in(talb,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag);
% [dvd,idv]=min([dvperi dvd dvsoi]);
% tavec=[taub tafb talb];
[dvd,idv]=min([ dvd dvsoi]);
tavec=[ tafb talb];
tafb=tavec(idv); %select min dv solution

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
% DV=(vinfinmag+vinfoutmag)*sin(delta);
rfbhat=rotuv(vinfinhat,hhat,-2*F);
costhetar=(rfbhat'*vinfouthat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;
rfb=rfbhat*rfbmag;
evec=((vout'*vout)*rfb-(rfb'*vout)*vout)/mu-rfbhat;
rpfbmag=aout*(1-sqrt(evec'*evec));

if rpfbmag<rpmag || rpfbmag>rsoi || ~optr %check if RF Case or fixed R option
    if (rpfbmag>rsoi)&&(optr)
       rpmag=rsoi; %adjust periapse to rsoi
       rsoib=2*rsoi; %adjust upper bound of transfer radius
       eout=1+rpmag*vinfoutmag2/mu;
    end
    %Optimized constrained RF Case
    hhmag=sqrt(aout*(1-eout*eout)*mu);
    hh=hhat*hhmag;
%     rpouthat = rotuv(-vinfouthat,hhat,acos(1/eout));
    nuinf=acos(-1/eout);
    rpouthat= [vinfouthat(2)*hhat(3) - vinfouthat(3)*hhat(2);...%cross product 
            vinfouthat(3)*hhat(1) - vinfouthat(1)*hhat(3);...
            vinfouthat(1)*hhat(2) - vinfouthat(2)*hhat(1)];
    rpouthat=cos(nuinf)*vinfouthat+sin(nuinf)*rpouthat;

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

    [tafb,dvd]=fminbnd(@(ta)dvout2in(ta,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag),...
    talb,taub,optimset('Display','off','TolFun',tol,'TolX',tol));

%check boundary conditions for edge solutions
% dvperi = dvout2in(0,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag);
dvsoi = dvout2in(talb,Rijkpqw,vinfinmag,vinfinhat,hhat(3),pout,eout,mu,phir,rpmag);
% [dvd,idv]=min([dvperi dvd dvsoi]);
% tavec=[0 tafb talb];
[dvd,idv]=min([dvd dvsoi]);
tavec=[tafb talb];
tafb=tavec(idv); %select min dv solution

else
    costhetar=(-rfbhat'*vinfinhat);
    D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag2/4);
    vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;  
    dvd = sqrt((vout - vin)'*(vout - vin));
%     DVD = (vinfinmag+vinfoutmag)*sin(delta)

end
    elseif vinfoutmag<vinfinmag %F case or Constrained RF Case
%     4
sigma=(pi-phir)/4;
delta=atan((vinfinmag-vinfoutmag)/(vinfinmag+vinfoutmag)*tan(sigma));
F=(sigma-delta);
% G=(sigma+delta);
L=sqrt(2*mu/rpmag);
rfbmag=rpmag*L*L*sin(F)^2/(vinfoutmag2*cos(delta)^2*(2*cos(sigma)^2-cos(delta)^2));
% DV=(vinfinmag+vinfoutmag)*sin(delta);
rfbhat=rotuv(vinfinhat,hhat,-2*F);
costhetar=(-rfbhat'*vinfinhat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;
rfb=rfbhat*rfbmag;
evec=((vin'*vin)*rfb-(rfb'*vin)*vin)/mu-rfbhat;
rpfbmag=ain*(1-sqrt(evec'*evec));

if rpfbmag<rpmag || rpfbmag>rsoi || ~optr %check if RF Case or fixed R option
        
    if (rpfbmag>rsoi)&&(optr)
       rpmag=rsoi; %adjust periapse to rsoi
       rsoib=2*rsoi; %adjust upper bound of transfer radius
       ein=1+rpmag*vinfinmag2/mu;
    end
    %Optimized constrained RF Case
    hhmag=sqrt(ain*(1-ein*ein)*mu);
    hh=hhat*hhmag;
%     rpinhat = rotuv(vinfinhat,hhat,-acos(1/ein));
    nuinf=acos(-1/ein);
    rpinhat= [vinfinhat(2)*hhat(3) - vinfinhat(3)*hhat(2);...%cross product 
            vinfinhat(3)*hhat(1) - vinfinhat(1)*hhat(3);...
            vinfinhat(1)*hhat(2) - vinfinhat(2)*hhat(1)];
    rpinhat=-cos(nuinf)*vinfinhat+sin(nuinf)*rpinhat;
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

    [tafb,dvd]=fminbnd(@(ta)dvin2out(ta,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag),...
    talb,taub,optimset('Display','off','TolFun',tol,'TolX',tol));

%check boundary conditions for edge solutions
% dvperi = dvin2out(0,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag);
dvsoi=dvin2out(taub,Rijkpqw,vinfoutmag,vinfouthat,hhat(3),pin,ein,mu,phir,rpmag);
% [dvd,idv]=min([dvperi dvd dvsoi]);
% tavec=[0 tafb taub];
[dvd,idv]=min([dvd dvsoi]);
tavec=[tafb taub];
tafb=tavec(idv); %select min dv solution

else
    costhetar=(rfbhat'*vinfouthat);
    D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag2/4);
    vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;  
    dvd = sqrt((vout - vin)'*(vout - vin));
end
    end
end
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
D=sqrt(mu/(rfbmag*(1+rfbhat'*vinfouthat))+vinfoutmag*vinfoutmag/4);
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
D=sqrt(mu/(rfbmag*(1-rfbhat'*vinfinhat))+vinfinmag*vinfinmag/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    
h3=rfb(1)*vin(2)-rfb(2)*vin(1);
if (sign(hhat3)~=sign(h3))&&((abs(ta-acos(-1/e))-pi)<phir)
D=-D;
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;    
end

        f=sqrt((vfb-vin)'*(vfb-vin));

end