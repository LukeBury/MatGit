function [dvd]=flybypOPT(rpmag,rP,vinfinhat,vinfinmag,vinfouthat,vinfoutmag,mu,rsoi,optr)
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
tatol=1e-6;
dvtol=1e-6;
ftol=1e-7;
zero=1e-12;
tazero=1e-6;
iter=0;
vinfinmag2=vinfinmag*vinfinmag;
vinfoutmag2=vinfoutmag*vinfoutmag;
ein=1+rpmag*vinfinmag2/mu;
eout=1+rpmag*vinfoutmag2/mu;
ain = -mu/vinfinmag2;    %semi-transverse axis of hyperbolic trajectory
aout= -mu/vinfoutmag2;
pin=ain*(1-ein^2);
pout=aout*(1-eout^2);
phig=asin(1/eout)+asin(1/ein); %maximum natural gravitational turn angle 
phir=acos(vinfouthat'*vinfinhat); %dot prod



    
if phir>phig
    if vinfoutmag>=vinfinmag %Constrained RF Case
%      1   
     
% Get DV and DV/dta at rSOI (max true anomaly / max radius)

r1=rsoi;
tamax=acos((pin-r1)/(ein*r1));
ta1=tamax;
vin=sqrt(2*mu/r1+vinfinmag2);
vout=sqrt(2*mu/r1+vinfoutmag2);
fpain=acos(sqrt(mu*pin)/(r1*vin));
theta=(phir-acos(-1/ein)-ta1)+pi;
c1=mu/(r1*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv1=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));


ddv1=0.5*(-4*mu/r1^2+(8*mu^2/r1^3+2*mu*(vinfinmag2+vinfoutmag2)/r1^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r1+pin*c1^2*sin(theta)/(mu*r1*ein*sin(ta1)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r1^2*sin(ta1))+mu*cos(fpaout)/(vout^2*r1^2))/sin(fpaout)+(r1*vinfinmag2+mu)*cot(fpain)/(r1*vin)^2))/dv1;
%###################convert dDV/dr to dDV/dta###########################
ddv1dta=ddv1*r1*ein*sin(ta1)/(1+ein*cos(ta1));
%#######################################################################        
        if sign(ddv1dta)<0 %Check if boundary case
                dvd=dv1; %Boundary case solution
        else %interior case solution algorithm
            
f=Inf; %initialize f root zero-function
tata=phir+pi-acos(-1/ein)-acos(-1/eout);
H=aout/ain*(1-eout^2)/(1-ein^2)/eout;
ta0=tata/2; %initialize true anomaly
while (abs(f)>ftol)  %Determine minimum true anomaly boundary at rmin (rpmag)
    K=H*(1+ein*cos(ta0))-1/eout;
    f=ta0-tata+acos(K);
    df=1+H*ein*sin(ta0)/sqrt(1-K^2);
    ddf=(-K*(df-1)^2+H*ein*cos(ta0))/sqrt(1-K^2);
%        tamin=ta0-f/df;
    tamin=ta0-2*f*df/(2*df^2-f*ddf);
    
    if ~isreal(tamin)
        tamin=ta0;
        break
    end
    ta0=tamin;
end            


% rmin=pin/(1+ein*cos(tamin)); 
tamin=max([tamin tazero]); %adjust min ta to avoid singularity at ta near zero

% ta2=(tamin+ta)/2; %Compute next point between min and max true anomaly 
ta2=tamin; %Compute next point at min ta, bracketing the ta @ opt solution

r2=pin/(1+ein*cos(ta2)); 
vin=sqrt(2*mu/r2+vinfinmag2);
vout=sqrt(2*mu/r2+vinfoutmag2);
fpain=acos(sqrt(mu*pin)/(r2*vin));
theta=(phir-acos(-1/ein)-ta2)+pi;
c1=mu/(r2*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv2=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));

ddv2=0.5*(-4*mu/r2^2+(8*mu^2/r2^3+2*mu*(vinfinmag2+vinfoutmag2)/r2^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r2+pin*c1^2*sin(theta)/(mu*r2*ein*sin(ta2)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r2^2*sin(ta2))+mu*cos(fpaout)/(vout^2*r2^2))/sin(fpaout)+(r2*vinfinmag2+mu)*cot(fpain)/(r2*vin)^2))/dv2;

%###################convert dDV/dr to dDV/dta###########################
ddv2dta=ddv2*r2*ein*sin(ta2)/(1+ein*cos(ta2));
%#######################################################################

%Determine cubic interpolation estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

%Check if unacceptable solutions 
% if (ta3<tamin)
%       ta3=tamin; %replace ta estimate with lower bound ta
% elseif (ta3>tamax)
%       ta3=tamax; %replace ta estimate with upper bound ta
% end
% tacheck=abs(ta2-ta3);

r=pin/(1+ein*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
% ta3=acos((pin-r)/(ein*r));
fpain=acos(sqrt(mu*pin)/(r*vin));
theta=(phir-acos(-1/ein)-ta3)+pi;
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r+pin*c1^2*sin(theta)/(mu*r*ein*sin(ta3)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r^2*sin(ta3))+mu*cos(fpaout)/(vout^2*r^2))/sin(fpaout)+(r*vinfinmag2+mu)*cot(fpain)/(r*vin)^2))/dv3;
% dvcheck=abs(dv2-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*ein*sin(ta3)/(1+ein*cos(ta3));
%#######################################################################
tacheck=Inf;
dvcheck=Inf;

        while (tacheck>tatol)&&(dvcheck>dvtol)
%update point-set
if ddv3dta>0
ta1=ta3;
dv1=dv3;
ddv1dta=ddv3dta;
else
ta2=ta3;
dv2=dv3;
ddv2dta=ddv3dta;
end
taold=ta3;
dvold=dv3;

%Determine cubic interpolation next estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

%Check if unacceptable solutions 
% if (ta3<tamin) 
%       ta3=tamin; %replace ta estimate with lower bound ta
% elseif (ta3>tamax)
%       ta3=tamax; %replace ta estimate with upper bound ta
% end
tacheck=abs(taold-ta3);

r=pin/(1+ein*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
fpain=acos(sqrt(mu*pin)/(r*vin));
theta=(phir-acos(-1/ein)-ta3)+pi;
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r+pin*c1^2*sin(theta)/(mu*r*ein*sin(ta3)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r^2*sin(ta3))+mu*cos(fpaout)/(vout^2*r^2))/sin(fpaout)+(r*vinfinmag2+mu)*cot(fpain)/(r*vin)^2))/dv3;
dvcheck=abs(dvold-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*ein*sin(ta3)/(1+ein*cos(ta3));
%#######################################################################
iter=iter+1;
        end
    dvd=dv3;

        end

    elseif vinfoutmag<vinfinmag %Constrained RF Case
%     2        

% Get DV and DV/dta at rSOI (max true anomaly / max radius)
r1=rsoi;
tamax=-acos((pout-r1)/(eout*r1));
ta1=tamax;
vin=sqrt(2*mu/r1+vinfinmag2);
vout=sqrt(2*mu/r1+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r1*vout));
theta=(phir-acos(-1/eout)+ta1+pi);
c1=mu/(r1*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv1=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv1=0.5*(-4*mu/r1^2+(8*mu^2/r1^3+2*mu*(vinfinmag2+vinfoutmag2)/r1^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r1-pout*c1^2*sin(theta)/(mu*r1*eout*sin(ta1)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r1^2*sin(ta1))+mu*cos(fpain)/(vin^2*r1^2))/sin(fpain)+(r1*vinfoutmag2+mu)*cot(fpaout)/(r1*vout)^2))/dv1;
%###################convert dDV/dr to dDV/dta###########################
ddv1dta=ddv1*r1*eout*sin(ta1)/(1+eout*cos(ta1));
%#######################################################################         
        if sign(ddv1dta)>0 %Check if boundary case
                dvd=dv1; %Boundary case solution
        else %interior case solution algorithm
f=Inf; %initialize f root zero-function
tata=phir+pi-acos(-1/ein)-acos(-1/eout);
H=ain/aout*(1-ein^2)/(1-eout^2)/ein;
ta0=tata/2; %initialize true anomaly

while abs(f)>ftol %Determine minimum true anomaly boundary at rmin (rpmag)
    K=H*(1+eout*cos(ta0))-1/ein;
    f=ta0-tata+acos(K);
    df=1+H*eout*sin(ta0)/sqrt(1-K^2);
    ddf=(-K*(df-1)^2+H*eout*cos(ta0))/sqrt(1-K^2);
%        tamin=ta0-f/df;
    tamin=ta0-2*f*df/(2*df^2-f*ddf);
    
    if ~isreal(tamin)
        tamin=ta0;
        break
    end
    ta0=tamin;
end

% rmin=pout/(1+eout*cos(tamin)); 
tamin=min([-tamin -tazero]); %adjust min ta to avoid singularity at ta near zero

% ta2=(tamin+ta)/2 %Compute next point between min and max true anomaly 
ta2=tamin; %Compute next point at min ta, bracketing the ta @ opt solution
r2=pout/(1+eout*cos(ta2));
vin=sqrt(2*mu/r2+vinfinmag2);
vout=sqrt(2*mu/r2+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r2*vout));
theta=(phir-acos(-1/eout)+ta2+pi);
c1=mu/(r2*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv2=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv2=0.5*(-4*mu/r2^2+(8*mu^2/r2^3+2*mu*(vinfinmag2+vinfoutmag2)/r2^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r2-pout*c1^2*sin(theta)/(mu*r2*eout*sin(ta2)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r2^2*sin(ta2))+mu*cos(fpain)/(vin^2*r2^2))/sin(fpain)+(r2*vinfoutmag2+mu)*cot(fpaout)/(r2*vout)^2))/dv2;
%###################convert dDV/dr to dDV/dta###########################
ddv2dta=ddv2*r2*eout*sin(ta2)/(1+eout*cos(ta2));
%####################################################################### 

%Determine cubic interpolation estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

%Check if unacceptable solutions 
% if (ta>tamin)
%       ta=tamin; %replace ta estimate with lower bound ta
% elseif (ta<tamax)
%       ta=tamax; %replace ta estimate with upper bound ta
% end
% tacheck=abs(ta2-ta3);

r=pout/(1+eout*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r*vout));
theta=(phir-acos(-1/eout)+ta3+pi);
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r-pout*c1^2*sin(theta)/(mu*r*eout*sin(ta3)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r^2*sin(ta3))+mu*cos(fpain)/(vin^2*r^2))/sin(fpain)+(r*vinfoutmag2+mu)*cot(fpaout)/(r*vout)^2))/dv3;
% dvcheck=abs(dv2-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*eout*sin(ta3)/(1+eout*cos(ta3));
%#######################################################################        
tacheck=Inf;
dvcheck=Inf;       
       
        while (tacheck>tatol)&&(dvcheck>dvtol)
%update point-set
if ddv3dta>0            
ta2=ta3;
dv2=dv3;
ddv2dta=ddv3dta;
else
ta1=ta3;
dv1=dv3;
ddv1dta=ddv3dta;
end
taold=ta3;
dvold=dv3;

%Determine cubic interpolation next estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

%Check if unacceptable solutions 
% if (ta>tamin)
%       ta=tamin; %replace ta estimate with lower bound ta
% elseif (ta<tamax)
%       ta=tamax; %replace ta estimate with upper bound ta
% end
tacheck=abs(taold-ta3);

r=pout/(1+eout*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r*vout));
theta=(phir-acos(-1/eout)+ta3+pi);
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r-pout*c1^2*sin(theta)/(mu*r*eout*sin(ta3)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r^2*sin(ta3))+mu*cos(fpain)/(vin^2*r^2))/sin(fpain)+(r*vinfoutmag2+mu)*cot(fpaout)/(r*vout)^2))/dv3;
dvcheck=abs(dvold-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*eout*sin(ta3)/(1+eout*cos(ta3));
%#######################################################################
iter=iter+1;
        end

    dvd=dv3;
        end
    end
elseif phir<=phig
    if vinfoutmag>=vinfinmag %F case or Constrained RF Case
%     3     
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
sigma=(pi-phir)/4;
delta=atan((vinfinmag-vinfoutmag)/(vinfinmag+vinfoutmag)*tan(sigma));
F=(sigma-delta);
% G=(sigma+delta);
L=sqrt(2*mu/rpmag);
rfbmag=rpmag*L*L*sin(F)^2/(vinfoutmag2*cos(delta)^2*(2*cos(sigma)^2-cos(delta)^2));
dvd=abs((vinfinmag+vinfoutmag)*sin(delta));
rfbhat=rotuv(vinfinhat,hhat,-2*F);
costhetar=(rfbhat'*vinfouthat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rfbhat;
rfb=rfbhat*rfbmag;
evec=((vout'*vout)*rfb-(rfb'*vout)*vout)/mu-rfbhat;
rpfbmag=aout*(1-sqrt(evec'*evec));

if rpfbmag<rpmag || rpfbmag>rsoi || ~optr %check if RF Case or fixed R option
%     30
    %Optimized constrained RF Case
    tamax=-acos((pout-rsoi)/(eout*rsoi));
    tamin= -tazero;
    ta1=tamin; %arbitrary small angle near rp @ rmin
    ta2=-acos((aout*(1-evec'*evec)-rfbmag)/(sqrt(evec'*evec)*rfbmag)); %Compute next point using prior infeasible solution ta
  
   if (rpfbmag>rsoi)&&(optr)
       rpmag=rsoi; %adjust periapse to rsoi
       eout=1+rpmag*vinfoutmag2/mu;
       pout = aout*(1-eout*eout);
       tamax=-acos((pout-2*rsoi)/(eout*2*rsoi)); %tamax at 2*rsoi
       ta2=tamax;
   end
   
% Get DV and DV/dta @ ta1
r1=pout/(1+eout*cos(ta1));
vin=sqrt(2*mu/r1+vinfinmag2);
vout=sqrt(2*mu/r1+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r1*vout));
theta=(phir-acos(-1/eout)+ta1+pi);
c1=mu/(r1*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv1=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv1=0.5*(-4*mu/r1^2+(8*mu^2/r1^3+2*mu*(vinfinmag2+vinfoutmag2)/r1^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r1-pout*c1^2*sin(theta)/(mu*r1*eout*sin(ta1)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r1^2*sin(ta1))+mu*cos(fpain)/(vin^2*r1^2))/sin(fpain)+(r1*vinfoutmag2+mu)*cot(fpaout)/(r1*vout)^2))/dv1;
%###################convert dDV/dr to dDV/dta###########################
ddv1dta=ddv1*r1*eout*sin(ta1)/(1+eout*cos(ta1));
%#######################################################################         


% Get DV and DV/dta @ ta2
r2=pout/(1+eout*cos(ta2));
vin=sqrt(2*mu/r2+vinfinmag2);
vout=sqrt(2*mu/r2+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r2*vout));
theta=(phir-acos(-1/eout)+ta2+pi);
c1=mu/(r2*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv2=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv2=0.5*(-4*mu/r2^2+(8*mu^2/r2^3+2*mu*(vinfinmag2+vinfoutmag2)/r2^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r2-pout*c1^2*sin(theta)/(mu*r2*eout*sin(ta2)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r2^2*sin(ta2))+mu*cos(fpain)/(vin^2*r2^2))/sin(fpain)+(r2*vinfoutmag2+mu)*cot(fpaout)/(r2*vout)^2))/dv2;
%###################convert dDV/dr to dDV/dta###########################
ddv2dta=ddv2*r2*eout*sin(ta2)/(1+eout*cos(ta2));
%####################################################################### 

%Determine cubic interpolation estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

tacheck=abs(ta2-ta3);
%Check if unacceptable solutions 
if (ta3>tamin)||~isreal(ta3)||isnan(ta3)||(ddv1dta<0)
      ta3=tamin; %replace ta estimate with lower bound ta
      tacheck=1e10; %set to arbitrary large value
elseif (ta3<tamax)
      ta3=tamax; %replace ta estimate with upper bound ta
      tacheck=1e10; %set to arbitrary large value
end


r=pout/(1+eout*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r*vout));
theta=(phir-acos(-1/eout)+ta3+pi);
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r-pout*c1^2*sin(theta)/(mu*r*eout*sin(ta3)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r^2*sin(ta3))+mu*cos(fpain)/(vin^2*r^2))/sin(fpain)+(r*vinfoutmag2+mu)*cot(fpaout)/(r*vout)^2))/dv3;
dvcheck=abs(dv2-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*eout*sin(ta3)/(1+eout*cos(ta3));
%#######################################################################        
           
   
        while (tacheck>tatol)&&(dvcheck>dvtol)
%update point-set
if ddv3dta>0
ta1=ta3;
dv1=dv3;
ddv1dta=ddv3dta;
else
ta2=ta3;
dv2=dv3;
ddv2dta=ddv3dta;  
end
taold=ta3;
dvold=dv3;

%Determine cubic interpolation next estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

%Check if unacceptable solutions 
if (ta3>tamin)||~isreal(ta3)||isnan(ta3)
      ta3=tamin; %replace ta estimate with lower bound ta
elseif (ta3<tamax)
      ta3=tamax; %replace ta estimate with upper bound ta
end
tacheck=abs(taold-ta3);

r=pout/(1+eout*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
fpaout=-acos(sqrt(mu*pout)/(r*vout));
theta=(phir-acos(-1/eout)+ta3+pi);
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpain=-((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*cos(theta)+(stheta*sqrt(c1+vinfinmag2/4)-vinfinmag/2)); %determine sign for fpain
fpain=sign(sfpain)*acos((stheta*sqrt(c1+vinfinmag2/4)+vinfinmag/2)*sin(theta)/vin);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    +2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r-pout*c1^2*sin(theta)/(mu*r*eout*sin(ta3)))*sin(theta)/(vin*sqrt(4*c1+vinfinmag2))...
    +cos(fpain)*cot(theta)*pout/(eout*r^2*sin(ta3))+mu*cos(fpain)/(vin^2*r^2))/sin(fpain)+(r*vinfoutmag2+mu)*cot(fpaout)/(r*vout)^2))/dv3;
dvcheck=abs(dvold-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*eout*sin(ta3)/(1+eout*cos(ta3));
%#######################################################################
iter=iter+1;
        end
dvd=dv3;

end
    elseif vinfoutmag<vinfinmag %F case or Constrained RF Case
%     4
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
sigma=(pi-phir)/4;
delta=atan((vinfinmag-vinfoutmag)/(vinfinmag+vinfoutmag)*tan(sigma));
F=(sigma-delta);
% G=(sigma+delta);
L=sqrt(2*mu/rpmag);
rfbmag=rpmag*L*L*sin(F)^2/(vinfoutmag2*cos(delta)^2*(2*cos(sigma)^2-cos(delta)^2));
dvd=abs((vinfinmag+vinfoutmag)*sin(delta));
rfbhat=rotuv(vinfinhat,hhat,-2*F);
costhetar=(-rfbhat'*vinfinhat);
D=sqrt(mu/(rfbmag*(1+costhetar))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rfbhat;
rfb=rfbhat*rfbmag;
evec=((vin'*vin)*rfb-(rfb'*vin)*vin)/mu-rfbhat;
rpfbmag=ain*(1-sqrt(evec'*evec));

if rpfbmag<rpmag || rpfbmag>rsoi || ~optr %check if RF Case or fixed R option
%     40
    %Optimized constrained RF Case
    tamax=acos((pin-rsoi)/(ein*rsoi));
    tamin=tazero;
    ta1=tamin; %arbitrary small angle near rp @ rmin
    ta2=acos((ain*(1-evec'*evec)-rfbmag)/(sqrt(evec'*evec)*rfbmag)); %Compute next point using prior infeasible solution ta

    if (rpfbmag>rsoi)&&(optr)
       rpmag=rsoi; %adjust periapse to rsoi
       ein=1+rpmag*vinfinmag2/mu;
       pin = ain*(1-ein*ein);
       tamax=acos((pin-2*rsoi)/(ein*2*rsoi)); %tamax at 2*rsoi
    end
    
% Get DV and DV/dta @ ta1  
r1=pin/(1+ein*cos(ta1)); 
vin=sqrt(2*mu/r1+vinfinmag2);
vout=sqrt(2*mu/r1+vinfoutmag2);
fpain=acos(sqrt(mu*pin)/(r1*vin));
theta=(phir-acos(-1/ein)-ta1)+pi;
c1=mu/(r1*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv1=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv1=0.5*(-4*mu/r1^2+(8*mu^2/r1^3+2*mu*(vinfinmag2+vinfoutmag2)/r1^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r1+pin*c1^2*sin(theta)/(mu*r1*ein*sin(ta1)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r1^2*sin(ta1))+mu*cos(fpaout)/(vout^2*r1^2))/sin(fpaout)+(r1*vinfinmag2+mu)*cot(fpain)/(r1*vin)^2))/dv1;
%###################convert dDV/dr to dDV/dta###########################
ddv1dta=ddv1*r1*ein*sin(ta1)/(1+ein*cos(ta1));
%#######################################################################        
 
% Get DV and DV/dta @ ta2  
r2=pin/(1+ein*cos(ta2)); 
vin=sqrt(2*mu/r2+vinfinmag2);
vout=sqrt(2*mu/r2+vinfoutmag2);
fpain=acos(sqrt(mu*pin)/(r2*vin));
theta=(phir-acos(-1/ein)-ta2)+pi;
c1=mu/(r2*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv2=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv2=0.5*(-4*mu/r2^2+(8*mu^2/r2^3+2*mu*(vinfinmag2+vinfoutmag2)/r2^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r2+pin*c1^2*sin(theta)/(mu*r2*ein*sin(ta2)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r2^2*sin(ta2))+mu*cos(fpaout)/(vout^2*r2^2))/sin(fpaout)+(r2*vinfinmag2+mu)*cot(fpain)/(r2*vin)^2))/dv2;

%###################convert dDV/dr to dDV/dta###########################
ddv2dta=ddv2*r2*ein*sin(ta2)/(1+ein*cos(ta2));
%#######################################################################

%Determine cubic interpolation estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

%Check if unacceptable solutions 
if (ta3<tamin)||~isreal(ta3)||isnan(ta3)||(ddv1dta>0)
      ta3=tamin; %replace ta estimate with lower bound ta
elseif (ta3>tamax)
      ta3=tamax; %replace ta estimate with upper bound ta
end
tacheck=abs(ta2-ta3);

r=pin/(1+ein*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
% ta3=acos((pin-r)/(ein*r));
fpain=acos(sqrt(mu*pin)/(r*vin));
theta=(phir-acos(-1/ein)-ta3)+pi;
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r+pin*c1^2*sin(theta)/(mu*r*ein*sin(ta3)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r^2*sin(ta3))+mu*cos(fpaout)/(vout^2*r^2))/sin(fpaout)+(r*vinfinmag2+mu)*cot(fpain)/(r*vin)^2))/dv3;
dvcheck=abs(dv2-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*ein*sin(ta3)/(1+ein*cos(ta3));
%#######################################################################
        
    
        while (tacheck>tatol)&&(dvcheck>dvtol)
%update point-set
if ddv3dta>0
ta2=ta3;
dv2=dv3;
ddv2dta=ddv3dta;
else
ta1=ta3;
dv1=dv3;
ddv1dta=ddv3dta;  
end
taold=ta3;
dvold=dv3;

%Determine cubic interpolation next estimated ta @ dv_min  
dta=(ta2-ta1);
B1=ddv1dta+ddv2dta+3*(dv1-dv2)/(dta);
B2=sign(dta)*sqrt(B1*B1-ddv1dta*ddv2dta);
ta3=ta2-(dta)*(ddv2dta+B2-B1)/(ddv2dta-ddv1dta+2*B2);

%Check if unacceptable solutions 
if (ta3<tamin)||~isreal(ta3)||isnan(ta3)
      ta3=tamin; %replace ta estimate with lower bound ta
elseif (ta3>tamax)
      ta3=tamax; %replace ta estimate with upper bound ta
end
tacheck=abs(taold-ta3);

r=pin/(1+ein*cos(ta3));
vin=sqrt(2*mu/r+vinfinmag2);
vout=sqrt(2*mu/r+vinfoutmag2);
fpain=acos(sqrt(mu*pin)/(r*vin));
theta=(phir-acos(-1/ein)-ta3)+pi;
c1=mu/(r*(1+cos(theta)));
stheta=sign(pi-theta); %determine sign for D term in hyperbolic bvp
sfpaout=((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*cos(theta)+(stheta*sqrt(c1+vinfoutmag2/4)-vinfoutmag/2)); %determine sign for fpaout
fpaout=sign(sfpaout)*acos((stheta*sqrt(c1+vinfoutmag2/4)+vinfoutmag/2)*sin(theta)/vout);
dv3=sqrt(vin^2+vout^2-2*vin*vout*cos(fpaout-fpain));
ddv3=0.5*(-4*mu/r^2+(8*mu^2/r^3+2*mu*(vinfinmag2+vinfoutmag2)/r^2)*cos(fpaout-fpain)/(vin*vout)...
    -2*vin*vout*sin(fpaout-fpain)*((-stheta*(c1/r+pin*c1^2*sin(theta)/(mu*r*ein*sin(ta3)))*sin(theta)/(vout*sqrt(4*c1+vinfoutmag2))...
    -cos(fpaout)*cot(theta)*pin/(ein*r^2*sin(ta3))+mu*cos(fpaout)/(vout^2*r^2))/sin(fpaout)+(r*vinfinmag2+mu)*cot(fpain)/(r*vin)^2))/dv3;
dvcheck=abs(dvold-dv3);
%###################convert dDV/dr to dDV/dta###########################
ddv3dta=ddv3*r*ein*sin(ta3)/(1+ein*cos(ta3));
%#######################################################################
iter=iter+1;
        end
    dvd=dv3;  
    
end
    end
end

end
