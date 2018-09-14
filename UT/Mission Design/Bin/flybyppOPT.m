function [dvd]=flybyppOPT(rpmag,rsoi,rP,vinfinhat,vinfinmag,vinfouthat,vinfoutmag,...
    mu,optr)
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
% optr: optimization of rp
%------OUTPUTS-----
% dvd: Delta-V required for hyperbolic depature trajectory
% 
zero=1e-12;
vinfinmag2=vinfinmag*vinfinmag;
vinfoutmag2=vinfoutmag*vinfoutmag;
ain = -mu/vinfinmag2;    %semi-transverse axis of hyperbolic trajectory
aout= -mu/vinfoutmag2;
ein=1-rpmag/ain;
eout=1-rpmag/aout;
phig=asin(1/eout)+asin(1/ein); %maximum natural gravitational turn angle 
phir=acos(vinfouthat'*vinfinhat); %dot prod

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
       end
   end
end
    
    
    
if phir>phig
    if vinfoutmag>=vinfinmag %Constrained RF Case
%      1   
% rpinhat = rotuv(vinfinhat,hhat,-acos(1/ein));
nuinf=acos(-1/ein);
rpinhat= [vinfinhat(2)*hhat(3) - vinfinhat(3)*hhat(2);...%cross product 
        vinfinhat(3)*hhat(1) - vinfinhat(1)*hhat(3);...
        vinfinhat(1)*hhat(2) - vinfinhat(2)*hhat(1)];
rpinhat=-cos(nuinf)*vinfinhat+sin(nuinf)*rpinhat;
rfb = rpmag*rpinhat;
% vpinhat=rotuv(vinfinhat,hhat,pi/2-betain);
vpinhat=[hhat(2)*rpinhat(3) - hhat(3)*rpinhat(2);...%cross product 
          hhat(3)*rpinhat(1) - hhat(1)*rpinhat(3);...
          hhat(1)*rpinhat(2) - hhat(2)*rpinhat(1)];
vin=vpinhat*sqrt(2*mu/rpmag+vinfinmag2);
D=sqrt(mu/(rpmag*(1+rpinhat'*vinfouthat))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat;  
h3=rfb(1)*vout(2)-rfb(2)*vout(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat; 
end

    elseif vinfoutmag<vinfinmag %Constrained RF Case
%     2        
% rpouthat = rotuv(-vinfouthat,hhat,acos(1/eout));
nuinf=acos(-1/eout);
rpouthat= [vinfouthat(2)*hhat(3) - vinfouthat(3)*hhat(2);...%cross product 
        vinfouthat(3)*hhat(1) - vinfouthat(1)*hhat(3);...
        vinfouthat(1)*hhat(2) - vinfouthat(2)*hhat(1)];
rpouthat=cos(nuinf)*vinfouthat+sin(nuinf)*rpouthat;
rfb = rpmag*rpouthat;
% vpouthat=rotuv(vinfouthat,hhat,betaout-pi/2);
vpouthat=[hhat(2)*rpouthat(3) - hhat(3)*rpouthat(2);...%cross product 
          hhat(3)*rpouthat(1) - hhat(1)*rpouthat(3);...
          hhat(1)*rpouthat(2) - hhat(2)*rpouthat(1)];
vout=vpouthat*sqrt(2*mu/rpmag+vinfoutmag2);
D=sqrt(mu/(rpmag*(1-rpouthat'*vinfinhat))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
h3=rfb(1)*vin(2)-rfb(2)*vin(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
end
    
    end
elseif phir<=phig
    if vinfoutmag>=vinfinmag %F case or Constrained RF Case
%     3     
% rpouthat = rotuv(-vinfouthat,hhat,acos(1/eout));
nuinf=acos(-1/eout);
rpouthat= [vinfouthat(2)*hhat(3) - vinfouthat(3)*hhat(2);...%cross product 
        vinfouthat(3)*hhat(1) - vinfouthat(1)*hhat(3);...
        vinfouthat(1)*hhat(2) - vinfouthat(2)*hhat(1)];
rpouthat=cos(nuinf)*vinfouthat+sin(nuinf)*rpouthat;
rfb = rpmag*rpouthat;
% vpouthat=rotuv(vinfouthat,hhat,betaout-pi/2);
vpouthat=[hhat(2)*rpouthat(3) - hhat(3)*rpouthat(2);...%cross product 
          hhat(3)*rpouthat(1) - hhat(1)*rpouthat(3);...
          hhat(1)*rpouthat(2) - hhat(2)*rpouthat(1)];
vout=vpouthat*sqrt(2*mu/rpmag+vinfoutmag2);
D=sqrt(mu/(rpmag*(1-rpouthat'*vinfinhat))+vinfinmag2/4);
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
h3=rfb(1)*vin(2)-rfb(2)*vin(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vin=(D+0.5*vinfinmag)*vinfinhat-(D-0.5*vinfinmag)*rpouthat;    
end
    elseif vinfoutmag<vinfinmag %F case or Constrained RF Case
%     4
% rpinhat = rotuv(vinfinhat,hhat,-acos(1/ein));
nuinf=acos(-1/ein);
rpinhat= [vinfinhat(2)*hhat(3) - vinfinhat(3)*hhat(2);...%cross product 
        vinfinhat(3)*hhat(1) - vinfinhat(1)*hhat(3);...
        vinfinhat(1)*hhat(2) - vinfinhat(2)*hhat(1)];
rpinhat=-cos(nuinf)*vinfinhat+sin(nuinf)*rpinhat;
rfb = rpmag*rpinhat;
% vpinhat=rotuv(vinfinhat,hhat,pi/2-betain);
vpinhat=[hhat(2)*rpinhat(3) - hhat(3)*rpinhat(2);...%cross product 
          hhat(3)*rpinhat(1) - hhat(1)*rpinhat(3);...
          hhat(1)*rpinhat(2) - hhat(2)*rpinhat(1)];
vin=vpinhat*sqrt(2*mu/rpmag+vinfinmag2);
D=sqrt(mu/(rpmag*(1+rpinhat'*vinfouthat))+vinfoutmag2/4);
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat;  
h3=rfb(1)*vout(2)-rfb(2)*vout(1);
if sign(hhat(3))~=sign(h3)
D=-D;
vout=(D+0.5*vinfoutmag)*vinfouthat+(D-0.5*vinfoutmag)*rpinhat; 
end
    end
end
dvd = sqrt((vout - vin)'*(vout - vin));
end