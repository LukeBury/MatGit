function [x0,rnb,nmcm,nsb]=mcmguess(type,event,JD0,tof,nplanet,gradet,rev,revig,...
    coenp,nodes,trajs,x0,rnb,musun)
Rj=71492;   %km

JD(1)=JD0;
for i=1:trajs
JD(i+1)=JD(i)+tof(i); 
%Determine julian dates from initial JD and times of flight  
end

%Check & Determine Midcourse Maneuvers Initial Guess
n=2; %initialize counter
while n<=nodes
   nmcm=0;
if event(n)==40
      while event(n)==40
         n=n+1;
         nmcm=nmcm+1; 
         if n>=nodes, break, end
      end
n1=n-nmcm-1;
n2=n;

%sum all trajectory segment's mulit-revolutions 
revs=sum([rev(n1+1:n2-1),revig(n1+1:n2-1)]); 

%Collect ephemeris and state data for other nodes
[rP1,vP1]=ephstate(type,event(n1),nplanet(n1),JD(n1),coenp(:,n1),rnb(:,n1));
[rP2,vP2]=ephstate(type,event(n2),nplanet(n2),JD(n2),coenp(:,n2),rnb(:,n2));

%Hybrid lambert targeter with Dr. Izzo (ESA) fast convergence solution and 
%Lancaster & Blancard's robust solution 
 [vd,va,exitflag]=lambert(rP1,rP2,(JD(n2)-JD(n1))*86400,revs,mode(gradet(n1+1:n2-1)),musun);
 if (exitflag ~=1)||any(isnan(vd))||any(isnan(va))
     if exitflag==1 %Catch any NaN cases that exitflag didn't
     exitflag=-1;
     end
     disp(['!! Lambert Targeter failed between planet ',num2str(nplanet(n1)),' and planet ',num2str(nplanet(n2)),'!!'])
 end

hhat=cross(rP1,rP2)/norm(cross(rP1,rP2));
theta = acos(hhat(3));
psi = atan2(hhat(2),hhat(1));
rot2=[cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)]; 
%frame rotation about y-axis
rot3 = [cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0; 0 0 1]; 
%frame rotation about z-axis
rot2n=[cos(-theta) 0 -sin(-theta);0 1 0; sin(-theta) 0 cos(-theta)]; 
%negative frame rotation about y-axis
rot3n = [cos(-psi) sin(-psi) 0;-sin(-psi) cos(-psi) 0; 0 0 1]; 
%negative frame rotation about z-axis
rP1=rot2*rot3*rP1;
rP2=rot2*rot3*rP2;
phi1=atan2(rP1(2),rP1(1));
phi2=atan2(rP2(2),rP2(1));

if hhat(3)<0
    if mode(gradet(n1+1:n2-1))==1
        if (phi2-phi1)>0
        phi=phi1+((phi2-phi1)+2*pi*revs)*(1:nmcm)/(nmcm+1);
        else
        phi=phi1+(2*pi+(phi2-phi1)+2*pi*revs)*(1:nmcm)/(nmcm+1);
        end
    else
        if (phi2-phi1)<0
        phi=phi1+((phi2-phi1)-2*pi*revs)*(1:nmcm)/(nmcm+1);
        else
        phi=phi1-(2*pi-(phi2-phi1)-2*pi*revs)*(1:nmcm)/(nmcm+1);
        end
    end
else
    if mode(gradet(n1+1:n2-1))==1
        if (phi2-phi1)<0
        phi=phi1+((phi2-phi1)-2*pi*revs)*(1:nmcm)/(nmcm+1);
        else
        phi=phi1+(-2*pi+(phi2-phi1)-2*pi*revs)*(1:nmcm)/(nmcm+1);  
        end
    else
        if (phi2-phi1)>0
        phi=phi1+((phi2-phi1)+2*pi*revs)*(1:nmcm)/(nmcm+1);
        else
        phi=phi1+(2*pi+(phi2-phi1)+2*pi*revs)*(1:nmcm)/(nmcm+1);
        end
    end
end
rP1mag=norm(rP1);
    for i=1:nmcm
rmcm(1:3,i)= (rP1mag+(norm(rP2)-rP1mag)*(i)/(nmcm+1))/Rj*rot3n*rot2n*...
    [cos(phi(i)); sin(phi(i)); 0];
rnb(1:3,n1+i)=rmcm(1:3,i);
    end
rmcm=reshape(rmcm,1,nmcm*3);
x0=[x0 rmcm];
rmcm=[];
end
n=n+1;
end

nmcm=(length(x0)-nodes)/3;
nsb=0;  %initialize spaceburn counter
for n=1:nodes
    if event(n)==41
        x0=[x0 rnb(1:3,n)'];
        nsb=nsb+1;
    end
end
end