function [c,ceq] = nonlconfunc(x,event,nplanet,grade,rev,rfbc,coep,coenp,...
    JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,dvopt,optnode,...
    ephtype,Rj,musun,eps,JD0scale,tofscale,MCMscale,MCMRUB,MCMRLB,MCMREQ,DVUB,DVLB,DVEQ,vinfinUB,vinfinLB,...
	vinfinEQ,vinfoutUB,vinfoutLB,vinfoutEQ,vinfmatch,maxstepJD0,maxsteptof,maxstepMCM,tofmin)
%fmincon nonlinear constraint function
%   function structure mirrors that of trajcon primary algorithm, such that
%   the exact computed parameters for the optimization can be used for
%   constraint evaluation.
global xint
c=[];
ceq=[];
    x=x(:)'; %check for x orientation from fminsearchcon
    %Check for optimizer stepsize exceeding any specified lengths
if maxstepJD0~=0
    dx(1)=x(1)-xint(1);
    overstep(1)=(abs(dx(1))>maxstepJD0/JD0scale);
    if overstep(1)
       x(1)=xint(1)+sign(dx(1))*maxstepJD0/JD0scale; 
    end
end
if any(maxsteptof)~=0
    dx(2:nodes)=x(2:nodes)-xint(2:nodes);
    overstep(2:nodes)=(abs(dx(2:nodes))>=maxsteptof/tofscale)&(maxsteptof~=0);
    indtof=find(overstep(2:nodes)); 
    if ~isempty(indtof)
    x(indtof+1)=xint(indtof+1)+sign(dx(indtof+1)).*maxsteptof(indtof)./tofscale(indtof);
    end
end
nx=1;
if any(maxstepMCM)~=0
    if nmcm||nsb
    dx=x(nodes+1:end)-xint(nodes+1:end);
        for n=1:nodes
            if (event(n)==40)||(event(n)==41)
    overstep(nx:nx+2)=(abs(dx(nx:nx+2))>...
        abs(maxstepMCM(:,n)./MCMscale(:,n)))&norm(maxstepMCM(:,n));
    indMCM=find(overstep(nx:nx+2)); 
        if ~isempty(indMCM)
            ndx=nx+indMCM-1;
    x(nodes+ndx)=xint(nodes+ndx)+sign(dx(ndx)).*maxstepMCM(ndx,n)./abs(MCMscale(ndx,n));
    nx=nx+3;
        end
            end
        end
    end
end
% x;
 xint=x; %reset x0 vector for next optimizer iteration

rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end 

nc=0;   %initialize inequality constraint counter
nceq=0; %initialize equality constraint counter

if ~isempty(MCMRUB)||~isempty(MCMRLB)||~isempty(MCMREQ)   
    if ~isempty(MCMRUB)
        for ni=1:length(MCMRUB)
            if MCMRUB(ni)
       nc=nc+1;
       c(nc)=norm(rnb(:,ni))-MCMRUB(ni);
            end
        end
    end
    if ~isempty(MCMRLB)
        for ni=1:length(MCMRLB)
            if MCMRLB(ni)
       nc=nc+1;
       c(nc)=MCMRLB(ni)-norm(rnb(:,ni));
            end
        end
    end
    if ~isempty(MCMREQ)
        for ni=1:length(MCMREQ)
            if MCMREQ(ni)
       nceq=nceq+1;
       ceq(nceq)=MCMREQ(ni)-norm(rnb(:,ni));
            end
        end
    end 
end

if ~isempty(vinfinUB)||~isempty(vinfinLB)||~isempty(vinfinEQ)||...
   ~isempty(vinfoutUB)||~isempty(vinfoutLB)||~isempty(vinfoutEQ)||...
   ~isempty(vinfmatch)||~isempty(DVUB)||~isempty(DVLB)||~isempty(DVEQ)
    
JD=zeros(nodes,1);            %initialize/preallocate vector for efficiency
rP=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vd=zeros(3,nodes);          %initialize/preallocate vector for efficiency
va=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vinfout=zeros(3,nodes);     %initialize/preallocate vector for efficiency
vinfin=zeros(3,nodes);      %initialize/preallocate vector for efficiency
vinfouthat=zeros(3,nodes);  %initialize/preallocate vector for efficiency
vinfinhat=zeros(3,nodes);   %initialize/preallocate vector for efficiency
vinfoutmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
vinfinmag=zeros(nodes,1);     %initialize/preallocate vector for efficiency

         
JD(1)=JD0+x(1)*JD0scale;

%Compute first node state
[rP(:,1),vP(:,1)]=ephstate(ephtype,event(1),nplanet(1),JD(1),coenp(:,1),rnb(:,1));

for n=1:trajs
    np1=n+1;
    if abs(x(np1))*tofscale(n)<(tofmin*0.99) %check if tof = 0 and add small perturbation 
        disp(['tof ',num2str(n),' is near-zero and is set to 0.1 days'])
        x(np1)=tofmin/tofscale(n); %arbitrarily small perturbation (days)
    end
 JD(np1)=JD(n)+abs(x(np1))*tofscale(n); %Determine julian dates from initial JD and times of flight 
 
    %Collect ephemeris and state data for other nodes
[rP(:,np1),vP(:,np1)]=ephstate(ephtype,event(np1),nplanet(np1),JD(np1),coenp(:,np1),rnb(:,np1));

%Hybrid lambert targeter with Dr. Izzo (ESA) fast convergence solution and 
%Lancaster & Blancard's robust solution 
 [vd(:,n),va(:,np1),exitflag]=lambert(rP(:,n),rP(:,np1),(JD(np1)-JD(n))*86400,rev(n),grade(n),musun);
 if (exitflag ~=1)||any(isnan(vd(:,n)))||any(isnan(va(:,np1)))
     if exitflag==1 %Catch any NaN cases that exitflag didn't
     exitflag=-1;
     end
     disp(['!! Lambert Targeter failed between planet ',num2str(nplanet(n)),' and planet ',num2str(nplanet(np1)),'!!'])
 end
 
nex=1; %initialize exitflag loop counter
while exitflag<1
    if abs(rev(n))>0
        rev(n)=sign(rev(n))*(abs(rev(n))-1);
        disp(['!! Recomputing arc for leg ',num2str(n),' with ',num2str(rev(n)),' revs !!'])
    else
        disp(['!! Recomputing arc for leg ',num2str(n),' with tof increased by 0.1 days !!'])
        x(np1)=abs(x(np1))+0.1/tofscale(n); %small perturbation (days)
        JD(np1)=JD(n)+x(np1); %adjust JD to new tof
    end
    
    %recompute lambert arcs with adjusted parameters
    %Collect ephemeris and state data for other nodes
[rP(:,np1),vP(:,np1)]=ephstate(ephtype,event(np1),nplanet(np1),JD(np1),coenp(:,np1),rnb(:,np1));

%Hybrid lambert targeter with Dr. Izzo (ESA) fast convergence solution and 
%Lancaster & Blancard's robust solution 
 [vd(:,n),va(:,np1),exitflag]=lambert(rP(:,n),rP(:,np1),(JD(np1)-JD(n))*86400,rev(n),grade(n),musun);
 if (exitflag ~=1)||any(isnan(vd(:,n)))||any(isnan(va(:,np1)))
     if exitflag==1 %Catch any NaN cases that exitflag didn't
     exitflag=-1;
     end
     disp(['!! Lambert Targeter failed between planet ',num2str(nplanet(n)),' and planet ',num2str(nplanet(np1)),'!!'])
 end
nex=nex+1;

    if exitflag<1 && nex>=10
        disp(['!! LAMBERT ARC FOR LEG ',num2str(n),' HAS NO SOLUTION FOR ',num2str(rev(n)),' REVS AND ',num2str(x(np1)*tofscale(n)),' days TOF !!'])
        break
    end
end

end
%Determine hyperbolic excess velocity (v-infinity) at each node

for n=1:trajs
    np1=n+1;
vinfout(:,n) = vd(:,n) - vP(:,n);
vinfin(:,np1) = va(:,np1) - vP(:,np1);
vinfoutmag(n) = norm(vinfout(:,n));
vinfinmag(np1) = norm(vinfin(:,np1));
if (vinfoutmag(n)==0) %check for vinf=0 case
vinfouthat(:,n)=vP(:,n)/sqrt(vP(:,n)'*vP(:,n));
else
vinfouthat(:,n) = vinfout(:,n)/vinfoutmag(n);
end
if (vinfinmag(np1)==0) %check for vinf=0 case
vinfinhat(:,np1)=vP(:,np1)/sqrt(vP(:,np1)'*vP(:,np1));
else
vinfinhat(:,np1) = vinfin(:,np1)/vinfinmag(np1);
end
end

rsoi=zeros(nodes,1);          %initialize/preallocate vector for efficiency
[mupl,rpl,rap,dp,wp]=planetinfo(nplanet,JD,nodes); %Acquire planetary data 
rpmag=rfbc.*rpl;
%Determine sphere of influence and planetocentric periapsis radii
for n=1:nodes
rsoi(n) = (mupl(n)/musun)^(2/5)*norm(rP(:,n));      
end


if ~isempty(vinfinUB)
    for ni=1:length(vinfinUB)
        if vinfinUB(ni)
    nc=nc+1;
    c(nc)=vinfinmag(ni)-vinfinUB(ni);
        end
    end
end

if ~isempty(vinfinLB)
    for ni=1:length(vinfinLB)
        if vinfinLB(ni)
    nc=nc+1;
    c(nc)=vinfinLB(ni)-vinfinmag(ni);
        end
    end
end

if ~isempty(vinfinEQ)
    for ni=1:length(vinfinEQ)
        if vinfinEQ(ni)
    nceq=nceq+1;
    ceq(nceq)=vinfinEQ(ni)-vinfinmag(ni);
        end
    end
end

if ~isempty(vinfoutUB)
    for ni=1:length(vinfoutUB)
        if vinfoutUB(ni)
    nc=nc+1;
    c(nc)=vinfoutmag(ni)-vinfoutUB(ni);
        end
    end
end

if ~isempty(vinfoutLB)
    for ni=1:length(vinfoutLB)
        if vinfoutLB(ni)
    nc=nc+1;
    c(nc)=vinfoutLB(ni)-vinfoutmag(ni);
        end
    end
end

if ~isempty(vinfoutEQ)
    for ni=1:length(vinfoutEQ)
        if vinfoutEQ(ni)
    nceq=nceq+1;
    ceq(nceq)=vinfoutEQ(ni)-vinfoutmag(ni);
        end
    end
end

if ~isempty(vinfmatch)
    for ni=1:length(vinfmatch)
        if vinfmatch(ni)
    nceq=nceq+1;
    ceq(nceq)=vinfoutmag(ni)-vinfinmag(ni); %match vinf magnitudes
    
    nc=nc+1;  % determine if required turning angle is greater than max turning angle
    c(nc)=acos(vinfouthat(:,ni)'*vinfinhat(:,ni))-2*asin(1/(1+rpmag(ni)*vinfinmag(ni)*vinfoutmag(ni)/mupl(ni))); 
        end
    end
end

if ~isempty(DVUB)||~isempty(DVLB)||~isempty(DVEQ)  

dv=zeros(3,nodes);          %initialize/preallocate vector for efficiency

f=0;
for n=1:nodes;

%--------------------------Departure Trajectory----------------------------
    if event(n) == 10 %parking orbit with variable i

%     [dv(n)]= park2hypoptOPT(vinfoutmag(n),coep(1:2,n),mupl(n));
    vcmag=(mupl(n)/(coep(1,n)*(1-coep(2,n)))); %velocity^2 of circular orbit with same periapsis
    dv(n)=abs(sqrt(2*vcmag+vinfoutmag(n)*vinfoutmag(n))-sqrt(vcmag*(1+coep(2,n))));
       
    f=f+dv(n)*optnode(n);
       
        
    elseif event(n) == 11 %parking orbit with fixed i
 
     [dv(n)]= park2hypdOPT(vinfouthat(:,n),vinfoutmag(n),rap(n),dp(n),coep(1:3,n),mupl(n));    
        
     f=f+dv(n)*optnode(n); 
     
     
    elseif event(n) == 12 %parking orbit with variable true anomaly
     
        
        
    elseif event(n) == 13 %fixed parking orbit
     
     [dv(n)]= park2hypfixOPT(vinfouthat(:,n),vinfoutmag(n),rap(n),dp(n),coep(1:6,n),mupl(n));    
        
     f=f+dv(n)*optnode(n); 
        
    elseif event(n) == 14 %launch vehicle trajectory
 
     [dv(n)]= vinf2hypdOPT(vinfouthat(:,n),vinfoutmag(n),rpl(n),rap(n),dp(n),...
         vinfavmag(n),ilv(n),altp(n),mupl(n),dvopt(n));
        
     f=f+dv(n)*optnode(n); 
        
            
    elseif (event(n) == 17) %Redezvous with User-Defined Body Trajectory
    
    dv(n)=sqrt((vd(:,n)-vP(:,n))'*(vd(:,n)-vP(:,n)));
        
    f=f+dv(n)*optnode(n);        
  
    
    elseif (event(n) == 18) % Orbit about User-Defined Body 
        
    f=f+dv(n)*optnode(n);         
         
        
%---------------------------Arrival Trajectory-----------------------------        
    elseif event(n) == 20 %parking orbit with variable i
    
%     [dv(n)]= hyp2parkoptOPT(vinfinmag(n),coep(1:2,n),mupl(n)); 
    vcmag=(mupl(n)/(coep(1,n)*(1-coep(2,n)))); %velocity^2 of circular orbit with same periapsis 
    dv(n)=abs(sqrt(2*vcmag+vinfinmag(n)*vinfinmag(n))-sqrt(vcmag*(1+coep(2,n))));
    
    f=f+dv(n)*optnode(n); 
        
        
    elseif event(n) == 21 %parking orbit with fixed i
   
    [dv(n)]= hypa2parkOPT(vinfinhat(:,n),vinfinmag(n),rap(n),dp(n),coep(1:3,n),mupl(n)); 
        
    f=f+dv(n)*optnode(n);
        
    elseif event(n) == 22 %parking orbit with variable true anomaly
     
        
        
    elseif event(n) == 23 %fixed parking orbit
     
     [dv(n)]=  hyp2parkfixOPT(vinfinhat(:,n),vinfinmag(n),rap(n),dp(n),coep(1:6,n),mupl(n)); 
        
     f=f+dv(n)*optnode(n);
     
    
    elseif event(n) == 24 %entry trajectory
        
     [dv(n)]= hypa2vinfOPT(vinfinhat(:,n),vinfinmag(n),rpl(n),rap(n),dp(n),...
         vinfavmag(n),ilv(n),altp(n),mupl(n),dvopt(n));        
        
    f=f+dv(n)*optnode(n); 
        
 
    elseif (event(n) == 27) %Redezvous with User-Defined Body Trajectory
    
    dv(n)=sqrt((vP(:,n)-va(:,n))'*(vP(:,n)-va(:,n)));      

    f=f+dv(n)*optnode(n);      
        
        
    elseif (event(n) == 28) % Orbit about User-Defined Body
        
        
    f=f+dv(n)*optnode(n);       
          
%----------------------------Flyby Trajectory------------------------------        
    elseif event(n) == 30 %unpowered gravity assist
    
        dv(n)=0; %assume no deltaV
%   f=f;
        
    elseif event(n) == 31 %unpowered gravity assist w/ dv @ Rinf  
        
    [dv(n)]=flybyinfOPT(rpmag(n),rsoi(n),vinfinhat(:,n),vinfinmag(n),vinfouthat(:,n),...
    vinfoutmag(n),mupl(n),optr(n));   

    f=f+dv(n)*optnode(n); 
        
    
    elseif event(n) == 33 %optimal powered gravity assist
    
    [dv(n)]=flybypOPT(rpmag(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),...
    vinfouthat(:,n),vinfoutmag(n),mupl(n),rsoi(n),optr(n));
    
    f=f+dv(n)*optnode(n); 
        
        
    elseif event(n) == 35 %periapse powered gravity assist (residual delta-V at periapsis)
    
    [dv(n)]=flybyppOPT(rpmag(n),rsoi(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),...
    vinfouthat(:,n),vinfoutmag(n),mupl(n),optr(n));  
    
    f=f+dv(n)*optnode(n); 
        
%----------------------------Midcourse Maneuver Trajectory--------------------------    
    elseif (event(n) == 40)||(event(n) == 41)   
        
    dv(n)=sqrt((vd(:,n)-va(:,n))'*(vd(:,n)-va(:,n))); 

    f=f+dv(n)*optnode(n);
       
    end

end
        if ~isempty(DVUB)
            for ni=1:length(DVUB)
                if DVUB(ni)
            nc=nc+1;
            c(nc)=dv(ni)-DVUB(ni);
                end
            end
        end
        
        if ~isempty(DVLB)
            for ni=1:length(DVLB)
                if DVLB(ni)
            nc=nc+1;
            c(nc)=DVLB(ni)-dv(ni); 
                end
            end
        end
        
        if ~isempty(DVEQ)
            for ni=1:length(DVEQ)
                if DVEQ(ni)
            nceq=nceq+1;
            ceq(nceq)=DVEQ(ni)-dv(ni);  
                end
            end
        end
end

end