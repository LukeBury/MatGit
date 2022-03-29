%Patched Conic Trajectory Configuration Tool
function [dv,dve,JD,tof,rsoi,vinfout,vinfin,vinfpl,vinfoutmag,vinfinmag,...
    rP,vP,vPplhat,vd,va,vplhat,coet,rev,dvRA,dvDEC,rm,vout,...
    vin,rpout,rpin,vpout,vpin,rpfb,rpmag,betaout,betain,deltaout,deltain,tsoiout,...
    tsoiin,Tp,mupl,rpl,rap,dp,wp,coep,coenp,rnb,dvmag,rpfbmag,rfbc,rpfbmagr,output]=...
    trajcon(event,optimization,plotting,JD0,tof,nplanet,grade,rev,revig,branch,rfbc,...
    coep,coenp,rnb0,vinfavmag,ilv,altp,optr,optoptions,optnode,dvopt,optscale,...
    branchopt,maxstepJD0,maxsteptof,maxstepMCM,...
    tofTOTUB,tofTOTLB,tofTOTEQ,tofUB,tofLB,tofEQ,JD0UB,JD0LB,JDUB,JDLB,JDEQ,...
    MCMUB,MCMLB,MCMEQ,MCMRUB,MCMRLB,MCMREQ,DVUB,DVLB,DVEQ,vinfinUB,vinfinLB,...
	vinfinEQ,vinfoutUB,vinfoutLB,vinfoutEQ,vinfmatch,ephtype,savefile,loadfile,...
    excelfile,writeout,pts,dt,movie,fps)
% warning off
global xint

if ~isempty(loadfile)
   load(loadfile,'JD0','tof','rnb0')
end

strictness=[]; % default is loose, choose 'strict' or 'superscrict' :fminsearchcon
algorithm=[]; %default is internal Nealder-Mead, choose 'fminsearch'

nptsP=1000; %number of plotting points for Planets
nptsUD=1000; %number of plotting points for user defined bodies/orbits
nptsSC=500; %number of plotting points for spacecraft in heliocentric space
nptsSCa=500; %number of plotting points for spacecraft at planet arrival 
nptsparka=1000; %number of plotting points for arrival parking orbit 
nptsSCd=500; %number of plotting points for spacecraft at planet departure
nptsparkd=1000; %number of plotting points for departure parking orbit 
nptsflyby=500; %number of plotting points for gravity assist legs

nodes=length(event);
trajs=length(tof);
dxtol=1e-7;
tofmin=0.1; %set minimum time of flight to avoid undefined lambert arcs

muj = 1.26712764e+008;      %km3/s2  
Rj=71492;   %km
ecl= 0; % Obliquity of ecliptic at J2000

tof=abs(tof); %ensure positive tof inputs
x0=[0 tof] ; %initial date and times of flight wrt reference julian date

%--------------Check & Determine Midcourse Maneuvers Initial Guess---------

rnb=rnb0; %initialize maneuver location
[x0,rnb,nmcm,nsb]=mcmguess(ephtype,event,JD0,tof,nplanet,grade,rev,revig,coenp,nodes,trajs,x0,rnb,muj);
lenx0=length(x0);

    %Scale optimization variables
if optscale==0
    JD0scale=1; %days (no scaling)
    tofscale=ones(1,length(tof)); %days (no scaling)
    MCMscale=ones(size(rnb)); %AUs (no scaling)
%     x0=x0; %(no scaling)
elseif optscale==1
    JD0scale=tof(1); %days (scale by initial tof(1) value)
    tofscale=tof; %days (scale by initial tofs values)
    MCMscale=rnb; %AUs (scale by intiial MCM values)
%     optoptions.FinDiffRelStep=sqrt(eps)./[JD0scale x0(2:end)];   
    x0=[0 ones(1,lenx0-1)];
elseif optscale==2
    JD0scale=sum(tof); %days (scale by total tof value)
    tofscale=ones(1,length(tof)).*sum(tof); %days (scale by total tof value)
    MCMscale=rnb; %AUs (scale by intiial MCM values)
%     optoptions.FinDiffRelStep=sqrt(eps)./[ones(1,nodes).*sum(tof) x0(nodes+1:end)];   
    x0(2:nodes)=x0(2:nodes)./tofscale;
    x0(nodes+1:lenx0)=ones(1,lenx0-nodes); 
elseif optscale==3
    JD0scale=mean(tof); %days (scale by total tof value)
    tofscale=ones(1,length(tof)).*mean(tof); %days (scale by total tof value)
    MCMscale=rnb; %AUs (scale by intiial MCM values)
%     optoptions.FinDiffRelStep=sqrt(eps)./[ones(1,nodes).*mean(tof) x0(nodes+1:end)];   
    x0(2:nodes)=x0(2:nodes)./tofscale;
    x0(nodes+1:lenx0)=ones(1,lenx0-nodes); 
end

if optimization>=3
    Acon=[];    %initialize inequality matrix
    bcon=[];    %initialize inequality vector
    Aeq=[];     %initialize equality matrix
    beq=[];     %initialize equality vector
    lb=[];      %initialize lower bound vector
    ub=[];      %initialize upper bound vector
    nonlcon=[]; %initialize nonlinear constraint function

    %Apply positive TOF constraints
    lb(1:lenx0)=-Inf; %Initialize lb vector
    ub(1:lenx0)=Inf;  %Initialize ub vector
    %Apply positive time of flight constraints
    lb(2:nodes)=tofmin./tofscale;  % positive nonzero tof constraint (arbitrary 0.1-day min)
    
    JDscale=[JD0scale tofscale]; %convenient scale vector
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   

    nc=0;       %initialize constraint counter
    if ~isempty(tofTOTUB) %Total time of flight upper bound
        nc=nc+1;
    Acon(nc,:)=zeros(1,lenx0);
    Acon(nc,2:nodes)=tofscale;
    bcon(nc,1)=tofTOTUB;
    end
    
    if ~isempty(tofTOTLB) %Total time of flight lower bound
        nc=nc+1;
    Acon(nc,:)=zeros(1,lenx0);
    Acon(nc,2:nodes)=-tofscale;
    bcon(nc,1)=-tofTOTLB;
    end
    
    if ~isempty(JDUB) % node Julian Date upper bound
        icon=find(JDUB);
        for ni=1:length(icon)
        nc=nc+1;
    Acon(nc,:)=zeros(1,lenx0);
    Acon(nc,1:icon(ni))=JDscale(1:icon(ni));
    bcon(nc,1)=JDUB(icon(ni))-JD0;
        end
    end
    
    if ~isempty(JDLB) % node Julian Date upper bound
        icon=find(JDLB);
        for ni=1:length(icon)
        nc=nc+1;
    Acon(nc,:)=zeros(1,lenx0);
    Acon(nc,1:icon(ni))=-JDscale(1:icon(ni));
    bcon(nc,1)=-JDLB(icon(ni))+JD0;
        end
    end
    
    if ephtype==1 %check for DE405 or 2nd order ephemeris
        nc=nc+1;
        Acon(nc,:)=zeros(1,lenx0);
        Acon(nc,1:nodes)=-JDscale;
        bcon(nc,1)=-2414992+JD0; %apply DE405 lower bound
        nc=nc+1;
        Acon(nc,:)=zeros(1,lenx0);
        Acon(nc,1:nodes)=JDscale;
        bcon(nc,1)=2525008-JD0; %apply DE405 upper bound
    else
        nc=nc+1;
        Acon(nc,:)=zeros(1,lenx0);
        Acon(nc,1:nodes)=-JDscale;
        bcon(nc,1)=-625673+JD0; %apply 3000 BC lower bound
        nc=nc+1;
        Acon(nc,:)=zeros(1,lenx0);
        Acon(nc,1:nodes)=JDscale;
        bcon(nc,1)=2816787-JD0; %apply 3000 AD upper bound
    end
           %Apply positive time of flight constraints with Ax<b
%     Acon=[Acon; zeros(trajs,1) -eye(trajs)]; 
%     bcon=[bcon;zeros(trajs,1)];
           %alternatively use lb and up

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   
    
    nc=0;       %initialize constraint counter
    if ~isempty(tofTOTEQ) %Total time of flight equality
        nc=nc+1;
    Aeq(nc,:)=zeros(1,lenx0);
    Aeq(nc,2:nodes)=tofscale;
    beq(nc,1)=tofTOTEQ;
    end
        
    if ~isempty(tofEQ) %Each individual time of flight equality
        for ni=1:length(tofEQ)
            if tofEQ(ni)
        nc=nc+1;
    Aeq(nc,:)=zeros(1,lenx0);
    Aeq(nc,ni+1)=tofscale(ni);
    beq(nc,1)=tofEQ(ni);
            end
        end
    end
    
    if ~isempty(JDEQ) %Final node Julian Date upper bound
        icon=find(JDEQ);
        for ni=1:length(icon)
        nc=nc+1;
    Aeq(nc,:)=zeros(1,lenx0);
    Aeq(nc,1:icon(ni))=JDscale(1:icon(ni));
    beq(nc,1)=JDEQ(icon(ni))-JD0;
        end
    end
    
    if ~isempty(MCMEQ) %Midcourse Maneuver (X,Y,Z) equality
        nm=0;
        [row,col]=size(MCMEQ);
        for ni=1:col
            if (event(ni)==40)||(event(ni)==41)
                nm=nm+1;
                if any(MCMEQ(:,ni))
        nc=nc+1;
        nc3=nc+2;
        n3=nodes+3*nm;
    Aeq(nc:nc3,:)=zeros(3,lenx0);  %Set MCM X equality
    Aeq(nc:nc3,n3-2:n3)=[MCMscale(1,ni) 0 0; 0 MCMscale(2,ni) 0; 0 0 MCMscale(3,ni)];
    beq(nc:nc3,1)=MCMEQ(:,ni);    
        nc=nc3;
                end
            end
        end
    end
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   
        
    if ~isempty(JD0LB) %Initial node Julian Date lower bound    
        lb(1)= (JD0LB-JD0)/JD0scale;
    end
    if ~isempty(JD0UB) %Initial node Julian Date upper bound    
        ub(1)= (JD0UB-JD0)/JD0scale;
    end
    
    if ~isempty(tofLB) %Each individual time of flight lower bound
        for ni=1:length(tofLB)
           if tofLB(ni)
        lb(ni+1)= tofLB(ni)./tofscale(ni);
           end
        end
    end
    
    if ~isempty(tofUB) %Each individual time of flight upper bound   
        for ni=1:length(tofUB)
           if tofUB(ni)
        ub(ni+1)= tofUB(ni)./tofscale(ni);
           end
        end
    end
   
    if ~isempty(MCMLB)    %Set MCM (X,Y,Z) lower bound
        nm=0;
        [row,col]=size(MCMLB);
        for ni=1:col
            if (event(ni)==40)||(event(ni)==41)
                nm=nm+1;
                if any(MCMLB(:,ni))
    n3=nodes+3*nm;
    lb(n3-2:n3)=MCMLB(:,ni)./MCMscale(:,ni);
                end
            end
        end
    end
    if ~isempty(MCMUB)    %Set MCM (X,Y,Z) upper bound
        nm=0;
        [row,col]=size(MCMUB);
        for ni=1:col
            if (event(ni)==40)||(event(ni)==41)
                nm=nm+1;
                if any(MCMUB(:,ni))
    n3=nodes+3*nm;
    ub(n3-2:n3)=MCMUB(:,ni)./MCMscale(:,ni);
                end
            end
        end
    end
    
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   

	if ~isempty(MCMRUB)||~isempty(MCMRLB)||~isempty(MCMREQ)||~isempty(DVUB)||...
	~isempty(DVLB)||~isempty(DVEQ)||~isempty(vinfinUB)||~isempty(vinfinLB)||...
	~isempty(vinfinEQ)||~isempty(vinfoutUB)||~isempty(vinfoutLB)||...
    ~isempty(vinfoutEQ)||~isempty(vinfmatch)
	
	nonlcon=@(x)nonlconfunc(x,event,nplanet,grade,rev,rfbc,coep,coenp,JD0,nodes,...
	trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,dvopt,optnode,ephtype,Rj,muj,ecl,...
    JD0scale,tofscale,MCMscale,MCMRUB,MCMRLB,MCMREQ,DVUB,DVLB,DVEQ,vinfinUB,...
    vinfinLB,vinfinEQ,vinfoutUB,vinfoutLB,vinfinEQ,vinfmatch,maxstepJD0,...
    maxsteptof,maxstepMCM,tofmin);
	end

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    
end

xint=x0; %set global xint for optimization checks/options

rev0=rev; %record original revs
if any(rev)
    if branch  
      rev=abs(rev);
      col=find(rev);  %find all nonzeros index
      nbranch=nnz(rev); %count all nonzeros
      revs=repmat(rev,2^nbranch,1); % intitialize rev combination matrix
      branches=ones(2,nbranch);% intitialize branch combination matrix
      pbranch0=ones(1,nbranch);% intitialize default branch vector
      revs(2,:)=-revs(2,:);
      if nbranch>1
         branches(3:2+nbranch,1:nbranch)=~eye(nbranch)-eye(nbranch); 
      end
      if nbranch>2
          for ib=2:nbranch-1;
              pbranch=pbranch0;
              pbranch(1:ib)=-1;
              branches=[branches;uniqueperms(pbranch)];
          end
      end
      if nbranch>1
          for ib=1:nbranch
              revs(:,col(ib))=revs(:,col(ib)).*branches(:,ib);
          end 
      end

ibmax=2^nbranch; %set counter for number of permutations of branches
fvalbest=inf; %initialize current best performance
optoptions0=optoptions; %save original options settings
for ib=1:ibmax
    if branchopt>0
    optoptions.MaxIter=branchopt;
    end
    rev(1,:)=revs(ib,:)
    if optimization==1 %optimize intial date and times of flight for min. delta-V 
                       %UNCONSTRAINED GRADIENT METHOD
    [x,fval,flag,output]=fminunc(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,optoptions);
    elseif optimization==2 %optimize intial date and times of flight for min. delta-V 
                       %SIMPLEX METHOD
    [x,fval,flag,output]=fminsearch(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,optoptions);
    elseif optimization==3 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GRADIENT METHOD (INTERIOR POINT/SQP/ect.)
    [x,fval,flag,output]=fmincon(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
    elseif optimization==4 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED SIMPLEX METHOD 
    [x,fval,flag,output]=fminsearchcon(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0',lb,ub,Acon,bcon,Aeq,beq,nonlcon,strictness,optoptions);
        x=x'; %switch x orientation from fminsearchcon
    elseif optimization==5 %optimize intial date and times of flight for min. (delta-V)^2 
                       %UNCONSTRAINED NONLINEAR SYSTEM OF EQUATIONS SOLVER
    [x,fval,flag,output]=fsolve(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,optoptions);    
     elseif optimization==6 %optimize intial date and times of flight for min. (delta-V)^2 
                       %CONSTRAINED NONLINEAR LEAST-SQUARES 
    [x,fval,flag,output]=lsqnonlin(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,lb,ub,optoptions);   
    elseif optimization==7 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED MINIMIZATION OF THE MAXIMUM DELTA-V MANEUVER
    [x,fval,flag,output]=fminimax(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
    elseif optimization==10 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED PATTERN SEARCH (DIRECTED GLOBAL OPTIMIZATION)
    [x,fval,flag,output]=patternsearch(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
    elseif optimization==11 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED SIMULATED ANNEALING (DIRECTED GLOBAL OPTIMIZATION)
    [x,fval,flag,output]=simulannealbnd(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,lb,ub,optoptions);
    elseif optimization==12 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GENETIC ALGORITHM (STOCHASTIC GLOBAL OPTIMIZATION)
        if branchopt>0
        optoptions.Generations=branchopt;
        end
    [x,fval,flag,output]=ga(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),length(x0),Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
    elseif optimization==0
     fval=fmin(x0,event,nplanet,grade,rev,rfbc,coep,coenp,JD0,nodes,trajs,...
     vinfavmag,ilv,altp,nmcm,nsb,optr,dvopt,optnode,ephtype,Rj,muj,...
     ecl,JD0scale,tofscale,MCMscale,maxstepJD0,maxsteptof,maxstepMCM,optimization);
     x=x0;
    end
    
    if fval<fvalbest
       xbest=x;
       revbest=rev;
       fvalbest=fval;
    end   
end
revbest %display the resulting best found rev combination

optoptions=optoptions0; %reset original optimization options;
    if branchopt==0
        x=xbest;
        rev=revbest;
        optimization=99;
        rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
        JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
        JD(1)=JD0+x(1)*JD0scale;
        tof=x(2:(nodes)).*tofscale;
                tof=abs(tof); 
        nx=1;            %initialize counting vector   
        if nmcm||nsb
           for n=1:nodes
              if (event(n)==40)||(event(n)==41)
                  rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
                  nx=nx+3;
              end
           end
        end
    else
    x0=xbest;
    rev=revbest;
    end      

else
  nbranch=0; 
end

end

%--------------------------------------------------------------------------
tic
if optimization==1 %optimize intial date and times of flight for min. delta-V 
                   %UNCONSTRAINED GRADIENT METHOD
    [x,fval,flag,output]=fminunc(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale;
            tof=abs(tof); 
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end    
    
elseif optimization==2 %optimize intial date and times of flight for min. delta-V 
                       %SIMPLEX METHOD
    [x,fval,flag,output]=fminsearch(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,optoptions);
    rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale;
            tof=abs(tof);  
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end      
    
elseif optimization==3 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GRADIENT METHOD (ACTIVE SET/SQP/ect.)
                       
    [x,fval,flag,output]=fmincon(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end    

elseif optimization==4 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED SIMPLEX METHOD 
                       
    [x,fval,flag,output]=fminsearchcon(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0',lb,ub,Acon,bcon,Aeq,beq,nonlcon,strictness,optoptions);
    
    x=x'; %switch x orientation from fminsearchcon
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end  
    
elseif optimization==5 %optimize intial date and times of flight for min. (delta-V)^2 
                       %CONSTRAINED NONLINEAR LEAST-SQUARES 
                       
    [x,fval,flag,output]=fsolve(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end      
    
elseif optimization==6 %optimize intial date and times of flight for min. (delta-V)^2 
                       %CONSTRAINED NONLINEAR LEAST-SQUARES 
                       
    [x,fval,flag,output]=lsqnonlin(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,lb,ub,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end  
    
elseif optimization==7 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GRADIENT METHOD (ACTIVE SET/SQP/ect.)
                       
    [x,fval,flag,output]=fminimax(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end        
    
elseif optimization==10 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GRADIENT METHOD (ACTIVE SET/SQP/ect.)
                       
    [x,fval,flag,output]=patternsearch(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end      
    
elseif optimization==11 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GRADIENT METHOD (ACTIVE SET/SQP/ect.)
                       
    [x,fval,flag,output]=simulannealbnd(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),x0,lb,ub,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end  
    
elseif optimization==12 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GRADIENT METHOD (ACTIVE SET/SQP/ect.)
                       
   [x,fval,flag,output]=ga(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),length(x0),Acon,bcon,Aeq,beq,lb,ub,nonlcon,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end     
    
elseif optimization==13 %optimize intial date and times of flight for min. delta-V 
                       %CONSTRAINED GRADIENT METHOD (ACTIVE SET/SQP/ect.)
                       
    [x,fval,flag,output]=particleswarm(@(x)fmin(x,event,nplanet,grade,rev,...
     rfbc,coep,coenp,JD0,nodes,trajs,vinfavmag,ilv,altp,nmcm,nsb,optr,...
     dvopt,optnode,ephtype,Rj,muj,ecl,JD0scale,tofscale,MCMscale,maxstepJD0,...
     maxsteptof,maxstepMCM,optimization),length(x0),lb,ub,optoptions);
     rnb(1:3,1:nodes)=0; %initialize/preallocate matrix for efficiency 
    JD=zeros(nodes,1);  %initialize/preallocate vector for efficiency     
    JD(1)=JD0+x(1)*JD0scale;
    tof=x(2:(nodes)).*tofscale; 
        
    nx=1;            %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end      
    
elseif optimization==0 %single trajectory determination (no optimization)
    JD=zeros(nodes,1); 
    JD(1)=JD0; 
    output={'None';1;'None';'None'};
    flag = 'None';
    
elseif optimization==99 %Selected fully optimized branch permutation
    %Nothing Necessary Here YET
else
    disp('Selected optimization method is NOT an available choice.');
    disp('Please choose between methods: 0,1,2, and 3');
end
extime=toc

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%-------------------Post Optimization Calculations-------------------------
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
rP=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vd=zeros(3,nodes);          %initialize/preallocate vector for efficiency
va=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vinfout=zeros(3,nodes);     %initialize/preallocate vector for efficiency
vinfin=zeros(3,nodes);      %initialize/preallocate vector for efficiency
vinfouthat=zeros(3,nodes);  %initialize/preallocate vector for efficiency
vinfinhat=zeros(3,nodes);   %initialize/preallocate vector for efficiency
vinfoutmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
vinfinmag=zeros(nodes,1);     %initialize/preallocate vector for efficiency
rsoi=zeros(nodes,1);          %initialize/preallocate vector for efficiency
rpmag=zeros(nodes,1);         %initialize/preallocate vector for efficiency
dv=zeros(3,nodes);          %initialize/preallocate vector for efficiency
dve=zeros(3,nodes);         %initialize/preallocate vector for efficiency
rm=zeros(3,nodes);          %initialize/preallocate vector for efficiency
rme=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vin=zeros(3,nodes);         %initialize/preallocate vector for efficiency
vout=zeros(3,nodes);        %initialize/preallocate vector for efficiency
vine=zeros(3,nodes);         %initialize/preallocate vector for efficiency
voute=zeros(3,nodes);        %initialize/preallocate vector for efficiency
vinfpl=zeros(3,nodes);      %initialize/preallocate vector for efficiency
vPplhat=zeros(3,nodes);     %initialize/preallocate vector for efficiency    
vplhat=zeros(3,nodes);      %initialize/preallocate vector for efficiency
betaout=zeros(nodes,1);       %initialize/preallocate vector for efficiency
betain=zeros(nodes,1);        %initialize/preallocate vector for efficiency
deltaout=zeros(nodes,1);      %initialize/preallocate vector for efficiency
deltain=zeros(nodes,1);       %initialize/preallocate vector for efficiency
Tp=zeros(nodes,1);            %initialize/preallocate vector for efficiency
tsoiout=zeros(nodes,1);       %initialize/preallocate vector for efficiency
tsoiin=zeros(nodes,1);       %initialize/preallocate vector for efficiency
rpfb=zeros(3,nodes);        %initialize/preallocate vector for efficiency
rpin=zeros(3,nodes);        %initialize/preallocate vector for efficiency
rpout=zeros(3,nodes);       %initialize/preallocate vector for efficiency
vpin=zeros(3,nodes);        %initialize/preallocate vector for efficiency
vpout=zeros(3,nodes);       %initialize/preallocate vector for efficiency
dvmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
rpfbmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
rpfbmagr=zeros(nodes,1);    %initialize/preallocate vector for efficiency
rmmagr=zeros(nodes,1);    %initialize/preallocate vector for efficiency
deltainr=zeros(nodes,1);    %initialize/preallocate vector for efficiency
deltaoutr=zeros(nodes,1);    %initialize/preallocate vector for efficiency
vinmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
voutmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
vpinmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
vpoutmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
delin=zeros(nodes,1);    %initialize/preallocate vector for efficiency
delout=zeros(nodes,1);    %initialize/preallocate vector for efficiency
dvRA=zeros(nodes,1);    %initialize/preallocate vector for efficiency
dvDEC=zeros(nodes,1);    %initialize/preallocate vector for efficiency
npp=zeros(3,nodes);    %initialize/preallocate vector for efficiency
hpo=zeros(3,nodes);    %initialize/preallocate vector for efficiency
hhat=zeros(3,nodes);    %initialize/preallocate vector for efficiency
coet=zeros(6,nodes);    %initialize/preallocate vector for efficiency
exitflag=[];                %initialize lambert error check

%Compute first node state
[rP(:,1),vP(:,1)]=ephstate(ephtype,event(1),nplanet(1),JD(1),coenp(:,1),rnb(:,1));

for n=1:trajs
    np1=n+1;
    if tof(n)<=(tofmin*0.99) %check if tof = 0 and add small perturbation 
        disp(['tof ',num2str(n),' is near-zero and is set to 0.1 days'])
        tof(n)=tofmin; %arbitrarily small perturbation (days)
    end
 JD(np1)=JD(n)+tof(n); %Determine julian dates from initial JD and times of flight  

%Collect ephemeris and state data for other nodes
[rP(:,np1),vP(:,np1)]=ephstate(ephtype,event(np1),nplanet(np1),JD(np1),coenp(:,np1),rnb(:,np1));

%Hybrid lambert targeter with Dr. Izzo (ESA) fast convergence solution and 
%Lancaster & Blancard's robust solution 
 [vd(:,n),va(:,np1),exitflag]=lambert(rP(:,n),rP(:,np1),tof(n)*86400,rev(n),grade(n),muj);
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
        tof(n)=tof(n)+0.1; %small perturbation (days)
        JD(n+1)=JD(n)+tof(n); %adjust JD to new tof
    end
    
    %Collect ephemeris and state data for other nodes
[rP(:,np1),vP(:,np1)]=ephstate(ephtype,event(np1),nplanet(np1),JD(np1),coenp(:,np1),rnb(:,np1));

%Hybrid lambert targeter with Dr. Izzo (ESA) fast convergence solution and 
%Lancaster & Blancard's robust solution 
 [vd(:,n),va(:,np1),exitflag]=lambert(rP(:,n),rP(:,np1),tof(n)*86400,rev(n),grade(n),muj);
 if (exitflag ~=1)||any(isnan(vd(:,n)))||any(isnan(va(:,np1)))
     if exitflag==1 %Catch any NaN cases that exitflag didn't
     exitflag=-1;
     end
     disp(['!! Lambert Targeter failed between planet ',num2str(nplanet(n)),' and planet ',num2str(nplanet(np1)),'!!'])
 end

nex=nex+1;

    if exitflag~=1 && nex>=10
        disp(['!! LAMBERT ARC FOR LEG ',num2str(n),' HAS NO SOLUTION FOR ',num2str(rev(n)),' REVS AND ',num2str(tof(n)),' days TOF !!'])
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
vinfouthat(1:3,n)=vP(:,n)/sqrt(vP(:,n)'*vP(:,n));
else
vinfouthat(:,n) = vinfout(:,n)/vinfoutmag(n);
end
if (vinfinmag(np1)==0) %check for vinf=0 case
vinfinhat(1:3,np1)=vP(:,np1)/sqrt(vP(:,np1)'*vP(:,np1));
else
vinfinhat(:,np1) = vinfin(:,np1)/vinfinmag(np1);
end
end

[mupl,rpl,rap,dp,wp]=planetinfo(nplanet,JD,nodes); %Acquire planetary data 
%Determine sphere of influence and planetocentric periapsis radii
for n=1:nodes
rsoi(n) = (mupl(n)/muj)^(2/5)*norm(rP(:,n));      
rpmag(n)=rfbc(n)*rpl(n);
end
for n=1:nodes;
%--------------------------Departure Trajectory----------------------------
    if event(n) == 10 %parking orbit with variable i

    [dv(:,n),dve(:,n),rm(:,n),vin(:,n),rme(:,n),vine(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vout(:,n),voute(:,n),betaout(n),deltaout(n),Tp(n),tsoiout(n),coep(1:6,n)]...
    = park2hypopt(vP(:,n),vd(:,n),vinfouthat(:,n),vinfoutmag(n),rap(n),dp(n),coep(1:3,n),...
    mupl(n),rsoi(n));
    
    elseif event(n) == 11 %parking orbit with fixed i
 
     [dv(:,n),dve(:,n),rm(:,n),vin(:,n),rme(:,n),vine(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vout(:,n),voute(:,n),betaout(n),deltaout(n),Tp(n),tsoiout(n),coep(1:6,n)]...
    = park2hypd(vP(:,n),vd(:,n),vinfouthat(:,n),vinfoutmag(n),rap(n),dp(n),coep(1:3,n),...
    mupl(n),rsoi(n));

    elseif event(n) == 12 %parking orbit with variable true anomaly

        
        
    elseif event(n) == 13 %fixed parking orbit
        
    [dv(:,n),dve(:,n),rm(:,n),vin(:,n),rme(:,n),vine(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vout(:,n),voute(:,n),betaout(n),deltaout(n),Tp(n),tsoiout(n)]...
    = park2hypfix(vP(:,n),vd(:,n),vinfouthat(:,n),vinfoutmag(n),rap(n),dp(n),coep(1:6,n),...
    mupl(n),rsoi(n));

    elseif event(n) == 14 %launch vehicle trajectory

     [dv(:,n),dve(:,n),rm(:,n),vin(:,n),rme(:,n),vine(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vout(:,n),voute(:,n),betaout(n),deltaout(n),Tp(n),tsoiout(n),coep(1:6,n)]...
    = vinf2hypd(vP(:,n),vd(:,n),vinfouthat(:,n),vinfoutmag(n),rpl(n),rap(n),dp(n),...
    vinfavmag(n),ilv(n),altp(n),mupl(n),rsoi(n),dvopt(n));  

    elseif (event(n) == 17) %Depart User-Defined Body Trajectory
    
    rm(:,n)=rP(:,n);    
    rme(:,n)=rm(:,n);
    vin(:,n)=vP(:,n);
    vine(:,n)=vin(:,n);
    vout(:,n)=vd(:,n);
    voute(:,n)=vout(:,n);
    dv(:,n)=vout(:,n)-vin(:,n);
    dve(:,n)=dv(:,n);
      
    elseif (event(n) == 18) %Depart from Parking orbit about User-Defined Body 
        
%---------------------------Arrival Trajectory-----------------------------        
    elseif event(n) == 20 %parking orbit with variable i
   
    [dv(:,n),dve(:,n),rm(:,n),vout(:,n),rme(:,n),voute(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vin(:,n),vine(:,n),betain(n),deltain(n),Tp(n),tsoiin(n),coep(1:6,n)]...
    = hyp2parkopt(vP(:,n),va(:,n),vinfinhat(:,n),vinfinmag(n),rap(n),dp(n),...
    coep(1:3,n),mupl(n),rsoi(n));   
 
    elseif event(n) == 21 %parking orbit with fixed i
  
    [dv(:,n),dve(:,n),rm(:,n),vout(:,n),rme(:,n),voute(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vin(:,n),vine(:,n),betain(n),deltain(n),Tp(n),tsoiin(n),coep(1:6,n)]...
    = hypa2park(vP(:,n),va(:,n),vinfinhat(:,n),vinfinmag(n),rap(n),dp(n),coep(1:3,n),...
    mupl(n),rsoi(n)); 
   
    elseif event(n) == 22 %parking orbit with variable true anomaly
        
    
    elseif event(n) == 23 %fixed parking orbit
    
    [dv(:,n),dve(:,n),rm(:,n),vout(:,n),rme(:,n),voute(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vin(:,n),vine(:,n),betain(n),deltain(n),Tp(n),tsoiin(n)]...
    = hyp2parkfix(vP(:,n),va(:,n),vinfinhat(:,n),vinfinmag(n),rap(n),dp(n),coep(1:6,n),...
    mupl(n),rsoi(n));     

     elseif event(n) == 24 %entry trajectory
        
     [dv(:,n),dve(:,n),rm(:,n),vout(:,n),rme(:,n),voute(:,n),vinfpl(:,n),vPplhat(:,n),...
     vplhat(:,n),vin(:,n),vine(:,n),betain(n),deltain(n),Tp(n),tsoiin(n),coep(1:6,n)]...
    = hypa2vinf(vP(:,n),va(:,n),vinfinhat(:,n),vinfinmag(n),rpl(n),rap(n),dp(n),...
    vinfavmag(n),ilv(n),altp(n),mupl(n),rsoi(n),dvopt(n)); 
 
    elseif event(n) == 27 %Redezvous with User-Defined Body Trajectory
        
    rm(:,n)=rP(:,n);    
    rme(:,n)=rm(:,n);
    vin(:,n)=va(:,n);
    vine(:,n)=vin(:,n);
    vout(:,n)=vP(:,n);
    voute(:,n)=vout(:,n);
    dv(:,n)=vout(:,n)-vin(:,n);    
    dve(:,n)=dv(:,n);  
    
    elseif event(n) == 28 %Arrival at Parking orbit about User-Defined Body 
        
        
%----------------------------Flyby Trajectory------------------------------     
    elseif event(n) == 30 %unpowered flyby (account for residual DV @ periapse)
        
    [dv(:,n),rm(:,n),rpfb(:,n),vin(:,n),vout(:,n),rpin(:,n),rpout(:,n),vpin(:,n),...
     vpout(:,n),betain(n),betaout(n),deltain(n),deltaout(n),tsoiin(n),tsoiout(n)]...
    =flybyunp(rpmag(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),vinfouthat(:,n),vinfoutmag(n),...
    mupl(n),rsoi(n),optr(n));

    rme(:,n)=rm(:,n);    
    vine(:,n)=vin(:,n);
    voute(:,n)=vout(:,n);
    dve(:,n)=dv(:,n);
    
    elseif event(n) == 31 %unpowered flyby (account for residual DV @ periapse)
        
    [dv(:,n),rm(:,n),rpfb(:,n),vin(:,n),vout(:,n),rpin(:,n),rpout(:,n),vpin(:,n),...
     vpout(:,n),betain(n),betaout(n),deltain(n),deltaout(n),tsoiin(n),tsoiout(n)]...
    =flybysoi(rpmag(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),vinfouthat(:,n),vinfoutmag(n),...
    mupl(n),rsoi(n),optr(n));

    rme(:,n)=rm(:,n);    
    vine(:,n)=vin(:,n);
    voute(:,n)=vout(:,n);
    dve(:,n)=dv(:,n);
        
    elseif event(n) == 33 %optimal powered flyby
   
    [dv(:,n),rm(:,n),rpfb(:,n),vin(:,n),vout(:,n),rpin(:,n),rpout(:,n),vpin(:,n),...
     vpout(:,n),betain(n),betaout(n),deltain(n),deltaout(n),tsoiin(n),tsoiout(n)]...
    =flybyp(rpmag(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),vinfouthat(:,n),vinfoutmag(n),...
    mupl(n),rsoi(n),optr(n));
   
    rme(:,n)=rm(:,n);    
    vine(:,n)=vin(:,n);
    voute(:,n)=vout(:,n);
    dve(:,n)=dv(:,n);
   
    elseif event(n) == 35 %periapse powered flyby (residual delta-V at periapsis)
   
    [dv(:,n),rm(:,n),rpfb(:,n),vin(:,n),vout(:,n),rpin(:,n),rpout(:,n),vpin(:,n),...
     vpout(:,n),betain(n),betaout(n),deltain(n),deltaout(n),tsoiin(n),tsoiout(n)]...
    =flybypp(rpmag(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),vinfouthat(:,n),vinfoutmag(n),...
    mupl(n),rsoi(n),optr(n));
    
    rme(:,n)=rm(:,n);    
    vine(:,n)=vin(:,n);
    voute(:,n)=vout(:,n);
    dve(:,n)=dv(:,n);
    
%----------------------------Spaceburn Trajectory--------------------------
    elseif (event(n) == 40)||(event(n) == 41)    

    rm(:,n)=rP(:,n);
    rme(:,n)=rm(:,n);
    vin(:,n)=va(:,n);
    vine(:,n)=vin(:,n);
    vout(:,n)=vd(:,n);
    voute(:,n)=vout(:,n);
    dv(:,n)=vout(:,n)-vin(:,n); 
    dve(:,n)=dv(:,n);
    end
dvmag(n)=norm(dv(:,n));
rpfbmag(n)=norm(rpfb(:,n));
rpfbmagr(n)=rpfbmag(n)/rpl(n);
rmmagr(n)=norm(rm(:,n))/rpl(n);
deltainr(n)=deltain(n)/rpl(n);
deltaoutr(n)=deltaout(n)/rpl(n);
vinmag(n)=norm(vin(:,n));
voutmag(n)=norm(vout(:,n));
vpinmag(n)=norm(vpin(:,n));
vpoutmag(n)=norm(vpout(:,n));
delin(n)=2*asin(cos(betain(n)));
delout(n)=2*asin(cos(betaout(n)));
dvRA(n)=atan2(dv(2,n),dv(1,n))*180/pi;
dvDEC(n)=asin(dv(3,n)/dvmag(n))*180/pi;

        %Planets' north pole vector
npp(:,n)=[cosd(rap(n))*cosd(dp(n));
cosd(ecl)*sind(rap(n))*cosd(dp(n))+sind(ecl)*sind(dp(n)); 
-sind(ecl)*sind(rap(n))*cosd(dp(n))+cosd(ecl)*sind(dp(n))];

        %Parking Orbit angular momentum vector direction h
hpo(:,n)=cross(rme(:,n),voute(:,n)); 
hhat(:,n)=hpo(:,n)/norm(hpo(:,n));

end

for n=1:trajs
coet(1:6,n)=rv2orbel(rP(:,n),vd(:,n),muj);
end



% intmbody=0; %Juno mission integration test
% 
% if intmbody==1
% %    for n=1:trajs
%    
%    Zint(:,n)=[rme(:,n)+vP(:,n) voute(:,n)];
% [tsc(:,n),yint(n,1:500,1:6)]=ode45(@(t,y)EOM3body(t,y,mup,mus,JD(n),nplanet(n)),linspace(0,tof(n)*86400,500),...
%     Zint(:,n),options);
%    
% %    end
% end




%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%---------------------------Write Output File------------------------------
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if writeout==1  
writeoutput(dv,dvmag,dve,JD,tof,rsoi,vinfout,vinfin,vinfpl,rP,vP,vPplhat,vd,va,...
    vplhat,coet,dvRA,dvDEC,rm,rmmagr,vout,vin,voutmag,vinmag,rpout,rpin,vpout,vpin,vpoutmag,...
    vpinmag,rpfb,rpmag,rfbc,rpfbmagr,betaout,betain,delout,delin,deltaoutr,deltainr,tsoiout,...
    tsoiin,Tp,mupl,rpl,rap,dp,wp,muj,coep,coenp,rnb,rnb0,event,optimization,nplanet,...
    grade,rev,revig,vinfavmag,ilv,altp,optr,optnode,dvopt,optoptions,...
    ephtype,excelfile,nodes,trajs,extime,vinfoutmag,vinfinmag,output,flag)
end

if ~isempty(savefile)
    JD0=JD(1);
    rnb0=rnb;
    save(savefile,'x','JD0','tof','rnb0','event','nplanet','grade','rev','rfbc','coep',...
        'coenp','vinfavmag','ilv','altp','optr','dvopt','optnode','ephtype')
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%------------------------Optimization Function-----------------------------
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function f=fmin(x,event,nplanet,grade,rev,rfbc,coep,coenp,JD0,nodes,trajs,...
        vinfavmag,ilv,altp,nmcm,nsb,optr,dvopt,optnode,ephtype,Rj,muj,...
        ecl,JD0scale,tofscale,MCMscale,maxstepJD0,maxsteptof,maxstepMCM,optimization)
    x=x(:)'; %check for x orientation from fminsearchcon
JD=zeros(nodes,1);            %initialize/preallocate vector for efficiency
rP=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vP=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vd=zeros(3,nodes);          %initialize/preallocate vector for efficiency
va=zeros(3,nodes);          %initialize/preallocate vector for efficiency
vinfout=zeros(3,nodes);     %initialize/preallocate vector for efficiency
vinfin=zeros(3,nodes);      %initialize/preallocate vector for efficiency
vinfouthat=zeros(3,nodes);  %initialize/preallocate vector for efficiency
vinfinhat=zeros(3,nodes);   %initialize/preallocate vector for efficiency
vinfoutmag=zeros(nodes,1);    %initialize/preallocate vector for efficiency
vinfinmag=zeros(nodes,1);     %initialize/preallocate vector for efficiency
rsoi=zeros(nodes,1);          %initialize/preallocate vector for efficiency
rpmag=zeros(nodes,1);         %initialize/preallocate vector for efficiency
dv=zeros(nodes,1);          %initialize/preallocate vector for efficiency
rnb=zeros(3,nodes);         %initialize/preallocate matrix for efficiency 
exitflag=[];                %initialize lambert error check
f=zeros(length(x),1);                        %initialize performance index
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


    nx=1; %initialize counting vector   
    if nmcm||nsb
       for n=1:nodes
          if (event(n)==40)||(event(n)==41)
              rnb(:,n)=x(nodes+nx:nodes+nx+2)'.*MCMscale(:,n);
              nx=nx+3;
          end
       end
    end       
JD(1)=JD0+x(1)*JD0scale;
% jed2date(JD(1))

%Compute first node state
[rP(:,1),vP(:,1)]=ephstate(ephtype,event(1),nplanet(1),JD(1),coenp(:,1),rnb(:,1));

for n=1:trajs
    np1=n+1;
    if abs(x(np1))*tofscale(n)<(tofmin*0.99) %check if tof = 0 and add small perturbation 
        disp(['tof ',num2str(n),' is near-zero and is set to 0.1 days'])
        x(np1)=tofmin/tofscale(n); %arbitrarily small perturbation (days)
    end
 JD(np1)=JD(n)+abs(x(np1))*tofscale(n); %Determine julian dates from initial JD and times of flight  
%  jed2date(JD(np1))
%  x(1:nodes)
    %Collect ephemeris and state data for other nodes
[rP(:,np1),vP(:,np1)]=ephstate(ephtype,event(np1),nplanet(np1),JD(np1),coenp(:,np1),rnb(:,np1));

%Hybrid lambert targeter with Dr. Izzo (ESA) fast convergence solution and 
%Lancaster & Blancard's robust solution 
 [vd(:,n),va(:,np1),exitflag]=lambert(rP(:,n),rP(:,np1),(JD(np1)-JD(n))*86400,rev(n),grade(n),muj);
 if (exitflag ~=1)||any(isnan(vd(:,n)))||any(isnan(va(:,np1)))
     if exitflag==1 %Catch any NaN cases that exitflag didn't
     exitflag=-1;
     end
     disp(['!! Lambert Targeter failed between planet ',num2str(nplanet(n)),' and planet ',num2str(nplanet(np1)),'!!'])
 end

nex=1; %initialize exitflag loop counter
while exitflag~=1
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
 [vd(:,n),va(:,np1),exitflag]=lambert(rP(:,n),rP(:,np1),(JD(np1)-JD(n))*86400,rev(n),grade(n),muj);
 if (exitflag ~=1)||any(isnan(vd(:,n)))||any(isnan(va(:,np1)))
     if exitflag==1 %Catch any NaN cases that exitflag didn't
     exitflag=-1;
     end
     disp(['!! Lambert Targeter failed between planet ',num2str(nplanet(n)),' and planet ',num2str(nplanet(np1)),'!!'])
 end
nex=nex+1;

    if exitflag~=1 && nex>=10
        disp(['!! LAMBERT ARC FOR LEG ',num2str(n),' HAS NO SOLUTION FOR ',num2str(rev(n)),' REVS AND ',num2str(x(np1)*tofscale(n)),' days TOF !!'])
        break
    end
end

end
%  x,xint
xint=x; %reset x0 vector for next optimizer iteration

%Determine hyperbolic excess velocity (v-infinity) at each node
for n=1:trajs
    np1=n+1;
vinfout(:,n) = vd(:,n) - vP(:,n);
vinfin(:,np1) = va(:,np1) - vP(:,np1);
vinfoutmag(n) = sqrt(vinfout(:,n)'*vinfout(:,n));
vinfinmag(np1) = sqrt(vinfin(:,np1)'*vinfin(:,np1));
if (vinfoutmag(n)==0) %check for vinf=0 case
vinfouthat(1:3,n)=vP(:,n)/sqrt(vP(:,n)'*vP(:,n));    
else
vinfouthat(:,n) = vinfout(:,n)/vinfoutmag(n);
end
if (vinfinmag(np1)==0) %check for vinf=0 case
vinfinhat(1:3,np1)=vP(:,np1)/sqrt(vP(:,np1)'*vP(:,np1));
else
vinfinhat(:,np1) = vinfin(:,np1)/vinfinmag(np1);
end
end

[mupl,rpl,rap,dp,wp]=planetinfo(nplanet,JD,nodes); %Acquire planetary data 
%Determine sphere of influence and planetocentric periapsis radii
for n=1:nodes
rsoi(n) = (mupl(n)/muj)^(2/5)*sqrt(rP(:,n)'*rP(:,n));      
rpmag(n)=rfbc(n)*rpl(n);
end

for n=1:nodes;

%--------------------------Departure Trajectory----------------------------
    if event(n) == 10 %parking orbit with variable i

%     [dv(n)]= park2hypoptOPT(vinfoutmag(n),coep(1:2,n),mupl(n));
    vcmag=(mupl(n)/(coep(1,n)*(1-coep(2,n)))); %velocity^2 of circular orbit with same periapsis
    dv(n)=abs(sqrt(2*vcmag+vinfoutmag(n)*vinfoutmag(n))-sqrt(vcmag*(1+coep(2,n))));
       
    f(n)=dv(n)*optnode(n);
       
        
    elseif event(n) == 11 %parking orbit with fixed i
 
     [dv(n)]= park2hypdOPT(vinfouthat(:,n),vinfoutmag(n),rap(n),dp(n),coep(1:3,n),mupl(n));    
        
    f(n)=dv(n)*optnode(n);
        
    elseif event(n) == 12 %parking orbit with variable true anomaly
     
            f(n)=0;

        
    elseif event(n) == 13 %fixed parking orbit
     
     [dv(n)]= park2hypfixOPT(vinfouthat(:,n),vinfoutmag(n),rap(n),dp(n),coep(1:6,n),mupl(n));    
        
    f(n)=dv(n)*optnode(n);
        
    elseif event(n) == 14 %launch vehicle trajectory
 
     [dv(n)]= vinf2hypdOPT(vinfouthat(:,n),vinfoutmag(n),rpl(n),rap(n),dp(n),...
         vinfavmag(n),ilv(n),altp(n),mupl(n),dvopt(n));
        
    f(n)=dv(n)*optnode(n);
        
            
    elseif (event(n) == 17) %Redezvous with User-Defined Body Trajectory
    
    dv(n)=sqrt((vd(:,n)-vP(:,n))'*(vd(:,n)-vP(:,n)));
        
    f(n)=dv(n)*optnode(n);
        
    elseif (event(n) == 18) % Orbit about User-Defined Body
        
        
    f(n)=0;
         
        
%---------------------------Arrival Trajectory-----------------------------        
    elseif event(n) == 20 %parking orbit with variable i
    
%     [dv(n)]= hyp2parkoptOPT(vinfinmag(n),coep(1:2,n),mupl(n)); 
    vcmag=(mupl(n)/(coep(1,n)*(1-coep(2,n)))); %velocity^2 of circular orbit with same periapsis 
    dv(n)=abs(sqrt(2*vcmag+vinfinmag(n)*vinfinmag(n))-sqrt(vcmag*(1+coep(2,n))));
    
    f(n)=dv(n)*optnode(n);
        
        
    elseif event(n) == 21 %parking orbit with fixed i
   
    [dv(n)]= hypa2parkOPT(vinfinhat(:,n),vinfinmag(n),rap(n),dp(n),coep(1:3,n),mupl(n)); 
        
    f(n)=dv(n)*optnode(n);
    
    
    elseif event(n) == 22 %parking orbit with variable true anomaly
     
        
        
    elseif event(n) == 23 %fixed parking orbit
     
     [dv(n)]=  hyp2parkfixOPT(vinfinhat(:,n),vinfinmag(n),rap(n),dp(n),coep(1:6,n),mupl(n)); 
        
    f(n)=dv(n)*optnode(n);
        
    elseif event(n) == 24 %entry trajectory
        
     [dv(n)]= hypa2vinfOPT(vinfinhat(:,n),vinfinmag(n),rpl(n),rap(n),dp(n),...
         vinfavmag(n),ilv(n),altp(n),mupl(n),dvopt(n));        
        
    f(n)=dv(n)*optnode(n);
        
 
    elseif (event(n) == 27) %Redezvous with User-Defined Body Trajectory
    
    dv(n)=sqrt((vP(:,n)-va(:,n))'*(vP(:,n)-va(:,n)));      

    f(n)=dv(n)*optnode(n);
        
        
    elseif (event(n) == 28) % Orbit about User-Defined Body
        
        
    f(n)=dv(n)*optnode(n);
          
%----------------------------Flyby Trajectory------------------------------        
    elseif event(n) == 30 %unpowered gravity assist
        
%         dv(n)=0; %assume no deltaV with constraints satisfied
    f(n)=0;
    
        
    elseif event(n) == 31 %unpowered gravity assist w/ dv @ Rinf  
        
    [dv(n)]=flybyinfOPT(rpmag(n),rsoi(n),vinfinhat(:,n),vinfinmag(n),vinfouthat(:,n),...
    vinfoutmag(n),mupl(n),optr(n));   

    f(n)=dv(n)*optnode(n);

    elseif event(n) == 33 %powered optimal gravity assist 
    [dv(n)]=flybypOPT(rpmag(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),...
    vinfouthat(:,n),vinfoutmag(n),mupl(n),rsoi(n),optr(n));
    
    f(n)=dv(n)*optnode(n);
        
        
    elseif event(n) == 35 %powered periapse gravity assist (residual delta-V at periapsis)

    [dv(n)]=flybyppOPT(rpmag(n),rsoi(n),rP(:,n),vinfinhat(:,n),vinfinmag(n),...
    vinfouthat(:,n),vinfoutmag(n),mupl(n),optr(n));  
    
    f(n)=dv(n)*optnode(n);
        
%----------------------------Midcourse Maneuver Trajectory--------------------------    
    elseif (event(n) == 40)||(event(n) == 41)   
            
    dv(n)=sqrt((vd(:,n)-va(:,n))'*(vd(:,n)-va(:,n))); 

    f(n)=dv(n)*optnode(n);
       

    end
end

if (optimization==5)||(optimization==6)||(optimization==7)
    %keep f in vector form
else
    f=sum(f);
end

end








%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%---------------------------Plotting Functions-----------------------------
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>





if ~isempty(movie)
writerObj = VideoWriter(strcat(movie,'.avi'),'Motion JPEG AVI');
writerObj.FrameRate=fps;
open(writerObj);
end
pn=0; %initialize current plot number
nplanetot=length(plotting); % total number of plot commands
options=odeset('RelTol',1e-8,'AbsTol',1e-8);
for i=1:nplanetot
%-------------------Plot Heliocentric Trajectories-------------------------    
    if (plotting(i) >= 0)&&(plotting(i) < 10)
        %Spacecraft trajectories
       for n=1:trajs
Zsc(n,:)=[rP(:,n)' vd(:,n)'];
[tsc(n,:),ysc(n,1:nptsSC,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,muj),linspace(0,tof(n)*86400,nptsSC),...
    Zsc(n,:),options);
       end
      
       %Planetary Trajectories
       for n=1:nodes
       if nplanet(n)==0
            if (event(n)==40)||(event(n)==41)
 rP0(:,n)=rP(:,n);
 vP0(1:3,n)=0;
 TP(n)=0;
 tP(n,1)=0;
 yP(n,1,1:6)=[rP0(:,n)' vP0(:,n)'];
            elseif (event(n)==17)||(event(n)==27)
                if coenp(2,n)>=1.0
                 rP0(:,n)=rP(:,n);
                 vP0(:,n)=vP(:,n);
                 TP(n)=Inf;
                 tP(n,1)=0;
                 yP(n,1,1:6)=[rP0(:,n)' vP0(:,n)'];
                else
% coenp(6,n)=kepdt(coenp(1:7,n)',JD(1),muj);   
[rP0(:,n), vP0(:,n)]=orbel2rv(muj,coenp(1,n),coenp(2,n),coenp(3,n),coenp(4,n),...
    coenp(5,n),coenp(6,n));     

[rP0(:,n), vP0(:,n)]=kepuv(muj,(JD(1)-coenp(7,n))*86400,rP0(:,n),vP0(:,n),0,dxtol);
TP(n)=2*pi*sqrt(coenp(1,n)^3/muj); 
ZP=[rP(:,n)' vP(:,n)'];
[tP(n,:),yP(n,1:nptsUD,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,muj),linspace(0,TP(n),nptsUD),ZP,options);       

                end
            end
       else    
           %Determine planets' orbital elements and period
[aP(n),eP(n),iP(n),oP(n),wP(n),taP(n),CaseP(n)] = rv2orbel(rP(:,n),vP(:,n),muj);
TP(n)=2*pi*sqrt(aP(n)^3/muj);          
            %Determine location of planets at initial julian date
[rP0(:,n),vP0(:,n)]=ephstate(ephtype,event(n),nplanet(n),JD(1),coenp(:,n),rnb(:,n));
        
ZP=[rP(:,n)' vP(:,n)'];
[tP(n,:),yP(n,1:nptsP,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,muj),linspace(0,TP(n),nptsP),ZP,options);      
    
       end
       end
       
    if (plotting(i)>=0)&&(plotting(i) <=2)
pn=pn+1;       
ColorSet=varycolor(trajs);
figure (pn)
axis equal
hold on
for n=1:trajs
plot3(ysc(n,:,1)/Rj,ysc(n,:,2)/Rj,ysc(n,:,3)/Rj,'Color',ColorSet(n,:))%plot s/c trajectory
end

for n=1:nodes       %plot 
    if (event(n)==40)||(event(n)==41)||(((event(n)==17)||(event(n)==27))&&(coenp(2,n)>=1.0))
plot3(yP(n,1,1)/Rj,yP(n,1,2)/Rj,yP(n,1,3)/Rj,'ok','MarkerFaceColor','k')
%plot planets' path        
    else
plot3(yP(n,:,1)/Rj,yP(n,:,2)/Rj,yP(n,:,3)/Rj,'k') %plot planets' path
    end
end
            %----------Plot Locations of planets at JD0---------
    if (plotting(i)==1)
    plot3(rP0(1,1)/Rj,rP0(2,1)/Rj,rP0(3,1)/Rj,'oc','MarkerFaceColor','c')    
    %plot initial location of 1st planet
for n=2:nodes
    plot3(rP(1,n)/Rj,rP(2,n)/Rj,rP(3,n)/Rj,'or','MarkerFaceColor','r')         
    %plot location of planet at intersection    
    plot3(rP0(1,n)/Rj,rP0(2,n)/Rj,rP0(3,n)/Rj,'dm','MarkerFaceColor','m')       
    %plot initial location of planet     
end
    end
       %-----Plot SOI, North Poles, Planetocentric Orbit h-vectors----
    if plotting(i)==2
for n=1:nodes
        if (event(n)~=17)&&(event(n)~=27)&&(event(n)~=40)&&(event(n)~=41)
theta = -acos(npp(3,n));
psi = -atan2(npp(2,n),npp(1,n));
rot2=[cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)];
rot3 = [cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0; 0 0 1];
for j = 1:361
    xsoi = rsoi(n)*cosd(j);
    ysoi = rsoi(n)*sind(j);
    zsoi  = 0 ;
    Rsoi(1:3,j) = [xsoi;ysoi;zsoi];
    Rsoi(1:3,j) = rot3*rot2*Rsoi(1:3,j)/Rj;
end
%     plot3(Rsoi(1,:)+rP(n,1)/AU,Rsoi(2,:)+rP(n,2)/AU,Rsoi(3,:)+rP(n,3)/AU,'-r') %plot SOI
    patch(Rsoi(1,:)+rP(1,n)/Rj,Rsoi(2,:)+rP(2,n)/Rj,Rsoi(3,:)+rP(3,n)/Rj,'r','FaceAlpha',0.5); %patch object SOI
    quiver3(rP(1,n)/Rj,rP(2,n)/Rj,rP(3,n)/Rj,npp(1,n),npp(2,n),npp(3,n),rsoi(n)*2/Rj,'g') %plot North Pole
    quiver3(rP(1,n)/Rj,rP(2,n)/Rj,rP(3,n)/Rj,hhat(1,n),hhat(2,n),hhat(3,n),rsoi(n)*2/Rj,'k') %plot orbit ang. mom.
        end
end
    end
    
plot3(0,0,0,'or','MarkerFaceColor','y')         %plot sun
hold off
xlabel('x axis (Rj)')
ylabel('y axis (Rj)')
zlabel('z axis (Rj)')
    end
%--------------------Plot Heliocentric Trajectory Movie--------------------
    if plotting(i)==9
        
pn=pn+1;
% pts=500;
toft=sum(tof);
ppt=pts/toft;
% dt=0.1;
yscm=[];
tof;
for n=1:trajs
num=max([3 ceil(ppt*tof(n))]); %check that at least 3 pts are used in order 
                               % to match the defined frame rate
Zsc=[rP(:,n)' vd(:,n)'];
[tscm,yscm0]=ode45(@(t,y)EOMwrtbody(t,y,muj),linspace(0,tof(n)*86400,num),Zsc,options);
yscm=[yscm;yscm0];
end
pts=length(yscm);
for n=1:nodes
    if (event(n)==40)||(event(n)==41)||(((event(n)==17)||(event(n)==27))&&(coenp(2,n)>=1.0))
tPm(n,1)=0;
yPm(n,1,1:6)=[rP0(:,n)' vP0(:,n)'];        
    else
ZP=[rP0(:,n)' vP0(:,n)'];
[tPm(n,:),yPm(n,1:pts,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,muj),linspace(0,toft*86400,pts),ZP,options); 
    end
end

figure (pn)
axis equal
hold on
for n=1:trajs
plot3(ysc(n,:,1)/Rj,ysc(n,:,2)/Rj,ysc(n,:,3)/Rj,'k') %plot s/c trajectory path
end
for n=1:nodes
    if (event(n)==40)||(event(n)==41)||(((event(n)==17)||(event(n)==27))&&(coenp(2,n)>=1.0))
plot3(yP(n,1,1)/Rj,yP(n,1,2)/Rj,yP(n,1,3)/Rj,'ok','MarkerFaceColor','c')
%plot planets' path        
    else    
plot3(yP(n,:,1)/Rj,yP(n,:,2)/Rj,yP(n,:,3)/Rj,'b') %plot planets' path
    end
end
%adjust axis limits so all dots,lines,ect are within recordable area
xminmax=xlim;
yminmax=ylim;
zminmax=zlim;
axis([xminmax yminmax zminmax].*1.05);
for m=1:pts
SC=plot3(yscm(m,1)/Rj,yscm(m,2)/Rj,yscm(m,3)/Rj,'ok','MarkerFaceColor','k');        
%plot s/c trajectory    
    for n=1:nodes
        if (event(n)==40)||(event(n)==41)||(((event(n)==17)||(event(n)==27))&&(coenp(2,n)>=1.0))
NodeP(n)=plot3(yPm(n,1,1)/Rj,yPm(n,1,2)/Rj,yPm(n,1,3)/Rj,'ok','MarkerFaceColor','c');            
        else      
NodeP(n)=plot3(yPm(n,m,1)/Rj,yPm(n,m,2)/Rj,yPm(n,m,3)/Rj,'ok','MarkerFaceColor','r');
        end
    end
    
    if ~isempty(movie) %save frame to movie file if creating movie
frame=getframe;
writeVideo(writerObj,frame);
    else
pause(dt) %use pause to slow movie down but is "choppy"
    end
delete(SC)
    for n=1:nodes
delete(NodeP(n))        
    end
    
end

SC=plot3(yscm(end,1)/Rj,yscm(end,2)/Rj,yscm(end,3)/Rj,'ok','MarkerFaceColor','k');        
%plot s/c trajectory    
    for n=1:nodes
        if (event(n)==40)||(event(n)==41)||(((event(n)==17)||(event(n)==27))&&(coenp(2,n)>=1.0))
NodeP(n)=plot3(yPm(n,1,1)/Rj,yPm(n,1,2)/Rj,yPm(n,1,3)/Rj,'ok','MarkerFaceColor','c');            
        else  
NodeP(n)=plot3(yPm(n,end,1)/Rj,yPm(n,end,2)/Rj,yPm(n,end,3)/Rj,'ok','MarkerFaceColor','r');
        end
    end
hold off
    end
    
%--------------------Plot Planetary Departure Trajectory-------------------
    elseif (plotting(i) >= 10)&&(plotting(i) < 20)
npp=0;
for n=1:nodes
 if (event(n) >= 10)&&(event(n) <= 14)
    if (event(n) >= 10)&&(event(n) <= 13)
        npp=npp+1;
Zpp(n,:)=[rm(:,n)' vin(:,n)']; %integrate parking orbit
[tpp(npp,:),ypp(npp,1:nptsparkd,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),linspace(0,Tp(n),nptsparkd),Zpp(n,:),options);
    end
     if (event(n) == 14)
         npp=npp+1;
 [rmlv, vmlv]=orbel2rv(mupl(n),coep(1,n),coep(2,n),coep(3,n),coep(4,n),...
     coep(5,n),coep(6,n));        
Zpp(n,:)=[rmlv' vmlv']; %integrate parking orbit
[tpp(npp,:),ypp(npp,1:nptsparkd,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),linspace(0,Tp(n),nptsparkd),Zpp(n,:),options);

    end
  
Zp(n,:)=[rm(:,n)' vout(:,n)']; %integrate departure trajectory
[tp(npp,:),yp(npp,1:nptsSCd,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),logspace(0,log10(tsoiout(n)+1),nptsSCd),Zp(n,:),options);

hhdhat=cross(rm(:,n),vout(:,n))/norm(cross(rm(:,n),vout(:,n)));
theta = -acos(hhdhat(3));
psi = -atan2(hhdhat(2),hhdhat(1));
rot1=[1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
rot2=[cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)];
rot3 = [cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0; 0 0 1];
for j = 1:361
    xsoi = rsoi(n)*cosd(j);
    ysoi = rsoi(n)*sind(j);
    zsoi  = 0 ;
    Rsoi(1:3,j) = [xsoi;ysoi;zsoi];
    Rsoi(1:3,j) = rot3*rot2*Rsoi(1:3,j);
end
vinfplhat(:,npp)=vinfpl(:,n)/norm(vinfpl(:,n));
vPpplhat(:,npp)=vPplhat(:,n);
vpplhat(:,npp)=vplhat(:,n);
dvhat(:,npp)=dv(:,n)/norm(dv(:,n));
scale(npp) = rsoi(n);
rpld(npp)=rpl(n);
 end
end
[xpl,ypl,zpl]=sphere(50);
for j=1:npp
    pn=pn+1;
figure (pn)
axis equal
hold on
plot3(ypp(j,:,1), ypp(j,:,2), ypp(j,:,3), 'k')
plot3(yp(j,1,1), yp(j,1,2), yp(j,1,3), 'ok')
plot3(yp(j,:,1), yp(j,:,2), yp(j,:,3), 'b')
plot3(Rsoi(1,:),Rsoi(2,:),Rsoi(3,:),'--r')
quiver3(0,0,0,vinfplhat(1,j),vinfplhat(2,j),vinfplhat(3,j),scale(j),'--r')
quiver3(0,0,0,vPpplhat(1,j),vPpplhat(2,j),vPpplhat(3,j),scale(j),'k')
quiver3(0,0,0,vpplhat(1,j),vpplhat(2,j),vpplhat(3,j),scale(j),'c')
quiver3(yp(j,1,1),yp(j,1,2),yp(j,1,3),dvhat(1,j),dvhat(2,j),dvhat(3,j),scale(j)/4,'g')
surf(xpl*rpld(j),ypl*rpld(j),zpl*rpld(j))
hold off
legend('Departure Parking Orbit','Periapsis','Hyperbolic Trajectory','R_{SOI}',...
        'V_{\infty}-dir','V_{planet}-dir','V_{helio}-dir','\DeltaV')
xlabel('x axis (km)')
ylabel('y axis (km)')
zlabel('z axis (km)')
end
%--------------------Plot Planetary Arrival Trajectory---------------------
    elseif (plotting(i) >= 20)&&(plotting(i) < 30)
npp=0;
for n=1:nodes
  if (event(n) >= 20)&&(event(n) <= 24)
      if (event(n) >= 20)&&(event(n) <= 23)
        npp=npp+1;
Zpp(n,:)=[rm(:,n)' vout(:,n)']; %integrate parking orbit
[tpp(npp,:),ypp(npp,1:nptsparka,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),linspace(0,Tp(n),nptsparka),Zpp(n,:),options);
      end
      if (event(n) == 24)
         npp=npp+1;
  [rmlv, vmlv]=orbel2rv(mupl(n),coep(1,n),coep(2,n),coep(3,n),coep(4,n),...
      coep(5,n),coep(6,n));        
Zpp(n,:)=[rmlv' vmlv']; %integrate parking orbit
[tpp(npp,:),ypp(npp,1:nptsparka,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),linspace(0,Tp(n),nptsparka),Zpp(n,:),options);        
      end
Zp(n,:)=[rm(:,n)' vin(:,n)']; %integrate arrival trajectory
[tp(npp,:),yp(npp,1:nptsSCa,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),-logspace(0,log10(tsoiin(n)+1),nptsSCa),Zp(n,:),options);

hhdhat=cross(rm(:,n),vin(:,n))/norm(cross(rm(:,n),vin(:,n)));
theta = -acos(hhdhat(3));
psi = -atan2(hhdhat(2),hhdhat(1));
rot1=[1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
rot2=[cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)];
rot3 = [cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0; 0 0 1];
for j = 1:361
    xsoi = rsoi(n)*cosd(j);
    ysoi = rsoi(n)*sind(j);
    zsoi  = 0 ;
    Rsoi(1:3,j) = [xsoi;ysoi;zsoi];
    Rsoi(1:3,j) = rot3*rot2*Rsoi(1:3,j);
end
vinfplhat(:,npp)=vinfpl(:,n)/norm(vinfpl(:,n));
vPpplhat(:,npp)=vPplhat(:,n);
vpplhat(:,npp)=vplhat(:,n);
dvhat(:,npp)=dv(:,n)/norm(dv(:,n));
scale(npp) = rsoi(n);
rpla(npp)=rpl(n);
rvinf(:,npp)=-vinfplhat(:,npp)*scale(npp);
    end
end
[xpl,ypl,zpl]=sphere(50);

for j=1:npp
    pn=pn+1;
figure (pn)
axis equal
hold on
plot3(ypp(j,:,1), ypp(j,:,2), ypp(j,:,3), 'k')
plot3(ypp(j,1,1), ypp(j,1,2), ypp(j,1,3), 'ok')
plot3(yp(j,:,1), yp(j,:,2), yp(j,:,3), 'b')
plot3(Rsoi(1,:),Rsoi(2,:),Rsoi(3,:),'--r')
quiver3(rvinf(1,j),rvinf(2,j),rvinf(3,j),vinfplhat(1,j),vinfplhat(2,j),vinfplhat(3,j),scale(j)*0.75,'--r')
quiver3(0,0,0,vPpplhat(1,j),vPpplhat(2,j),vPpplhat(3,j),scale(j),'k')
quiver3(0,0,0,vpplhat(1,j),vpplhat(2,j),vpplhat(3,j),scale(j),'c')
quiver3(ypp(j,1,1),ypp(j,1,2),ypp(j,1,3),dvhat(1,j),dvhat(2,j),dvhat(3,j),scale(j)/4,'g')
surf(xpl*rpla(j),ypl*rpla(j),zpl*rpla(j))
quiver3(rvinf(1,j),rvinf(2,j),rvinf(3,j),vinfplhat(1,j),vinfplhat(2,j),vinfplhat(3,j),scale(j),'--r','ShowArrowHead','off')

hold off
legend('Arrival Parking Orbit','Periapsis','Hyperbolic Trajectory','R_{SOI}',...
    'V_{\infty}-dir','V_{planet}-dir','V_{helio}-dir','\DeltaV')
xlabel('x axis (km)')
ylabel('y axis (km)')
zlabel('z axis (km)')
end
%--------------------Plot Planetary Flyby Trajectory-----------------------
    elseif (plotting(i) >= 30)&&(plotting(i) < 40)
nfb=0;
for n=1:nodes
    if (event(n) == 30)||(event(n) == 31)||(event(n) == 33)||(event(n) == 35)
        nfb=nfb+1;
Zfb=[rpin(:,n)' vpin(:,n)'];
[tfb(nfb,:),yfbin(nfb,1:nptsflyby,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),-logspace(0,log10(tsoiin(n)+1),nptsflyby),Zfb,options);
Zfb=[rpout(:,n)' vpout(:,n)'];
[tfb(nfb,:),yfbout(nfb,1:nptsflyby,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),logspace(0,log10(tsoiout(n)+1),nptsflyby),Zfb,options);
Zfb=[rm(:,n)' vin(:,n)'];
[tfb(nfb,:),yfb1(nfb,1:nptsflyby,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),-logspace(0,log10(tsoiin(n)+1),nptsflyby),Zfb,options);
Zfb=[rm(:,n)' vout(:,n)'];
[tfb(nfb,:),yfb2(nfb,1:nptsflyby,1:6)]=ode45(@(t,y)EOMwrtbody(t,y,mupl(n)),logspace(0,log10(tsoiout(n)+1),nptsflyby),Zfb,options);
rnfb(:,nfb)=rm(:,n);
vinfinhatfb(:,nfb)=vinfinhat(:,n);
vinfouthatfb(:,nfb)=vinfouthat(:,n);
dvhat(:,nfb)=dv(:,n)/norm(dv(:,n));
scale(nfb)=rsoi(n);
rplf(nfb)=rpl(n);
rvinf(:,nfb)=-vinfinhatfb(:,nfb)*scale(nfb);
    end
end
[xpl,ypl,zpl]=sphere(50);

for j=1:nfb
    pn=pn+1;
figure (pn)
axis equal
hold on
plot3(yfbin(j,:,1), yfbin(j,:,2), yfbin(j,:,3), ':r')     
%plot theoretical incoming hyp trajectory
plot3(yfbout(j,:,1), yfbout(j,:,2), yfbout(j,:,3), ':b')  
%plot theoretical outgoing hyp trajectory
plot3(yfb1(j,:,1), yfb1(j,:,2), yfb1(j,:,3), 'k')        
%plot s/c incoming trajectory
plot3(yfb2(j,:,1), yfb2(j,:,2), yfb2(j,:,3), 'c')        
%plot s/c outgoing trajectory
surf(xpl*rplf(j),ypl*rplf(j),zpl*rplf(j))%plot gravitational body
quiver3(rnfb(1,j),rnfb(2,j),rnfb(3,j),dvhat(1,j),dvhat(2,j),dvhat(3,j),scale(j)/4,'g')
quiver3(rvinf(1,j),rvinf(2,j),rvinf(3,j),vinfinhatfb(1,j),vinfinhatfb(2,j),vinfinhatfb(3,j),scale(j)*.75,'--r')
quiver3(0,0,0,vinfouthatfb(1,j),vinfouthatfb(2,j),vinfouthatfb(3,j),scale(j),'--b')
quiver3(rvinf(1,j),rvinf(2,j),rvinf(3,j),vinfinhatfb(1,j),vinfinhatfb(2,j),vinfinhatfb(3,j),scale(j),'--r','ShowArrowHead','off')

hold off
legend('Theoretical Incoming Hyperbolic Leg','Theoretical Outgoing Hyperbolic Leg',...
    'Incoming Flyby Trajectory','Outgoing Flyby Trajectory','Gravitational Body',...
    '\DeltaV Dir.','V_{\infty}^{in} Dir.','V_{\infty}^{out} Dir.')
xlabel('x axis (km)')
ylabel('y axis (km)')
zlabel('z axis (km)')
end
    elseif plotting(i) ==40
        
    end
end
if ~isempty(movie) %close movie file if created movie
close(writerObj);
end
end
