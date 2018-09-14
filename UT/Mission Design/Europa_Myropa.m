%% ************************************************************************
%----------------------------Cassini Trajectory Driver---------------------
%**************************************************************************
%                       ______ ____   ___    ______ ______
%                      /_  __// __ \ /   |  / ____//_  __/
%                       / /  / /_/ // /| | / /      / /   
%                      / /  / _, _// ___ |/ /___   / /    
%                     /_/  /_/ |_|/_/  |_|\____/  /_/                                      
clear all
close all
clc
addpath('Bin')
Rj=71492;  %(km) provides generally accepted Rj measurement for convenience 
muj=1.26712764e+008;
excelfile=mfilename; %saves the name of this trajectory file for naming the Excel Output file
savefile=[]; %Saves the case variables to MAT-file with desired file name (use ' ')
             %Can use mfilename if appropriate (keep empty if not saving file)
loadfile=[]; %Loads the case variables to MAT-file with desired file name (use ' ')
             %Can use mfilename if appropriate (keep empty if not loading file)
             %Only loads and overwrites optimizing variables (JD0,tof,rnb)
             %NOTE: Make sure maneuvers and other settings below are chosen 
             %to agree with the loaded file (same # of nodes, revs, ect.)
ephtype=1;      %Choose one of the following ephemeris types
                %1: JPL's Horizons ephemeris (basic kepler propagation)
                

writeout=0;     %Choose to create an output excel file with all the results
                %of the optimization routine
                %0: Do not write file
                %1: Create Excel file with results
                
%----------------------Optimization Parameters-----------------------------
optimization=3;   %Choose to optimize the initial node's julian date and time's
                  %of flight between nodes or simply determine trajectory
                  %for the specified fixed times.
                  %Fixed parameters (optimization iterations disabled): 0
            %OPTIMIZATION TOOLBOX OPTIONS:      
                  %Optimize Unconstrained: 1 (fminunc: gradient method, numerical derivatives)
                  %Optimize Unconstrained: 2 (fminsearch: simplex method, gradient-free)
                  %Optimize Nonlinear Con: 3 (fmincon: gradient method, numerical derivatives)
                  %Optimize Nonlinear Con: 4 (fminsearchcon: simplex method, gradient-free)
                  %Optimize Unconstrained: 5 (fsolve: gradient method, numerical derivatives)
                  %Optimize L/U Bound Con: 6 (lsqnonlin: gradient method, numerical derivatives)
                  %Optimize Nonlinear Con: 7 (fminimax: gradient method, numerical derivatives)
            %GLOBAL OPTIMIZATION TOOLBOX OPTIONS:    
                  %Optimize Nonlinear Con: 10 (patternsearch: directed search, gradient-free)  
                  %Optimize L/U Bound Con: 11 (simulannealbnd: stochastic search, gradient-free)  
                  %Optimize Nonlinear Con: 12 (ga: genetic algorithm, stochastic search, gradient-free)
                        %See AI:6 for applying constraints to the procedure
                  
                  %Use the following optoptions for fminunc Optimization [1]     
% optoptions=optimoptions('fminunc','Display','iter','MaxFunEvals',9000,'MaxIter',9000,...
%                      'TolFun',1e-6,'TolX',1e-8,'Algorithm','quasi-newton');
                        
                  %Use the following optoptions for fminsearch Optimization [2,4]     
% optoptions=optimset('Display','iter','MaxFunEvals',3000,'MaxIter',3000,...
%                      'TolFun',1e-6,'TolX',1e-8);

                  %Use the following optoptions for fmincon Optimization [3]
                  %Choose 'interior-point', 'sqp', or 'active-set' algorithm
                  %Use active-set algorithm with 'RelLineSrchBnd',## to avoid
                  %large steps
                  %Can add in 'UseParallel','always'
optoptions=optimoptions('fmincon','Display','iter','MaxFunEvals',9000,'MaxIter',9000,...
                    'Algorithm','active-set','TolFun',1e-6,'TolX',1e-8,...
                    'TolCon',1e-8);             

                  %Use the following optoptions for fsolve Optimization [5]  
                  %Choose 'trust-region-dogleg','trust-region-reflective',or 'levenberg-marquardt'
% optoptions=optimoptions('fsolve','Display','iter','MaxFunEvals',9000,'MaxIter',9000,...
%                      'TolFun',1e-6,'TolX',1e-8,'Algorithm','levenberg-marquardt');
                 
                  %Use the following optoptions for lsqnonlin Optimization [6]  
                  %Choose 'trust-region-reflective', 'levenberg-marquardt'
% optoptions=optimoptions('lsqnonlin','Display','iter','MaxFunEvals',9000,'MaxIter',9000,...
%                      'TolFun',1e-6,'TolX',1e-8,'TolCon',1e-8,'Algorithm','levenberg-marquardt');                 
               
                  %Use the following optoptions for fminimax Optimization [7]     
% optoptions=optimoptions('fminimax','Display','iter','MaxFunEvals',9000,'MaxIter',9000,...
%                      'TolFun',1e-6,'TolX',1e-8,'TolCon',1e-8,'MeritFunction','singleobj');


                 %Use the following optoptions for patternsearch Optimization [10]     
                 %Can add in: 'SearchMethod',@searchneldermead or
                 %@searchga or @searchlhs or @GSSPositiveBasisNp1
                 %@MADSPositiveBasisNp1
                 %Can add in: ,'UseParallel','always'
% optoptions=psoptimset('Display','iter','MaxFunEvals',9000,'MaxIter',9000,...
%                      'TolFun',1e-6,'TolX',1e-8,'TolCon',1e-8,'TimeLimit',100); 

                 %Use the following optoptions for simulannealbnd Optimization [11]  
                 %Can add in: 'HybridFcn',@fminsearch or @patternsearch or
                 %@fminunc or @fmincon 
                 %Can add in: ,'InitialTemperature', 100 
% optoptions=saoptimset('Display','iter','MaxFunEvals',9000,'MaxIter',9000,...
%                      'TolFun',1e-6,'TimeLimit',100); 
                 
                 %Use the following optoptions for ga (genetic algorithm) Optimization [12]  
                 %Can add in: 'HybridFcn',@fminsearch or @patternsearch or
                 %@fminunc or @fmincon 
                 %Genetic algorithm options: 'CreationFcn','CrossoverFcn','CrossoverFraction', 
                 %'MigrationFraction','MigrationInterval','MigrationDirection',
                 %'MutationFcn','NonlinConAlgorithm','SelectionFcn',''
                 %Can add in: ,'UseParallel','always'
% optoptions=gaoptimset('Display','iter','Generations',400,...
%                      'TolFun',1e-8,'TolCon',1e-8,'TimeLimit',100); 

                  %See AI:6 for additional optimization options  
                  
%-------------------------Plotting Parameters------------------------------
plotting=[9];   %Choose preconfigured Plots to be generated
                 %0:    Heliocentric trajectory with planets' orbits
                 %1:    Heliocentric trajectory with planets' orbits and
                 %      positions at date of first node
                 %2:    Heliocentric trajectory with planets' orbits with
                 %      spheres of influence, North Poles, and angular
                 %      momentum dir of planetary orbits/encounters
                 %9:    Heliocentric trajectory movie with planets
                 %10:   Departure trajectory(s) at node(s)
                 %20:   Arrival trajectory(s) at node(s)
                 %30:   Flyby trajectory(s) at node(s) 
                 
 pts=200;   %number of pts to plot [~200-2000] pts
 dt=0.1;    %delay between movie frames [~0.1-1] sec
 movie=['orbitz'];  %Saves Heliocentric movie file with desired file name (use ' ')
            %Can use mfilename if appropriate (keep empty if not saving file)
 fps=20;    %Adjusts the saved movie frame rate (typically 20-30)          
            
%% -------------------------General Trajectory Parameters--------------------
npl = [3 3 3 3 4 3 4];    %moon's number at which each trajectory occurs 
         %(1=Io, 2=Europa, 3=Ganymede, 4=Callisto)
         %0=spaceburn/midcourse maneuver/Defined-Orbit)
                    
nn = length(npl);   %(Do not change)number of interplanetary nodes    

event(1:nn)=[31 31 31 31 31 31 24]; % Trajectory Event at each node:
                     % 10: Departure from a parking orbit about planet 
                     %     (optimization variables: i,o,w), i will be >=
                     %     trajectory declination, see AI:1
                     % 11: Departure from a parking orbit about planet 
                     %     (optimization variables: o,w), i is fixed, see AI:1
                     % 12: Departure from a parking orbit about planet 
                     %     (optimization variables: ta), i,o,w are fixed, see AI:1
                     % 13: Departure from a parking orbit about planet 
                     %     (optimization variables: None), fixed orbit, see AI:1
                     % 14: Departure from a launch vehicle trajectory  
                     %     (optimization variables: o,w), vinfav and i are 
                     %      fixed, see AI:1
                     % 17: Depart user defined body(with negligible
                     %     gravity) or heliocentric orbit, see AI:4
                     % 18: Depart from a parking orbit about user defined body (not
                     %     available yet)
                     % 20: Arrival to a parking orbit about planet
                     %     (optimization variables: i,o,w), i will be >=
                     %     trajectory declination 
                     % 21: Arrival to a parking orbit about planet 
                     %     (optimization variables: o,w), i is fixed 
                     % 22: Arrival to a parking orbit about planet 
                     %     (optimization variables: ta), i,o,w are fixed, see AI:1
                     % 23: Arrival to a parking orbit about planet 
                     %     (optimization variables: None), fixed orbit, see AI:1
                     % 24: Arrival to entry trajectory
                     %     (optimization variables: o,w), vinfav and i are 
                     %      fixed, see AI:1
                     % 27: Rendezvous with user defined body(with negligible
                     %     gravity) or heliocentric orbit, see AI:4
                     % 28: Arrival to a parking orbit about user defined body (not
                     %     available yet)
                     % 30: Unpowered gravity assist, see AI:2 & 6 
                     %      Requires nonlinear constraint vinfmatch set
                     %      for each unpowered gravity assist node
                     % 31: Unpowered gravity assist (DV at Rinf), see AI:2
                     %      Converges upon an unpowered gravity assist
                     %      without any constraints
                     % 33: Powered gravity assist (optimal), see AI:2 
                     %      Utilizes DV impulse at optimal location during
                     %      the gravity assist encounter
                     % 35: Powered gravity assist (@ periapse), see AI:2
                     %      Utilizes DV impulse at periapse of encounter
                     % 40: Midcourse Maneuver/(heliocentric Delta-V) 
                     %     Initial position guess provided by program,
                     %     If multiple consecutive midcourse maneuvers,
                     %     additional revs strictly to improve guess 
                     %     at AI:3 (Additional Revs...)
                     %     CANNOT BE FIRST OR LAST NODE!
                     %     CANNOT use consecutively with Space Burn
                     % 41: Space Burn/(heliocentric Delta-V) for ANY NODE
                     %     Initial position guess provided by user at
                     %     AI:3 (Non-body Node Location) 
                     %     CANNOT use consecutively with Midcourse Maneuver
                     
JD0 = date2jed([2028 4 4 12 0 0]); %initial Julian Date (yr,mon,d,h,min,s)
TG=30.42;
tof(1:nn-1) = [6.7*TG 1.9*TG 0.9*TG 0.6*TG 0.8*TG 0.5*TG];    %time of flight between nodes (days)
grade(1:nn-1) = [0 0 0 0 0 0];   %directio      n of trajectory between nodes 
                                 %(typically prograde)
                            %Prograde: 0, Retrograde: 1
rev(1:nn-1)= [0 0 0 0 0 0]; %number of orbital revolutions for trajectory between nodes
                            %rev value can be negative for left branch
                            %solution and positive for right branch
branch=0;              %Identify if use particular branch configuration above
                       %by setting to 0, or set to 1 if exploring ALL
                       %permutations (executes multiple optimization runs!)
                       
branchopt=1;           %Set the number of iterations to be run for each branch
                       %when evaluating branch permutations. 
                       %1 iteration, 10, 100, etc.
                       %Set to 0 for fully optimized permutation cases
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%% ************************************************************************
%**************************************************************************
%----------------------------Additional Input (AI)-------------------------
%**************************************************************************
%**************************************************************************


%% AI:1---------------------Departure and Arrival Parameters---------------
                                       
%      ___                      _        __    _             _           
%     /   \___ _ __   __ _ _ __| |_     / /   /_\  _ __ _ __(_)_   _____ 
%    / /\ / _ \ '_ \ / _` | '__| __|   / /   //_\\| '__| '__| \ \ / / _ \
%   / /_//  __/ |_) | (_| | |  | |_   / /   /  _  \ |  | |  | |\ V /  __/
%  /___,' \___| .__/ \__,_|_|   \__| /_/    \_/ \_/_|  |_|  |_| \_/ \___|
%             |_|                                                                                                                                                                                                                                                                                                           
                                      %%%
                                       %
coep(1:6,1:nn)=0; %(do not change)initialize classical orbital elements for each node
altp(1:nn)=0; %(do not change)initialize theoretical periapsis altitude
ilv(1:nn)=0; %(do not change)initialize launch inclination
vinfavmag(1:nn)=0; %(do not change)initialize vinf available from launch
dvopt(1:nn)=0; %(do not change)initialize delta-V optimization criteria

%Define departure and arrival parking orbital elements at intended nodes
%If additional departure/arrival parking orbits are required (such as for
%mid-trajectory/mission stopover at a planet) simply Copy & Paste the orbit
%information for a single node as many times as necessary and specify the
%appropriate orbital elements and intended node for your mission
%--------------------------------------------------------------------------

%+++++++++++++++++++++++++++Single Node Launch Information+++++++++++++++
%Departure launch elements
% node=1;     %particular node for launch elements
% altp(node) = 170; % altitude of periapsis (km),usually 185-200km
% ilv(node) = 28.5;    %latitude of launch site (deg), inclination of optimum launch trajectory
% vinfavmag(node) = sqrt(16.6); %hyperbolic veloicty at infinity available from launch vehicle
% dvopt(node)=1;    %defines if excess vinf is acceptable (no delta-V) or detrimental (delta-V)
            %vinfreq=vinfav: 0
            %vinfreq<=vinfav:1
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++Single Node Orbit Information++++++++++++++++
%Departure parking orbital elements
% node=1;      %particular node for orbital elements
% 
% ap = 10000; %semimajor axis (km)
% ep = 0.1;    %eccentricity
% ip = 60;     %inclination (deg)
% op = 0;      %right ascension of ascending node (deg)
% wp = 0;      %argument of periapsis (deg)
% tap = 0;     %true anomaly (deg)
% coep(1:6,node) = [ap ep ip op wp tap]; %classical orbital elements at particular node
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%+++++++++++++++++++++++++++Single Node Orbit Information++++++++++++++++
%Arrival parking orbital elements
% node=nn;      %particular node for orbital elements
% 
% ap = 4585959; %semimajor axis (km)
% ep = 0.98239;    %eccentricity
% ip = 23.534;     %inclination (deg)
% op = 200;      %right ascension of ascending node (deg)
% wp = 75;      %argument of periapsis (deg)
% tap = 10;     %true anomaly (deg)
% coep(1:6,node) = [ap ep ip op wp tap]; %classical orbital elements at particular node
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++Single Node Launch Information+++++++++++++++
%Departure launch elements
node=nn;     %particular node for launch elements
altp(node) = 100; % altitude of periapsis (km),usually 185-200km
ilv(node) =90;    %latitude of launch site (deg), inclination of optimum launch trajectory
vinfavmag(node) = 10; %hyperbolic veloicty at infinity available from launch vehicle
dvopt(node)=1;    %defines if excess vinf is acceptable (no delta-V) or detrimental (delta-V)
            %vinfreq=vinfav: 0
            %vinfreq<=vinfav:1
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%% AI:2----------------Gravity Assist (Flyby) Parameters-------------------
                                                            
%        ___                 _ _             _           _     _       
%       / _ \_ __ __ ___   _(_) |_ _   _    /_\  ___ ___(_)___| |_ ___ 
%      / /_\/ '__/ _` \ \ / / | __| | | |  //_\\/ __/ __| / __| __/ __|
%     / /_\\| | | (_| |\ V /| | |_| |_| | /  _  \__ \__ \ \__ \ |_\__ \
%     \____/|_|  \__,_| \_/ |_|\__|\__, | \_/ \_/___/___/_|___/\__|___/
%                                  |___/                                                           
                                      %%%
                                       %
rfbc(1:nn)=1.1; %(do not change)initialize flyby radius vector
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		%Specify specific node's gravity assist closest approach (periapse)
		%radius mag. with respect to planet, otherwise the above default value
		%will be used.

% node=2;       %particular node for specified flyby radius
% rV=6051.8; %km
% rfbc(node)=(rV+287)/rV; %flyby radius magnitude (radii of particular planet)
% %
% node=4;       %particular node for specified flyby radius
% rV=6051.8; %km
% rfbc(node)=(rV+603)/rV; %flyby radius magnitude (radii of particular planet)
% %
% node=5;       %particular node for specified flyby radius
% rE=6378.1366; %km
% rfbc(node)=(rE+1171)/rE; %flyby radius magnitude (radii of particular planet)
% %
% node=6;       %particular node for specified flyby radius
% rJ=71492; %km
% rfbc(node)=(rJ+9720000)/rJ; %flyby radius magnitude (radii of particular planet)
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%++++++++++++++++++++++++++Gravity Assist Optimization+++++++++++++++++++++
          
optr(1:nn)=1; %(do not change)initialize flyby radius optimization            
            %Choose to vary the flyby periapsis radius distance to better 
            %optimize the Delta-V
            %Use fixed rpmag value: 0 (integer)
            %Vary rpmag value: 1 (integer)
            
% node=2;
% optr(node)=1;
% 
% node=4;
% optr(node)=1;
% 
% node=5;
% optr(node)=1;
% 
% node=6;
% optr(node)=1;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%% AI:3---------------Define Non-body (MCM/DSM) Node Location--------------  
%                            ___             __     ___  __         
%                   /\/\    / __\ /\/\      / /    /   \/ _\  /\/\  
%                  /    \  / /   /    \    / /    / /\ /\ \  /    \ 
%                 / /\/\ \/ /___/ /\/\ \  / /    / /_// _\ \/ /\/\ \
%                 \/    \/\____/\/    \/ /_/    /___,'  \__/\/    \/
                                                  
                                      %%%
                                       %
rnb(1:3,1:nn)=0; %(do not change)initialize non-body node position vector
revig(1:nn)=0; %(do not change)initialize total trajectory rev initial guess 
                %through non-body node
%--------------------------------------------------------------------------

%++++++Additional Revs for Consecutive Maneuvers (only for event:40)+++++++
               %Define a node's "revig" function to be an integer of additional 
               %orbital revolutions before arriving at the next node
% node=2;
% revig(node)=1; %Trajectory rev through non-body node      
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%+++++++++++++++++Node Location Information (only for event:41)++++++++++++
node=nn-1;
rnb(1:3,node)=[100 100 0]; 
%non-body position in rectangular coord. (Rj) at node     
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%% AI:4-------------------Rendezvous with Defined Orbit--------------------

%                __                _                                
%               /__\ ___ _ __   __| | ___ ______   _____  _   _ ___ 
%              / \/// _ \ '_ \ / _` |/ _ \_  /\ \ / / _ \| | | / __|
%             / _  \  __/ | | | (_| |  __// /  \ V / (_) | |_| \__ \
%             \/ \_/\___|_| |_|\__,_|\___/___|  \_/ \___/ \__,_|___/
                                                      
                                      %%%
                                       %
coenp(1:8,1:nn)=0; %initialize non-planet (heliocentric) orbit
%--------------------------------------------------------------------------

%+++++++++++++++++++++++++++Single Node Orbit Information++++++++++++++++++
%Heliocentric orbital elements
node=1;      %particular node for orbital elements
rp = 12.96*Rj;
ra = 265.1*Rj;
ap = (rp+ra)/2; %semimajor axis (km)
ep = 1-rp/ap;    %eccentricity
ip = 5.4;     %inclination (deg)
op = 0;      %right ascension of ascending node (deg)
wp = 0;      %argument of periapsis (deg)
tap = 0;     %true anomaly (deg)
Tp=2*pi*sqrt(ap^3/muj)/86400;
JDp = JD0+Tp/2; %reference juilian date of orbital elements
mu = 0; %gravitational parameter (km^2/s^3)
coenp(1:8,node) = [ap ep ip op wp tap JDp mu]; %classical orbital elements at particular node
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%% AI:6-----------------------Optimization Options-------------------------

%            ___       _   _           _          _   _             
%           /___\_ __ | |_(_)_ __ ___ (_)______ _| |_(_) ___  _ __  
%          //  // '_ \| __| | '_ ` _ \| |_  / _` | __| |/ _ \| '_ \ 
%         / \_//| |_) | |_| | | | | | | |/ / (_| | |_| | (_) | | | |
%         \___/ | .__/ \__|_|_| |_| |_|_/___\__,_|\__|_|\___/|_| |_|
%               |_|                                                 

                                      %%%
                                       %
%++++++++++++++++++++++++++++++Optimization Nodes++++++++++++++++++++++++++
optnode(1:nn)=1; %(do not change)initialize optimization nodes vector
                 %Specify the weight for any node. If you wish a node's DV
                 %to be ignored in the optimization process, define that 
				 %node's "optnode" parameter to equal zero, such as:
%                       node=1;
%                       optnode(node)=0;
				%alternatively 
%                       optnode(1)=0;   
                
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%+++++++++++++++++++++++++++++Optimization Scaling+++++++++++++++++++++++++
optscale=0; %(do not change) initialize optimization scaling 
                %Specify which method is used for normalizing the
                %optimization variables (JD0, tofs, MCMs)
                %0: No scaling, variables are iterated in natural units
                %1: Scale variables by initial guess value 
                %2: Scale times by initial total tof (MCMs scaled as in [1])
                %3: Scale times by avg tof (MCMs scaled as in [1])
%                     optscale=3;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%+++++++++++++++++++++Optimization Variable MaxStep++++++++++++++++++++++++
maxstepJD0=0; %(do not change) initialize maxstep for initial julian date
maxsteptof(1:nn-1)=0; %(do not change) initialize maxstep for each tof
maxstepMCM(1:3,1:nn)=0; %(do not change) initialize maxstep for each MCM (X,Y,Z)
                %Specify a nonzero maxstep for JD0
%             maxstepJD0=0.9; %(days)

                %Specify a nonzero maxstep for tofs
%             leg=1;    
%             maxsteptof(leg)=1; %(days)
%             maxsteptof(3)=0.1; %(days)
%             maxsteptof(:)=1; %(days) for ALL legs

                %Specify a nonzero maxstep for MCMs (X,Y,Z)
%             node=3;    
%             maxstepMCM(:,node)=[0.1 0.1 0.1]; %(Rj)
%             maxstepMCM(1,5)=0.2; %(Rj)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%+++++++++++++++++++++++Optimization Constraints+++++++++++++++++++++++++++
tofTOTUB=[];%(do not change)initialize Total TOF lower bound
tofTOTLB=[];%(do not change)initialize Total TOF upper bound
tofTOTEQ=[];%(do not change)initialize Total TOF equality constraint
tofUB=[];%(do not change)initialize each TOF lower bound
tofLB=[];%(do not change)initialize each TOF upper bound
tofEQ=[];%(do not change)initialize each TOF equality constraint
JD0UB=[]; %(do not change)initialize Initial JD upper bound
JD0LB=[]; %(do not change)initialize Initial JD lower bound
JDUB=[]; %(do not change)initialize Final JD upper bound
JDLB=[]; %(do not change)initialize Final JD lower bound
JDEQ=[]; %(do not change)initialize Final JD equality constraint
MCMUB=[]; %(do not change)initialize  MCM/DSM location (X,Y,Z) upper bound
MCMLB=[]; %(do not change)initialize  MCM/DSM location (X,Y,Z) lower bound
MCMEQ=[]; %(do not change)initialize  MCM/DSM location (X,Y,Z) equality constraint
MCMRUB=[];%(do not change)initialize  MCM/DSM location radius mag. upper bound
MCMRLB=[];%(do not change)initialize  MCM/DSM location radius mag. lower bound
MCMREQ=[];%(do not change)initialize  MCM/DSM location radius mag. equality constraint
DVUB=[]; %(do not change)initialize DV mag. upper bound
DVLB=[]; %(do not change)initialize DV mag. lower bound
DVEQ=[]; %(do not change)initialize DV mag. equality constraint
vinfinUB=[]; %(do not change)initialize vinfin mag. upper bound
vinfinLB=[]; %(do not change)initialize vinfin mag. lower bound
vinfinEQ=[]; %(do not change)initialize vinfin mag. equality constraint
vinfoutUB=[]; %(do not change)initialize vinfout mag. upper bound
vinfoutLB=[]; %(do not change)initialize vinfout mag. lower bound
vinfoutEQ=[]; %(do not change)initialize vinfout mag. equality constraint
vinfmatch=[]; %(do not change)initialize vinf in/out mag. equality constraint

% Choose each particular constraint  to be implemented for a parameter
% Lower and upper bounds can be set for each paramter (and individual
% trajectory leg tofs, MCMs, or DVs). When an equality constraint  
% is used for a parameter DO NOT also implement an upper OR lower bound; 
% such conflicting constraints impairs the optimization process. 

%##########################################################################
%#######################----UPPER/LOWER BOUNDS----#########################
%##########################################################################

        %Specify each leg's TOF constraint upper/lower bounds OR equality
        %Define the upper/lower bound vectors or equality vector where
        %each element represents that leg's tof constraint
        %Any zero tof will be ignored (only nonzero values are constraints)
%                       leg = 1;
%                     tofUB(leg)= 500; %(days)
%                       leg = 3;
%                     tofLB(leg)= 300; %(days)
%                       leg = 5;
%                     tofEQ(leg)= 123.456; %(days)
        %alternatively:
%                     tofUB(1)= 500; %(days)
%                     tofLB(3)= 300; %(days)
%                     tofEQ(5)= 123.456; %(days)
% tofLB(1)=90;
% tofUB(1)=250;
% tofLB(2)=150;
% tofUB(2)=300;
% tofLB(3)=150;
% tofUB(3)=300;
% tofLB(4)=10;
% tofUB(4)=90;
% tofLB(5)=450;
% tofUB(5)=600;
% tofLB(6)=1000;
% tofUB(6)=1300;
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        %Specify Initial JD constraints upper/lower bounds
%                     JD0UB= date2jed([1998 12 31 0 0 0]); %(days)
%                     JD0LB= date2jed([1996 1 1 0 0 0]); %(days)
%                     

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        %Specify MCM/DSM location (X,Y,Z) constraints upper/lower bounds OR equality
        %Must have at least one nonzero element of the node's MCM constraint,
        %otherwise the node's MCM constraint is ignored
%                       node = 2;
%                     MCMUB(:,node)= [1.5 1.6 0.1]; %(Rj)
%                       node = 3;
%                     MCMLB(:,node)= [1.4 1.5 -0.1]; %(Rj)
        %alternatively:
%                     MCMUB(:,2)= [1.5 1.6 0.1]; %(Rj)
%                     MCMLB(:,3)= [0 1.1 -0.5]; %(Rj)


%##########################################################################
%#######################----LINEAR CONSTRAINTS----#########################
%##########################################################################
            

        %Specify Total TOF constraints upper/lower bounds OR equality
%                     tofTOTUB= 2800; %(days)
%                     tofTOTLB= 1000; %(days)
%                     tofTOTEQ= 1234.567; %(days)
%                     tofTOTEQ= 2451; %(days)

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        %Specify JD constraints upper/lower bounds OR equality
        %For best performance constrain initial JD(1) with JD0 options above 
%                       node = 1;            
%                     JDUB(1)= date2jed([2025 10 15 8 43 0]); %(days)
%                       node = 1;
%                     JDLB(1)= date2jed([2025 12 15 8 43 0]); %(days)
%                       node = 3;
%                     JDEQ(3)= date2jed([2025 11 15 8 43 0]); %(days)       
        		%alternatively:
%                     JDUB(1)= date2jed([2025 10 15 8 43 0]); %(days)
%                     JDLB(1)= date2jed([2025 12 15 8 43 0]); %(days)
%                     JDEQ(3)= date2jed([2025 11 15 8 43 0]); %(days)

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        %Specify MCM/DSM location (X,Y,Z) equality constraints  
        %Must have at least one nonzero element of the node's MCM constraint,
        %otherwise the node's MCM constraint is ignored

%                       node = 5;
%                     MCMEQ(:,node)= [1.435 1.567 -0.001278]; %(Rj)
        %alternatively:
        
%                     MCMEQ(:,5)= [1.435 1.567 -0.001278]; %(Rj)

%##########################################################################
%#####################----NONLINEAR CONSTRAINTS----########################
%##########################################################################

		%Specify MCM/DSM location radius mag. constraints upper/lower bounds OR equality
		%ONLY USE IF NECESSARY, this nonlinear constraint dramatically slows the
		%optimization process [MUST BE NONZERO]
%                       node = 2;
%                     MCMRUB(node)= 1.75; %(Rj)
%                       node = 3;
%                     MCMRLB(node)= 2.5; %(Rj)
%                       node = 5;
%                     MCMREQ(node)= norm([-0.117  1.570  0.028]); %(Rj)
		%alternatively:
%					  MCMRUB(2)= 1.75; %(Rj)
%					  MCMRLB(3)= 2.5; %(Rj)
%					  MCMREQ(5)= 2.6; %(Rj)

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		%Specify each node's DV mag. constraints upper/lower bounds OR equality
		%ONLY USE IF NECESSARY, this nonlinear constraint dramatically slows the
		%optimization process [MUST BE NONZERO]
%                       node = 1;
%                     DVUB(node)= 1.0; %(km/s)
%                       node = 2;
%                     DVLB(node)= 0.1; %(km/s)
%                       node = 4;
%                     DVEQ(node)= 0.0; %(km/s)
		%alternatively:
%					  DVUB(1)= 1.0; %(km/s)
% 					  DVLB(2)= 0.1; %(km/s)
%					  DVEQ(4)= 0.001; %(km/s)


% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		%Specify each node's vinfin mag. constraints upper/lower bounds OR equality
		%ONLY USE IF NECESSARY, this nonlinear constraint dramatically slows the
		%optimization process [MUST BE NONZERO]
%                       node = 1;
%                     vinfinUB(node)= 1.0; %(km/s)
%                       node = 2;
%                     vinfinLB(node)= 1.1; %(km/s)
%                       node = 4;
%                     vinfinEQ(node)= 2.0; %(km/s)
		%alternatively:
%					  vinfinUB(1)= 1.0; %(km/s)
%					  vinfinLB(2)= 1.1; %(km/s)
%					  vinfinEQ(4)= 2.0; %(km/s)

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		%Specify each node's vinfout mag. constraints upper/lower bounds OR equality
		%ONLY USE IF NECESSARY, this nonlinear constraint dramatically slows the
		%optimization process [MUST BE NONZERO]
%                       node = 1;
%                     vinfoutUB(node)= 1.0; %(km/s)
%                       node = 2;
%                     vinfoutLB(node)= 1.1; %(km/s)
%                       node = 4;
%                     vinfoutEQ(node)= 2.0; %(km/s)
		%alternatively:
%					  vinfoutUB(1)= 1.0; %(km/s)
%					  vinfoutLB(2)= 1.1; %(km/s)
%					  vinfoutEQ(4)= 2.0; %(km/s) 

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		%Specify each node to match vinf magnitudes (in and out) and ensure
		%the required periapse radius is greater than the minimum allowable
            % Set vinfmatch to 1 for each node intended to match 
%                       node = 2;
%                     vinfmatch(node)= 1;
%                       node = 3;
%                     vinfmatch(node)= 1;
%                       node = 4;
%                     vinfmatch(node)= 1; 
		%alternatively:
% 					  vinfmatch(2)= 1; 
% 					  vinfmatch(4)= 1;
% 					  vinfmatch(5)= 1;
% 					  vinfmatch(6)= 1;
%					  vinfmatch(3)= 1; 
%					  vinfmatch(4)= 1; 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%% ----------Call to Trajectory Configuration Master Function--------------
                
[dv,dve,JD,tof,rsoi,vinfout,vinfin,vinfpl,vinfoutmag,vinfinmag,...
    rP,vP,vPplhat,vd,va,vplhat,coet,rev,dvRA,dvDEC,rm,vout,...
 vin,rpout,rpin,vpout,vpin,rpfb,rpmag,betaout,betain,deltaout,deltain,tsoiout,...
 tsoiin,Tp,mupl,rpl,rap,dp,wp,coep,coenp,rnb,dvmag,rpfbmag,rfbc,rpfbmagr,output]=...
 trajcon(event,optimization,plotting,JD0,tof,npl,grade,rev,revig,branch,rfbc,...
 coep,coenp,rnb,vinfavmag,ilv,altp,optr,optoptions,optnode,dvopt,optscale,...
 branchopt,maxstepJD0,maxsteptof,maxstepMCM,...
 tofTOTUB,tofTOTLB,tofTOTEQ,tofUB,tofLB,tofEQ,JD0UB,JD0LB,JDUB,JDLB,JDEQ,...
 MCMUB,MCMLB,MCMEQ,MCMRUB,MCMRLB,MCMREQ,DVUB,DVLB,DVEQ,vinfinUB,vinfinLB,...
 vinfinEQ,vinfoutUB,vinfoutLB,vinfoutEQ,vinfmatch,ephtype,savefile,loadfile,...
 excelfile,writeout,pts,dt,movie,fps);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
jed2date(JD(1))
jed2date(JD(end))
