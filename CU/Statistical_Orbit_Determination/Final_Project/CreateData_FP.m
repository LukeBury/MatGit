clear
clc
% close all
addpath('../../Research/Europa_Hopper_Analysis/ProjectBin')
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Plot/Run Options; 0 = no, 1 = yes
% ------------------------------------------------------------------------
plotSystem          = 0;
plotEuropaFixed     = 0;
plotEuropaInertial  = 1;
plotEuropaAltitude  = 1;
plotJCPercentChange = 0;
plotVisibility      = 1;
plotRangeData       = 1;
% ------------------------------------------------------------------------
%%% System Parameters
% ------------------------------------------------------------------------
d2r = pi/180;
aEur = 671100; % semimajor axis (distance-normalizing value)
wEur = [0,0,2.047200349303344e-05]; % rad/s
w_n = [0, 0, wEur(3)/wEur(3)]; % angular velocity

%%% Normalizing Factors (real / xN = normalized)
wNorm = 1/wEur(3); % Angular Velocity (rad/s)
rNorm = aEur; % Position (km)
vNorm = rNorm/wNorm; % Velocity (km/s)
tNorm = wNorm;

%%% Gravitational Parameters
u1 = 126672520; % Jupiter, km^3 / s^2
u2 = 3203.413216; % Europa, km^3 / s^2
u = u2/u1;

% Body Radii
rad1_n = 71492/rNorm; % radius of Jupiter (km normalized)
rad2_n = 1560.8/rNorm; % radius of Europa (km normalized)
rEur = rad2_n.*rNorm;

%%% Rotating frame cooridinates
rB1_BCR_n = [-u, 0, 0];
rB2_BCR_n = [1-u, 0, 0];

%%% Initial Julian Date
[JD0] = JD(2025,3,5,0,0);
fprintf('***JD0 from buffinton paper\n')

%%% Earth & Jupiter Parameters
% Jupiter
rJ = rad1_n*rNorm; % km

% Earth
rE = 6378.1363; % km
wDay = [0,0,2*pi/86400]; % rad/s
fprintf('***Planetary Prms from Valado\n')
fprintf('***Assuming coplanar\n')
% ------------------------------------------------------------------------
%%% Defining Satellite State
% ------------------------------------------------------------------------
%%% Initial Particle Position
% Position of satellite wrt Europa (normalized)
% rH0_ECR_n = [rad2_n+50/rNorm, 0, 0]'
a = rad2_n*rNorm+50;
e = 0 *d2r;
i = 35 *d2r;
raan = 0 *d2r;
w = 0 *d2r;
ta = 180 *d2r;
[rH0_ECR_n, vH0_BCR_n] = OE2ECI(a, e, i, raan, w, ta, u2);
rH0_ECR_n = (rH0_ECR_n./rNorm)';
vH0_BCR_n = (vH0_BCR_n./vNorm)';

save('OrbitalInfo.mat','a','e','i','raan','w','ta')


% Position wrt Barycenter (normalized)
rH0_BCR_n = rH0_ECR_n + rB2_BCR_n;

% ------------------------------------------------------------------------
%%% Propagating State
% ------------------------------------------------------------------------
%%% Europa period (days)
pEur = 3.551181;
%%% Setting normalized time vector
day = 86400; % sec
ti = 0;
dt = 60;
tf = 1*day;
time = (ti:dt:tf)./tNorm;

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Events',@impact_FP,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (ECEF)
X0_n = [rH0_BCR_n, vH0_BCR_n]'; % km, km/s

%%% Propagating the State
[Times_n,States_BCR_n] = ode45(@integrator_FPdata,time,X0_n,options,u,rB1_BCR_n,rB2_BCR_n,rad2_n,tNorm,rNorm);
warning('Put in Jupiter J2 & J3\n')

%%% Un-Normalizing Time
Times = Times_n.*tNorm; % sec
if tf == Times(end)
    Times = [ti:dt:tf];
else
    warning('Weird collision or something\n')
end

%%% Assigning variables to state components
rH_BCR_n = States_BCR_n(:,1:3);
vH_BCR_n = States_BCR_n(:,4:6);

if plotSystem == 1
    figure; hold all
    bodySurface2(rad1_n, rB1_BCR_n, [1 .5 0],1.5)
    bodySurface2(rad2_n, rB2_BCR_n, [0 1 1], 1.5)
    plot3(States_BCR_n(:,1),States_BCR_n(:,2),States_BCR_n(:,3),'m')
    PlotBoi3('X','Y','Z',16)
    axis equal
end

% ------------------------------------------------------------------------
%%% ECEF Europa w/ Trajectory Plot
% ------------------------------------------------------------------------
if plotEuropaFixed == 1
    figure; hold all
    plot3(States_BCR_n(:,1),States_BCR_n(:,2),States_BCR_n(:,3),'m','linewidth',1.5)
    bodySurface3(rad2_n, rB2_BCR_n, [0 1 1])
    quiver3(rB2_BCR_n(1),rB2_BCR_n(2),rB2_BCR_n(3),-rB2_BCR_n(1), -rB2_BCR_n(2), -rB2_BCR_n(3),.004)
    PlotBoi3('X','Y','Z',16)
    title('ECEF')
    axis equal
end

% ------------------------------------------------------------------------
%%% Jacobi Constants
% ------------------------------------------------------------------------
[JCs_n] = JacobiConstantCalculator(u,rH_BCR_n,vH_BCR_n);

if plotJCPercentChange == 1
    figure
    plot(Times,(JCs_n-JCs_n(1))*100./(JCs_n(1)),'-o','linewidth',1.5)
    PlotBoi2('Time, sec','JC Percent Change',16)
end

% ------------------------------------------------------------------------
%%% Calculating ECI and EurCI States
% ------------------------------------------------------------------------
states_BCI = zeros(length(Times),6);
bc_JCI = zeros(length(Times),6);
europa_BCI = zeros(length(Times),6);
europa_JCI = zeros(length(Times),6); 
earth_HC = zeros(length(Times),6);
jupiter_HC = zeros(length(Times),6);
europaAltitudes = zeros(length(Times),1);
States_EurCI_truth = zeros(length(Times),6);
States_ECI = zeros(length(Times),6);
rJupiter_ECI = zeros(length(Times),3);
rEuropa_ECI = zeros(length(Times),3);
states_EurCR = zeros(length(Times),6);
states_EurCI2 = zeros(length(Times),6);

%%% Barycenter Inertial State (BCI)
for k = 1:length(Times)
    %%% Rotation Angle
    th = Times(k) * wEur(3); % rad
    
    %%% States in BCI
    states_BCI(k,1:3) = R3(States_BCR_n(k,1:3),th).*rNorm; % km
    states_BCI(k,4:6) = R3(States_BCR_n(k,4:6),th).*vNorm + cross(wEur,states_BCI(k,1:3)); % km/s
    
    %%% States in EurCR
    states_EurCR(k,1:3) = (States_BCR_n(k,1:3) - rB2_BCR_n).*rNorm; % km/s
    states_EurCR(k,4:6) = States_BCR_n(k,4:6).*vNorm; % km/s
    
    %%% EurCI 2!!!
    states_EurCI2(k,1:3) = R3(states_EurCR(k,1:3),th);
    states_EurCI2(k,4:6) = R3(states_EurCR(k,4:6),th) + cross(wEur,states_EurCI2(k,1:3));
    
    %%% Europa in BCI
    europa_BCI(k,1:3) = R3([1-u; 0; 0],th).*rNorm; % km
    europa_BCI(k,4:6) = cross(wEur,europa_BCI(k,1:3)); % km/s
    
    %%% Barycenter in Jupiter-Centered Inertial (JCI)
    bc_JCI(k,1:3) = R3([u; 0; 0],th).*rNorm; % km
    bc_JCI(k,4:6) = cross(wEur,bc_JCI(k,1:3)); % km/s
    
    %%% Europa in JCI
    europa_JCI(k,1:3) = bc_JCI(k,1:3) + europa_BCI(k,1:3); % km
    europa_JCI(k,4:6) = cross(wEur,europa_JCI(k,1:3)); % km/s
end

%%% Creating Earth and Jupiter Heliocentric Inertial States
for k = 1:length(Times)
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,3,'EME2000'); % km, km/s, km^3/s^2
    earth_HC(k,1:6) = [r',v'];
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,5,'EME2000'); % km, km/s, km^3/s^2
    jupiter_HC(k,1:6) = [r',v'];
end

%%% Creating EurCI States and Altitudes
for k = 1:length(Times)
    States_EurCI_truth(k,1:3) = states_BCI(k,1:3) - europa_BCI(k,1:3); % km
    States_EurCI_truth(k,4:6) = states_BCI(k,4:6) - europa_BCI(k,4:6); % km/s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    europaAltitudes(k) = norm(States_EurCI_truth(k,1:3)) - rad2_n*rNorm; % km
end

%%% Creating ECI States & Jupiter/Europa Positions
for k = 1:length(Times)
    States_ECI(k,1:3) = States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3); % km
    States_ECI(k,4:6) = States_EurCI_truth(k,4:6) + europa_JCI(k,4:6) + jupiter_HC(k,4:6) - earth_HC(k,4:6); % km/s
    
    rJupiter_ECI(k,:) = jupiter_HC(k,1:3) - earth_HC(k,1:3); % km
    rEuropa_ECI(k,:) = europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3); % km
end

if plotEuropaInertial == 1
    figure; hold all
    plot3(States_EurCI_truth(:,1),States_EurCI_truth(:,2),States_EurCI_truth(:,3),'k','linewidth',1.5)
    bodySurface3(rad2_n*rNorm, [0 0 0], [0 1 1])
    PlotBoi3('X','Y','Z',16)
    title('ECI')
    axis equal
end


if plotEuropaAltitude == 1
    figure; hold all
    plot(Times,europaAltitudes)
    PlotBoi2('Times, sec','Europa Altitude, km',16)
end

% ------------------------------------------------------------------------
%%% Station Positions and Velocities
% ------------------------------------------------------------------------
%%% Initial lat/lon/alts of Stations
latLonAlt = [-35.398333, 148.981944, 0.691750;
    40.427222, -4.250556, 0.834539;
    35.247164, 243.205, 1.07114904]; % deg | deg | km

%%% Calculating Intial Positions and Velocities of Stations in ECI
% statStates [6x3xTimes] = [[6x1]-Station 1 State | [6x1]-Station 2 State | [6x1]-Station 3 State] ... 3rd dimension is time
stations_ECI_truth = zeros(6,3,length(Times));

for k = 1:3
    %%% Acquiring ECEF/ECI station position
    [r] = latlon2surfECEF(latLonAlt(k,1), latLonAlt(k,2), latLonAlt(k,3) + rE); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correct way of handling altitude?

    %%% Storing initial station positions and velocities
    stations_ECI_truth(1:3,k,1) = r'; % km
    stations_ECI_truth(4:6,k,1) = cross(wDay',r'); % km/s
end; clear r

%%% Propagating Station States in ECI
for k = 2:length(Times)
    %%% Calculating current rotation angle
    th = wDay(3)*Times(k); % rad 
    for s = 1:3
        stations_ECI_truth(1:3,s,k) = R3(stations_ECI_truth(1:3,s,1),th); % km
        stations_ECI_truth(4:6,s,k) = R3(stations_ECI_truth(4:6,s,1),th); % km/s
    end
end

% ------------------------------------------------------------------------
%%% Checking Visibility, Creating Data, and Adding Noise
% ------------------------------------------------------------------------
%%% Setting Noise
% pError = (5e-3)*2; % km
% dpError = 5e-7; % km/s
pError = 0; % km
dpError = 0; % km/s
warning('Find better noise values\n')

%%% Setting Station Elevation Mask
elevMask = 10*d2r; % rad

%%% Initializing Data Vector
data = zeros(1,7);

dk       = 1;       % Overall measurement number
shadowAngles_J = zeros(length(Times),1);
shadowAngles_E = zeros(length(Times),1);
satAngles_J = zeros(length(Times),1);
satAngles_E = zeros(length(Times),1);
behindJ = NaN(length(Times),1);
behindE = NaN(length(Times),1);
ranges_truth = NaN(length(Times),3);
rangeRates_truth = NaN(length(Times),3);
rangeData_truth = zeros(length(Times),6,4);
rangeData_truth(:,:,1) = States_EurCI_truth(:,1:6);
rangeData_truth(:,:,2) = europa_JCI(:,1:6);
rangeData_truth(:,:,3) = jupiter_HC(:,1:6);
rangeData_truth(:,:,4) = earth_HC(:,1:6);
warning('For visibility, only looking at X-Y plane ... Ignoring Z since that''s messy\n')

for k = 1:length(Times)
    %%% Calculating current rotation angle
    th = wDay(3)*Times(k); % rad
    
    %%% Rotating current satellite states into ECEF
    rECEF = R3(States_ECI(k,1:3),-th); % (km)
    
    %%% Checking Visibility
    % Angle between center and rim of Jupiter in ECI
%     shadowAngles(k) = atan(rJ/norm(jupiterStates_ECI(k,1:3))); % rad
    shadowAngles_J(k) = atan(rJ/norm([rJupiter_ECI(k,1:2), 0])); % rad
    shadowAngles_E(k) = atan(rEur/norm(rEuropa_ECI(k,1:3))); % rad
    % Angle between center of Jupiter and satellite vector in ECI
%     satAngles(k) = atan2(norm(cross(jupiterStates_ECI(k,1:3),States_ECI(k,1:3))),dot(jupiterStates_ECI(k,1:3),States_ECI(k,1:3))); % rad
    satAngles_J(k) = atan2(norm(cross([rJupiter_ECI(k,1:2), 0],[States_ECI(k,1:2), 0])),dot([rJupiter_ECI(k,1:2), 0],[States_ECI(k,1:2), 0])); % rad
    satAngles_E(k) = atan2(norm(cross(rEuropa_ECI(k,1:3),States_ECI(k,1:3))),dot(rEuropa_ECI(k,1:3),States_ECI(k,1:3))); % rad
    if satAngles_J(k) > shadowAngles_J(k) || norm([States_ECI(k,1:2), 0]) < norm([rJupiter_ECI(k,1:2), 0]) % If not blocked by Jupiter
%         if satAngles_E(k) > shadowAngles_E(k) || norm(States_ECI(k,1:3)) < norm(rEuropa_ECI(k,1:3)) % If not blocked by Europa
            %%% Calculating elevations and storing ranges and range rates
            [az, el1, p] = ECEF2az_el_p(rECEF', latLonAlt(1,1)*d2r, latLonAlt(1,2)*d2r, latLonAlt(1,3));
            [az, el2, p] = ECEF2az_el_p(rECEF', latLonAlt(2,1)*d2r, latLonAlt(2,2)*d2r, latLonAlt(2,3));
            [az, el3, p] = ECEF2az_el_p(rECEF', latLonAlt(3,1)*d2r, latLonAlt(3,2)*d2r, latLonAlt(3,3));
            if el1 > elevMask
                data(dk,1) = Times(k); % sec
                ranges_truth(k,1) = sqrt(dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,1,k)',States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,1,k)')) + mvnrnd(zeros(1,1),pError^2,1); % km ;
                rangeRates_truth(k,1) = dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,1,k)',States_EurCI_truth(k,4:6) + europa_JCI(k,4:6) + jupiter_HC(k,4:6) - earth_HC(k,4:6) - stations_ECI_truth(4:6,1,k)')/ranges_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1); % km/s ;
                ranges_truth(k,2) = Times(k);
                rangeRates_truth(k,2) = Times(k);
                ranges_truth(k,3) = 1;
                rangeRates_truth(k,3) = 1;
                data(dk,2) = ranges_truth(k,1); % km
                data(dk,5) = rangeRates_truth(k,1); % km/s
                dk = dk + 1;

            elseif el2 > elevMask
                data(dk,1) = Times(k);
                ranges_truth(k,1) = sqrt(dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,2,k)',States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,2,k)')); % km ;
                ranges_truth(k,1) = ranges_truth(k,1) + mvnrnd(zeros(1,1),pError^2,1);
                rangeRates_truth(k,1) = dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,2,k)',States_EurCI_truth(k,4:6) + europa_JCI(k,4:6) + jupiter_HC(k,4:6) - earth_HC(k,4:6) - stations_ECI_truth(4:6,2,k)')/ranges_truth(k,1); % km/s;
                rangeRates_truth(k,1) = rangeRates_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1);
                ranges_truth(k,2) = Times(k);
                rangeRates_truth(k,2) = Times(k);
                ranges_truth(k,3) = 2;
                rangeRates_truth(k,3) = 2;
                data(dk,3) = ranges_truth(k,1); % km
                data(dk,6) = rangeRates_truth(k,1); % km/s
                dk = dk + 1;
                
            elseif el3 > elevMask
                data(dk,1) = Times(k);
                ranges_truth(k,1) = sqrt(dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,3,k)',States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,3,k)')) + mvnrnd(zeros(1,1),pError^2,1); % km ;
%                 rangeRates_truth(k,1) = dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,3,k)',States_EurCI_truth(k,4:6) + europa_JCI(k,4:6) + jupiter_HC(k,4:6) - earth_HC(k,4:6) - stations_ECI_truth(4:6,3,k)')/ranges_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1); % km/s;
                rangeRates_truth(k,1) = dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI_truth(1:3,3,k)',States_EurCI_truth(k,4:6) + europa_JCI(k,4:6) + jupiter_HC(k,4:6) - earth_HC(k,4:6) - stations_ECI_truth(4:6,3,k)')/ranges_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1); % km/s;
                ranges_truth(k,2) = Times(k);
                rangeRates_truth(k,2) = Times(k);
                ranges_truth(k,3) = 3;
                rangeRates_truth(k,3) = 3;
                data(dk,4) = ranges_truth(k,1); % km
                data(dk,7) = rangeRates_truth(k,1); % km/s
%                 data(dk,4) = norm(States_ECI(k,1:3) - statStates(1:3,3,k)') + mvnrnd(zeros(1,1),pError^2,1); % km
%                 data(dk,7) = dot(States_ECI(k,1:3) - statStates(1:3,3,k)',States_ECI(k,4:6) - statStates(4:6,3,k)')/data(dk,4);% + mvnrnd(zeros(1,1),dpError^2,1); % km/s
                dk = dk + 1;

            end
%         end
    end
end

for k = 1:length(Times)
    if satAngles_J(k) < shadowAngles_J(k) && norm([States_ECI(k,1:2), 0]) > norm([rJupiter_ECI(k,1:2), 0]) % If blocked by Jupiter
        behindJ(k) = satAngles_J(k);
    end
    if satAngles_E(k) < shadowAngles_E(k) && norm(States_ECI(k,1:3)) > norm(rEuropa_ECI(k,1:3)) % If blocked by Europa
        behindE(k) = satAngles_E(k);
    end
end

%%% Correcting/Configuring Data
firstTime = data(1,1);
data(data == 0) = NaN;
data(1,1) = firstTime;
fprintf('***Station 1 overwrites 2... 2 overwrites 3\n')

%%% Plotting Visibility 
if plotVisibility == 1
    figure; 
    subplot(2,1,1); hold all
    plot(Times./3600,shadowAngles_J,'k')
    plot(Times./3600,satAngles_J,'m')
    plot(Times./3600,behindJ,'r','linewidth',2)
    PlotBoi2('Times, hr','Visibility Angles, rads',16)
    title('Jupiter Visibility Angles in X-Y axis ONLY')
    legend('shadowAngle','satAngle','Behind Jupiter')
    
    subplot(2,1,2); hold all
    plot(Times./3600,shadowAngles_E,'k')
    plot(Times./3600,satAngles_E,'m')
    plot(Times./3600,behindE,'r','linewidth',2)
    PlotBoi2('Times, hr','Visibility Angles, rads',16)
    title('Europa Visibility Angles')
    legend('shadowAngle','satAngle','Behind Europa')
end

%%% Plotting Range Data
if plotRangeData == 1
    figure
    subplot(2,1,1); hold all
    plot(data(:,1)./3600,data(:,2),'.m')
    plot(data(:,1)./3600,data(:,3),'.r')
    plot(data(:,1)./3600,data(:,4),'.b')
    legend('Canberra','Madrid','Goldstone')
    PlotBoi2('','Range Data, km',16)
    
    subplot(2,1,2); hold all
    plot(data(:,1)./3600,data(:,5),'.r')
    plot(data(:,1)./3600,data(:,6),'.b')
    plot(data(:,1)./3600,data(:,7),'.m')
    legend('Canberra','Madrid','Goldstone')
    PlotBoi2('Time, hr','Range-Rate Data, km/s',16)
end

%%% Writing Data
dlmwrite('finalProjectData_ECI.txt',data,'precision',15)

%%% Saving Truth
save('truthTraj_EurCI.mat','States_EurCI_truth', 'ranges_truth', 'rangeRates_truth','rangeData_truth','stations_ECI_truth')
