clear
clc
close all
addpath('../../bin')

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
RE = 6378; % Radius of Earth (km)
uE = 398600.4415; % (km^3/s^2)
wE = [0, 0, 2*pi/86400]; % (rad/s)
Stat1 = [-35.398333, 148.981944]; % Lat, Lon (deg)
Stat2 = [40.427222, 355.749444]; % Lat, Lon (deg)
Stat3 = [35.247164, 243.205]; % Lat, Lon (deg)

d2r = pi/180; % Conversion from degrees to radians
r2d = 180/pi; % Conversion from radians to degrees
%%% LEO Satellite (ECI)
a_LEO = 7000; % (km)
e_LEO = 0.001;
i_LEO = 30*d2r; % (rad)
raan_LEO = 80*d2r; % (rad)
w_LEO = 40*d2r; % (rad)
ta_LEO = 0*d2r; % (rad)

%%% GEO Satellite (ECI)
a_GEO = 42163; % (km)
e_GEO = 1e-5;
i_GEO = 1*d2r; % (rad)
raan_GEO = 70*d2r; % (rad)
w_GEO = 0*d2r; % (rad)
ta_GEO = 180*d2r; % (rad)

%%% Hyperbolic Flyby (ECI)
r0_HV = [-7737.559071593195, -43881.87809094457, 0]; % (km)
v0_HV = [3.347424567061589, 3.828541915617483, 0]; % (km/s)

%%% Doppler Shift Variables
fT = 8.44; % (GHz)
c = 3e5; % km/s
% ------------------------------------------------------------------------
%%% Acquiring Initial States
% ------------------------------------------------------------------------
%%% LEO Satellite (ECI)
[r0_LE0, v0_LEO] = OE2ECI(a_LEO, e_LEO, i_LEO, raan_LEO, w_LEO, ta_LEO, uE);
X0_LEO = [r0_LE0, v0_LEO]; % (km, km/s)

%%% GEO Satellite (ECI)
[r0_GE0, v0_GEO] = OE2ECI(a_GEO, e_GEO, i_GEO, raan_GEO, w_GEO, ta_GEO, uE);
X0_GEO = [r0_GE0, v0_GEO]; % (km, km/s)

%%% Hyperbolic Flyby (ECI)
X0_HF = [r0_HV v0_HV]; % (km, km/s)

% ------------------------------------------------------------------------
%%% Propagating States
% ------------------------------------------------------------------------
%%% Setting time period
t0 = 0; % (sec)
dt = 10; % (sec)
tf = 24*3600; % (sec)
time = t0:dt:tf; % (sec)
fprintf('***Lower dt\n')
%%% Setting integrator options
tol = 1E-13;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating LEO
[Times,States_LEO] = ode45(@KeplarianOrbitIntegrator,time,X0_LEO,options,uE);

%%% Propagating GEO
[Times,States_GEO] = ode45(@KeplarianOrbitIntegrator,time,X0_GEO,options,uE);

%%% Propagating Hyperbolic Flyby
[Times,States_HF] = ode45(@KeplarianOrbitIntegrator,time,X0_HF,options,uE);

% ------------------------------------------------------------------------
%%% Determining Station Positions and Velocities
% ------------------------------------------------------------------------
%%% Finding station positions in ECEF
[rStat1_ECEF] = latlon2surfECEF(Stat1(1), Stat1(2), RE);
[rStat2_ECEF] = latlon2surfECEF(Stat2(1), Stat2(2), RE);
[rStat3_ECEF] = latlon2surfECEF(Stat3(1), Stat3(2), RE);

%%% Finding initial velocities of stations
vStat1_ECI_0 = cross(wE,rStat1_ECEF);
vStat2_ECI_0 = cross(wE,rStat2_ECEF);
vStat3_ECI_0 = cross(wE,rStat3_ECEF);

%%% Initializing ECI station state vectors
rStat1_ECI = zeros(length(Times),3); % (km)
rStat2_ECI = zeros(length(Times),3); % (km)
rStat3_ECI = zeros(length(Times),3); % (km)

vStat1_ECI = zeros(length(Times),3); % (km)
vStat2_ECI = zeros(length(Times),3); % (km)
vStat3_ECI = zeros(length(Times),3); % (km)

%%% Rotating station positions into ECI
for k = 1:length(Times)
    %%% Calculating current rotation angle
    th = wE(3)*Times(k); % rad
    
    %%% Rotating Positions
    rStat1_ECI(k,:) = R3(rStat1_ECEF,th);
    rStat2_ECI(k,:) = R3(rStat2_ECEF,th);
    rStat3_ECI(k,:) = R3(rStat3_ECEF,th);
    
    %%% Rotating Velocities
    vStat1_ECI(k,:) = R3(vStat1_ECI_0,th);
    vStat2_ECI(k,:) = R3(vStat2_ECI_0,th);
    vStat3_ECI(k,:) = R3(vStat3_ECI_0,th);    
end

% ------------------------------------------------------------------------
%%% Calculating Elevations, Ranges, and Range Rates
% ------------------------------------------------------------------------
%%% Initializing range vectors
range1_LEO = zeros(length(Times),1);
range2_LEO = zeros(length(Times),1);
range3_LEO = zeros(length(Times),1);

range1_GEO = zeros(length(Times),1);
range2_GEO = zeros(length(Times),1);
range3_GEO = zeros(length(Times),1);

range1_HF = zeros(length(Times),1);
range2_HF = zeros(length(Times),1);
range3_HF = zeros(length(Times),1);

%%% Initializing range-rate vectors
rangeRate1_LEO = zeros(length(Times),1);
rangeRate2_LEO = zeros(length(Times),1);
rangeRate3_LEO = zeros(length(Times),1);

rangeRate1_GEO = zeros(length(Times),1);
rangeRate2_GEO = zeros(length(Times),1);
rangeRate3_GEO = zeros(length(Times),1);

rangeRate1_HF = zeros(length(Times),1);
rangeRate2_HF = zeros(length(Times),1);
rangeRate3_HF = zeros(length(Times),1);

%%% Initializing Satellite ECEF Vectors
rLEO_ECEFs = zeros(length(Times),3);
rGEO_ECEFs = zeros(length(Times),3);
rHF_ECEFs = zeros(length(Times),3);

for k = 1:length(Times)
    %%% Calculating current rotation angle
    th = wE(3)*Times(k); % rad
    
    %%% Rotating current satellite states into ECEF
    rLEO_ECEF = R3(States_LEO(k,1:3),-th); % (km)
    rGEO_ECEF = R3(States_GEO(k,1:3),-th); % (km)
    rHF_ECEF  = R3(States_HF(k,1:3),-th); % (km)
    
    rLEO_ECEFs(k,:) = rLEO_ECEF;
    rGEO_ECEFs(k,:) = rGEO_ECEF;
    rHF_ECEFs(k,:) = rHF_ECEF;

    %%% Calculating elevations and storing ranges and range rates
    % LEO
    [az, el1, p] = ECEF2az_el_p(rLEO_ECEF', Stat1(1)*d2r, Stat1(2)*d2r, 0);
    [az, el2, p] = ECEF2az_el_p(rLEO_ECEF', Stat2(1)*d2r, Stat2(2)*d2r, 0);
    [az, el3, p] = ECEF2az_el_p(rLEO_ECEF', Stat3(1)*d2r, Stat3(2)*d2r, 0);
    if el1*r2d > 10
        range1_LEO(k) = norm(States_LEO(k,1:3)-rStat1_ECI(k,:)); % (km)
        rangeRate1_LEO(k) = dot(States_LEO(k,1:3)-rStat1_ECI(k,:), States_LEO(k,4:6)-vStat1_ECI(k,:))/range1_LEO(k); % (km/s)
    end
    if el2*r2d > 10
        range2_LEO(k) = norm(States_LEO(k,1:3)-rStat2_ECI(k,:)); % (km)
        rangeRate2_LEO(k) = dot(States_LEO(k,1:3)-rStat2_ECI(k,:), States_LEO(k,4:6)-vStat2_ECI(k,:))/range2_LEO(k); % (km/s)
    end
    if el3*r2d > 10
        range3_LEO(k) = norm(States_LEO(k,1:3)-rStat3_ECI(k,:)); % (km)
        rangeRate3_LEO(k) = dot(States_LEO(k,1:3)-rStat3_ECI(k,:), States_LEO(k,4:6)-vStat3_ECI(k,:))/range3_LEO(k); % (km/s)
    end
    
    % GEO
    [az, el1, p] = ECEF2az_el_p(rGEO_ECEF', Stat1(1)*d2r, Stat1(2)*d2r, 0);
    [az, el2, p] = ECEF2az_el_p(rGEO_ECEF', Stat2(1)*d2r, Stat2(2)*d2r, 0);
    [az, el3, p] = ECEF2az_el_p(rGEO_ECEF', Stat3(1)*d2r, Stat3(2)*d2r, 0);
    if el1*r2d > 10
        range1_GEO(k) = norm(States_GEO(k,1:3)-rStat1_ECI(k,:)); % (km)
        rangeRate1_GEO(k) = dot(States_GEO(k,1:3)-rStat1_ECI(k,:), States_GEO(k,4:6)-vStat1_ECI(k,:))/range1_GEO(k); % (km/s)
    end
    if el2*r2d > 10
        range2_GEO(k) = norm(States_GEO(k,1:3)-rStat2_ECI(k,:)); % (km)
        rangeRate2_GEO(k) = dot(States_GEO(k,1:3)-rStat2_ECI(k,:), States_GEO(k,4:6)-vStat2_ECI(k,:))/range2_GEO(k); % (km/s)
    end
    if el3*r2d > 10
        range3_GEO(k) = norm(States_GEO(k,1:3)-rStat3_ECI(k,:)); % (km)
        rangeRate3_GEO(k) = dot(States_GEO(k,1:3)-rStat3_ECI(k,:), States_GEO(k,4:6)-vStat3_ECI(k,:))/range3_GEO(k); % (km/s)
    end
    
    % HF
    [az, el1, p] = ECEF2az_el_p(rHF_ECEF', Stat1(1)*d2r, Stat1(2)*d2r, 0);
    [az, el2, p] = ECEF2az_el_p(rHF_ECEF', Stat2(1)*d2r, Stat2(2)*d2r, 0);
    [az, el3, p] = ECEF2az_el_p(rHF_ECEF', Stat3(1)*d2r, Stat3(2)*d2r, 0);
    if el1*r2d > 10
        range1_HF(k) = norm(States_HF(k,1:3)-rStat1_ECI(k,:)); % (km)
        rangeRate1_HF(k) = dot(States_HF(k,1:3)-rStat1_ECI(k,:), States_HF(k,4:6)-vStat1_ECI(k,:))/range1_HF(k); % (km/s)
    end
    if el2*r2d > 10
        range2_HF(k) = norm(States_HF(k,1:3)-rStat2_ECI(k,:)); % (km)
        rangeRate2_HF(k) = dot(States_HF(k,1:3)-rStat2_ECI(k,:), States_HF(k,4:6)-vStat2_ECI(k,:))/range2_HF(k); % (km/s)
    end
    if el3*r2d > 10
        range3_HF(k) = norm(States_HF(k,1:3)-rStat3_ECI(k,:)); % (km)
        rangeRate3_HF(k) = dot(States_HF(k,1:3)-rStat3_ECI(k,:), States_HF(k,4:6)-vStat3_ECI(k,:))/range3_HF(k); % (km/s)
    end
    
end

%%% Removing 'zeros' from vectors
range1_LEO(range1_LEO == 0) = NaN;
range2_LEO(range2_LEO == 0) = NaN;
range3_LEO(range3_LEO == 0) = NaN;
rangeRate1_LEO(rangeRate1_LEO == 0) = NaN;
rangeRate2_LEO(rangeRate2_LEO == 0) = NaN;
rangeRate3_LEO(rangeRate3_LEO == 0) = NaN;

range1_GEO(range1_GEO == 0) = NaN;
range2_GEO(range2_GEO == 0) = NaN;
range3_GEO(range3_GEO == 0) = NaN;
rangeRate1_GEO(rangeRate1_GEO == 0) = NaN;
rangeRate2_GEO(rangeRate2_GEO == 0) = NaN;
rangeRate3_GEO(rangeRate3_GEO == 0) = NaN;

range1_HF(range1_HF == 0) = NaN;
range2_HF(range2_HF == 0) = NaN;
range3_HF(range3_HF == 0) = NaN;
rangeRate1_HF(rangeRate1_HF == 0) = NaN;
rangeRate2_HF(rangeRate2_HF == 0) = NaN;
rangeRate3_HF(rangeRate3_HF == 0) = NaN;

% ------------------------------------------------------------------------
%%% Doppler Shift
% ------------------------------------------------------------------------
fDop1 = zeros(length(Times),1);
fDop2 = zeros(length(Times),1);
fDop3= zeros(length(Times),1);

for k = 1:length(Times)
    fDop1(k) = (-2*rangeRate1_HF(k)*fT/c)*1e9; % Hz
    fDop2(k) = (-2*rangeRate2_HF(k)*fT/c)*1e9; % Hz
    fDop3(k) = (-2*rangeRate3_HF(k)*fT/c)*1e9; % Hz
end


% ------------------------------------------------------------------------
%%% Adding Noise
% ------------------------------------------------------------------------
rangeRate1_HF_n = zeros(size(range1_HF));
rangeRate2_HF_n = zeros(size(range2_HF));
rangeRate3_HF_n = zeros(size(range3_HF));

for k = 1:length(Times)
    rangeRate1_HF_n(k) = rangeRate1_HF(k) + mvnrnd(zeros(1,1),0.0000005,1);
    rangeRate2_HF_n(k) = rangeRate2_HF(k) + mvnrnd(zeros(1,1),0.0000005,1);
    rangeRate3_HF_n(k) = rangeRate3_HF(k) + mvnrnd(zeros(1,1),0.0000005,1);
end
% ------------------------------------------------------------------------
%%% Plotting Ranges and Range Rates
% ------------------------------------------------------------------------
%%% LEO Ranges
figure
subplot(1,2,1); hold all
plot(Times,range1_LEO,'r','linewidth',2)
plot(Times,range2_LEO,'b','linewidth',2)
plot(Times,range3_LEO,'k','linewidth',2)
PlotBoi2('Time, sec','Range, km',16)
xlim([0 Times(end)]);

%%% LEO Range Rates
subplot(1,2,2); hold all
plot(Times,rangeRate1_LEO,'r','linewidth',2)
plot(Times,rangeRate2_LEO,'b','linewidth',2)
plot(Times,rangeRate3_LEO,'k','linewidth',2)
legend('Stat1','Stat2','Stat3')
PlotBoi2('Time, sec','Range Rate, km/s',16)
xlim([0 Times(end)]);

%%% GEO Ranges
figure
subplot(1,2,1); hold all
plot(Times,range1_GEO,'r','linewidth',2)
plot(Times,range2_GEO,'b','linewidth',2)
plot(Times,range3_GEO,'k','linewidth',2)
PlotBoi2('Time, sec','Range, km',16)
xlim([0 Times(end)]);

%%% GEO Range Rates
subplot(1,2,2); hold all
plot(Times,rangeRate1_GEO,'r','linewidth',2)
plot(Times,rangeRate2_GEO,'b','linewidth',2)
plot(Times,rangeRate3_GEO,'k','linewidth',2)
legend('Stat1','Stat2','Stat3')
PlotBoi2('Time, sec','Range Rate, km/s',16)
xlim([0 Times(end)]);

%%% HF Ranges
figure
subplot(1,2,1); hold all
plot(Times,range1_HF,'r','linewidth',2)
plot(Times,range2_HF,'b','linewidth',2)
plot(Times,range3_HF,'k','linewidth',2)
PlotBoi2('Time, sec','Range, km',16)
xlim([0 Times(end)]);

%%% HF Noisey Range Rates
subplot(1,2,2); hold all
plot(Times,rangeRate1_HF_n,'.m','linewidth',2)
plot(Times,rangeRate2_HF_n,'.m','linewidth',2)
plot(Times,rangeRate3_HF_n,'.m','linewidth',2)
xlim([0 Times(end)]);

%%% HF Range Rates
p1 = plot(Times,rangeRate1_HF,'r','linewidth',2);
p2 = plot(Times,rangeRate2_HF,'b','linewidth',2);
p3 =plot(Times,rangeRate3_HF,'k','linewidth',2);
legend([p1 p2 p3],'Stat1','Stat2','Stat3')
PlotBoi2('Time, sec','Range Rate, km/s',16)
xlim([0 Times(end)]);

%%% HF Doppler Shift
figure; hold all
plot(Times,fDop1,'r','linewidth',2);
plot(Times,fDop2,'b','linewidth',2);
plot(Times,fDop3,'k','linewidth',2);
legend('Stat1','Stat2','Stat3')
PlotBoi2('Time, sec','Doppler Shift, Hz',16)
xlim([0 Times(end)]);

%%% HF Range Rates & Noise
figure; hold all
plot(Times,rangeRate1_HF_n,'.m','linewidth',2)
plot(Times,rangeRate2_HF_n,'.m','linewidth',2)
p4 = plot(Times,rangeRate3_HF_n,'.m','linewidth',2);
p1 = plot(Times,rangeRate1_HF,'r','linewidth',2);
p2 = plot(Times,rangeRate2_HF,'b','linewidth',2);
p3 = plot(Times,rangeRate3_HF,'k','linewidth',2);
legend([p1 p2 p3 p4],'Stat1','Stat2','Stat3', 'Noisy Measurements')
PlotBoi2('Time, sec','Range, km',16)
xlim([0 Times(end)]);


% %%% Visual of simulation
% figure; hold all
% [xE,yE,zE] = sphere;
% surf(xE*RE,yE*RE,zE*RE)
% colormap([0 1 1])
% view(3);
% 
% plot3(rStat1_ECEF(1),rStat1_ECEF(2),rStat1_ECEF(3),'r.','markersize',25,'linewidth',2)
% plot3(rStat2_ECEF(1),rStat2_ECEF(2),rStat2_ECEF(3),'b.','markersize',25,'linewidth',2)
% plot3(rStat3_ECEF(1),rStat3_ECEF(2),rStat3_ECEF(3),'k.','markersize',25,'linewidth',2)
% 
% plot3(rLEO_ECEFs(:,1),rLEO_ECEFs(:,2),rLEO_ECEFs(:,3))
% plot3(rGEO_ECEFs(:,1),rGEO_ECEFs(:,2),rGEO_ECEFs(:,3))
% plot3(rHF_ECEFs(:,1),rHF_ECEFs(:,2),rHF_ECEFs(:,3))
% axis equal

