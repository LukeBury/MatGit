clear
clc
close all
addpath('../../bin')
tic

% ------------------------------------------------------------------------
%%% Plot/Run Options
% ------------------------------------------------------------------------
%%% Plot Range & Range-Rate results (P1)?
rangeRangeRatePlots = 0; % 0 = no,   1 = yes

%%% Debugging Filters
devOn             = 1; % Turn on deviation from ref trajectory
J3On              = 0; % Turn on J3 for truth trajectory
measNoise         = 1; % Turn on measurement noise
waitFirstMeas     = 1; % No filters until first measuremement 
processNoiseOn    = 0; % Turn on process noise
processNoiseLimit = 0; % # of no-measurements before turning off process noise
Q_RIC_On          = 0; % Form Q_RIC then Q (Q_ECI)?
dragOn            = 0; % Include drag in truth dynamics

%%% CKF Options
runCKF   = 1;
CKFPlots = 1;

%%% EKF Options
runEKF   = 1;
EKFPlots = 1;
ekfLimit = 0; % Measurements before switching to EKF

%%% RMS Options
rmsPlotsOn = 0;

%%% Setting time period
t0 = 0; % (sec)
dt = 10; % (sec)
tf = 24*3600; % (sec)

%%% Other
fs = 14; % Font size

% ------------------------------------------------------------------------
%%% Determining Equations of Motion
% ------------------------------------------------------------------------
syms x y z dx dy dz u RE J2 J3

r = sqrt(x^2 + y^2 + z^2);
U_J2J3 = -(3*u*RE*RE*J2*z*z)/(2*(r^5))...
    + (u*RE*RE*J2)/(2*(r^3))...
    - (5*u*RE*RE*RE*J3*z*z*z)/(2*(r^7))...
    + (3*u*RE*RE*RE*J3*z)/(2*(r^5));

Utot_uJ2J3 = -u/r + U_J2J3; % Total potential

EQM = [dx; dy; dz; diff(Utot_uJ2J3, x); diff(Utot_uJ2J3, y); diff(Utot_uJ2J3,z)];
% -EQM(4)
% -EQM(5)
% -EQM(6)
state = [x; y; z; dx; dy; dz];
Asym = jacobian(EQM, state);

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
RE = 6378; % Radius of Earth (km)
uE = 398600.4415; % (km^3/s^2)
wE = [0, 0, 2*pi/86400]; % (rad/s)
Stat1 = [-35.398333, 148.981944]; % Lat, Lon (deg)
Stat2 = [40.427222, 355.749444]; % Lat, Lon (deg)
Stat3 = [35.247164, 243.205]; % Lat, Lon (deg)
J2 = 0.0010826269;
if J3On == 1
    J3 = -2.52e-6;
else
    J3 = 0;
end

d2r = pi/180; % Conversion from degrees to radians
r2d = 180/pi; % Conversion from radians to degrees
%%% LEO Satellite (ECI)
a_Sat = 7000; % (km)
e_Sat = 0.001;
i_LEO = 30*d2r; % (rad)
raan_LEO = 80*d2r; % (rad)
w_LEO = 40*d2r; % (rad)
ta_LEO = 0*d2r; % (rad)

% ------------------------------------------------------------------------
%%% Acquiring Initial States
% ------------------------------------------------------------------------
%%% LEO Satellite (ECI)
[r0_Sat, v0_Sat] = OE2ECI(a_Sat, e_Sat, i_LEO, raan_LEO, w_LEO, ta_LEO, uE);
X0_Sat = [r0_Sat; v0_Sat]; % (km, km/s)

% ------------------------------------------------------------------------
%%% Propagating States
% ------------------------------------------------------------------------
%%% Setting time vector
time = t0:dt:tf; % (sec)

%%% Setting integrator options
tol = 1E-13;
options = odeset('RelTol',tol,'AbsTol',tol);

stm0 = eye(6);
IC_truth = [X0_Sat; reshape(stm0,36,1)];

%%% Building part of drag
if dragOn == 1
    p = 1e-12; % kg/m^3
    Cd = 2.2;
    A = 1; % m^2
    m = 200; % kg
    v = norm(v0_Sat); % km/s
    ad_RIC = [0; 500*p*v*v*Cd*A/m; 0];
else
    ad_RIC =0;
end

%%% Propagating the State (STM)
[Times,truthStateSTM] = ode45(@Hw3_STMInt,time,IC_truth,options,uE,RE,J2,J3,dragOn,ad_RIC);

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

%%% Initializing range-rate vectors
rangeRate1_LEO = zeros(length(Times),1);
rangeRate2_LEO = zeros(length(Times),1);
rangeRate3_LEO = zeros(length(Times),1);

%%% Initializing Satellite ECEF Vectors
rLEO_ECEFs = zeros(length(Times),3);

for k = 1:length(Times)
    %%% Calculating current rotation angle
    th = wE(3)*Times(k); % rad
    
    %%% Rotating current satellite states into ECEF
    rLEO_ECEF = R3(truthStateSTM(k,1:3),-th); % (km)
    
    rLEO_ECEFs(k,:) = rLEO_ECEF;

    %%% Calculating elevations and storing ranges and range rates
    % LEO
    [az, el1, p] = ECEF2az_el_p(rLEO_ECEF', Stat1(1)*d2r, Stat1(2)*d2r, 0);
    [az, el2, p] = ECEF2az_el_p(rLEO_ECEF', Stat2(1)*d2r, Stat2(2)*d2r, 0);
    [az, el3, p] = ECEF2az_el_p(rLEO_ECEF', Stat3(1)*d2r, Stat3(2)*d2r, 0);
    if el1*r2d > 10
        range1_LEO(k) = norm(truthStateSTM(k,1:3)-rStat1_ECI(k,:)); % (km)
        rangeRate1_LEO(k) = dot(truthStateSTM(k,1:3)-rStat1_ECI(k,:), truthStateSTM(k,4:6)-vStat1_ECI(k,:))/range1_LEO(k); % (km/s)
    end
    if el2*r2d > 10
        range2_LEO(k) = norm(truthStateSTM(k,1:3)-rStat2_ECI(k,:)); % (km)
        rangeRate2_LEO(k) = dot(truthStateSTM(k,1:3)-rStat2_ECI(k,:), truthStateSTM(k,4:6)-vStat2_ECI(k,:))/range2_LEO(k); % (km/s)
    end
    if el3*r2d > 10
        range3_LEO(k) = norm(truthStateSTM(k,1:3)-rStat3_ECI(k,:)); % (km)
        rangeRate3_LEO(k) = dot(truthStateSTM(k,1:3)-rStat3_ECI(k,:), truthStateSTM(k,4:6)-vStat3_ECI(k,:))/range3_LEO(k); % (km/s)
    end
end

%%% Removing 'zeros' from vectors
range1_LEO(range1_LEO == 0) = NaN;
range2_LEO(range2_LEO == 0) = NaN;
range3_LEO(range3_LEO == 0) = NaN;
rangeRate1_LEO(rangeRate1_LEO == 0) = NaN;
rangeRate2_LEO(rangeRate2_LEO == 0) = NaN;
rangeRate3_LEO(rangeRate3_LEO == 0) = NaN;

% ------------------------------------------------------------------------
%%% Adding Noise
% ------------------------------------------------------------------------
range1_LEO_n = zeros(size(range1_LEO));
range2_LEO_n = zeros(size(range2_LEO));
range3_LEO_n = zeros(size(range3_LEO));

rangeRate1_LEO_n = zeros(size(range1_LEO));
rangeRate2_LEO_n = zeros(size(range2_LEO));
rangeRate3_LEO_n = zeros(size(range3_LEO));

if measNoise == 1
    pError = 0.001;
    dpError = 0.000001;
else
    pError = 0;
    dpError = 0;
end

for k = 1:length(Times)
    range1_LEO_n(k) = range1_LEO(k) + mvnrnd(zeros(1,1),pError^2,1);
    range2_LEO_n(k) = range2_LEO(k) + mvnrnd(zeros(1,1),pError^2,1);
    range3_LEO_n(k) = range3_LEO(k) + mvnrnd(zeros(1,1),pError^2,1);
    
    rangeRate1_LEO_n(k) = rangeRate1_LEO(k) + mvnrnd(zeros(1,1),dpError^2,1);
    rangeRate2_LEO_n(k) = rangeRate2_LEO(k) + mvnrnd(zeros(1,1),dpError^2,1);
    rangeRate3_LEO_n(k) = rangeRate3_LEO(k) + mvnrnd(zeros(1,1),dpError^2,1);
end

% ------------------------------------------------------------------------
%%% Plotting Ranges and Range Rates
% ------------------------------------------------------------------------
if rangeRangeRatePlots == 1
    % Note: Remove deviation and J3 to get them to match
    %%% LEO Ranges w/ Noise
    figure
    subplot(1,2,1); hold all
    plot(Times,range1_LEO_n,'cx','markersize',10,'linewidth',2)
    plot(Times,range2_LEO_n,'cx','markersize',10,'linewidth',2)
    plot(Times,range3_LEO_n,'cx','markersize',10,'linewidth',2)
    plot(Times,range1_LEO,'ro','linewidth',2)
    plot(Times,range2_LEO,'bo','linewidth',2)
    plot(Times,range3_LEO,'ko','linewidth',2)
    PlotBoi2('Time, sec','Range, km',16)
    xlim([0 Times(end)]);

    %%% LEO Range Rates w/ Noise
    subplot(1,2,2); hold all
    plot(Times,rangeRate1_LEO_n,'cx','markersize',10,'linewidth',2)
    plot(Times,rangeRate2_LEO_n,'cx','markersize',10,'linewidth',2)
    p4 = plot(Times,rangeRate3_LEO_n,'cx','markersize',10,'linewidth',2);
    p1 = plot(Times,rangeRate1_LEO,'ro','linewidth',2);
    p2 = plot(Times,rangeRate2_LEO,'bo','linewidth',2);
    p3 = plot(Times,rangeRate3_LEO,'ko','linewidth',2);
    legend([p1 p2 p3 p4], 'Stat1','Stat2','Stat3', 'Noise')
    PlotBoi2('Time, sec','Range Rate, km/s',16)
    xlim([0 Times(end)]);
end

% ------------------------------------------------------------------------
%%% Time and Measurement Corrections
% ------------------------------------------------------------------------
%%% Reassigning "NaN" to "0"
range1_LEO_n(isnan(range1_LEO_n)==1) = 0;
range2_LEO_n(isnan(range2_LEO_n)==1) = 0;
range3_LEO_n(isnan(range3_LEO_n)==1) = 0;
rangeRate1_LEO_n(isnan(rangeRate1_LEO_n)==1) = 0;
rangeRate2_LEO_n(isnan(rangeRate2_LEO_n)==1) = 0;
rangeRate3_LEO_n(isnan(rangeRate3_LEO_n)==1) = 0;

%%% Combining Measurements (Assumes no more than one measurement per time)
noisyRangeMeas = range1_LEO_n + range2_LEO_n + range3_LEO_n;
noisyRangeRateMeas = rangeRate1_LEO_n + rangeRate2_LEO_n + rangeRate3_LEO_n;

%%% Generating vector of measurement time indices
measTimes = [];
arcTimes = zeros(13,2);
counter = 0;
for k = 1:length(Times)
   if noisyRangeMeas(k) ~= 0
       measTimes = [measTimes; Times(k)];
   end
   
   if k > 1 && noisyRangeMeas(k) ~=0 && noisyRangeMeas(k-1) == 0
       counter = counter + 1; 
       arcTimes(counter,1) = Times(k);
   end
   if k > 1 && noisyRangeMeas(k) == 0 && noisyRangeMeas(k-1) ~= 0
       arcTimes(counter,2) = Times(k-1);
   end
end

%%% If filters should wait for first measurement
if waitFirstMeas == 1
    mTimeIndex = find(Times==measTimes(1));
    noisyRangeMeas = noisyRangeMeas(mTimeIndex:end);
    noisyRangeRateMeas = noisyRangeRateMeas(mTimeIndex:end);
    
    rStat1_ECI = rStat1_ECI(mTimeIndex:end,:);
    rStat2_ECI = rStat2_ECI(mTimeIndex:end,:);
    rStat3_ECI = rStat3_ECI(mTimeIndex:end,:);
    vStat1_ECI = vStat1_ECI(mTimeIndex:end,:);
    vStat2_ECI = vStat2_ECI(mTimeIndex:end,:);
    vStat3_ECI = vStat3_ECI(mTimeIndex:end,:);
    
    range1_LEO_n = range1_LEO_n(mTimeIndex:end);
    range2_LEO_n = range2_LEO_n(mTimeIndex:end);
    range3_LEO_n = range3_LEO_n(mTimeIndex:end);

    truthStateSTM = truthStateSTM(mTimeIndex:end,:);
    
    Times = Times(mTimeIndex:end);
else
    mTimeIndex = 1;
end

% ------------------------------------------------------------------------
%%% Processing the State Transition Matrix
% ------------------------------------------------------------------------
%%% Reshaping ref states and STMs
States_truth = zeros(length(Times),6);
STMs_truth = zeros(6,6,length(Times));
for k = 1:length(Times)
    %%% Reference States
    States_truth(k,:) = truthStateSTM(k,1:6);
    
    %%% STMs
    STMs_truth(:,:,k) = reshape(truthStateSTM(k,7:end),6,6);
end

% ------------------------------------------------------------------------
%%% Preparing for Filters
% ------------------------------------------------------------------------
%%% Setting up Initial Conditions
if devOn == 1
    dev = [.003; .001; .002; .000003; .000001;.000002];
else
    dev = zeros(6,1);
end

XBar0 = States_truth(1,:)' + dev;
stm0 = eye(6);
IC_ref = [XBar0; reshape(stm0,36,1)];

%%% Measurement noise / uncertainty
R = [pError^2 0; 0 dpError^2];
% PBar0 = blkdiag(eye(3)*pError^2,eye(3)*dpError^2)*10000;
PBar0 = blkdiag(eye(3)*pError^2,eye(3)*dpError^2)*500;

pnSig = 10^-7;
%%% Process noise / uncertainty
if processNoiseOn == 1
    if Q_RIC_On == 0
        Q = [(pnSig)^2 0 0;
            0 (pnSig)^2 0;
            0 0 (pnSig)^2];
    else
        %%% Developing Q_RIC and Q
        Rbar = States_truth(1,1:3)./norm(States_truth(1,1:3));
        h = cross(States_truth(1,1:3),States_truth(1,4:6));
        Cbar = h./norm(h);
        Ibar = cross(Cbar,Rbar);
        B_RIC2ECI = [Rbar; Ibar; Cbar]';

        Q_RIC = [(1e-7)^2 0 0;
                0 (1e-7)^2 0;
                0 0 (1e-7)^2];
        Q = B_RIC2ECI'*Q_RIC*B_RIC2ECI;
    end
    gamma = dt*[eye(3)*(dt/2); eye(3)];
else
    Q = zeros(3,3);
    gamma = zeros(6,3);
end

% ========================================================================
%%% CKF (p203)
% ========================================================================
%%% Matrices to store results
refStatesSTM_CKF = zeros(length(Times),42);
xHats_CKF = zeros(length(Times),6);
P_est_CKF = zeros(6,6,length(Times));
yis_C = zeros(2,length(Times)); % Prefit residual
eis_C = zeros(2,length(Times)); % Posfit residual
est_CKF = zeros(length(Times),6); % for state estimates

P_est_CKF(:,:,1) = PBar0;
xHats_CKF(1,:) = dev;
Pi = PBar0;

refStatesSTM_CKF(1,:) = IC_ref';
est_CKF(1,:) = refStatesSTM_CKF(1,1:6) + xHats_CKF(1,:);

xhat = dev;
PBari = PBar0;

if runCKF == 1
fprintf('CKF Started: %3.1f sec\n',toc)

%%% Integrate ref trajectory and STM
[rT,refStatesSTM_CKF] = ode45(@Hw3_STMInt,Times,IC_ref,options,uE,RE,J2,0,0,0); % J3 = 0, drag off

noMeasCount = 0;
temp = 0; % For summing post fit residuals
for k = 2:length(Times)
        
    %%% reshape new STM
    refSTM_CKF = reshape(refStatesSTM_CKF(k,7:end),6,6);
    refSTM_CKF =  refSTM_CKF * inv(reshape(refStatesSTM_CKF(k-1,7:end),6,6));
        
    if noisyRangeMeas(k) ~= 0 % If there is a measurement
        if ismember(Times(k),measTimes) ~= 1
            warning('CKF Time index error\n')
        end
        noMeasCount = 0;
        %%% Renaming Current Measurements
        Yi = [noisyRangeMeas(k); noisyRangeRateMeas(k)];
        
        %%% Time Update
        xBari = refSTM_CKF*xHats_CKF(k-1,:)';
        PBari = refSTM_CKF*P_est_CKF(:,:,k-1)*refSTM_CKF';
        
        %%% Generating Htilde and Reading Observation
        Htilde = zeros(2,6); % Observation-State Matrix
        
        %%% If measurement is from station 1
        if range1_LEO_n(k) ~= 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) == 0
            %%% DP/DX
            Htilde(1,1) = DPDxyzCalculator(1, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:));
            Htilde(1,2) = DPDxyzCalculator(2, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:));
            Htilde(1,3) = DPDxyzCalculator(3, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:));
            %%% DdP/DX
            Htilde(2,1) = DdPDxyzCalculator(1, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat1_ECI(k,:));
            Htilde(2,2) = DdPDxyzCalculator(2, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat1_ECI(k,:));
            Htilde(2,3) = DdPDxyzCalculator(3, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat1_ECI(k,:));
            Htilde(2,4) = DdPDvxvyvzCalculator(1, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:));
            Htilde(2,5) = DdPDvxvyvzCalculator(2, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:));
            Htilde(2,6) = DdPDvxvyvzCalculator(3, refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:));
            %%% Prefit Residual
            yi = [Yi(1) - norm(refStatesSTM_CKF(k,1:3) - rStat1_ECI(k,:));...
                Yi(2) -  rangeRateCalculator(refStatesSTM_CKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat1_ECI(k,:))];
            
        %%% If measurement is from station 2
        elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) ~= 0 && range3_LEO_n(k) == 0
            %%% DP/DX
            Htilde(1,1) = DPDxyzCalculator(1, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:));
            Htilde(1,2) = DPDxyzCalculator(2, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:));
            Htilde(1,3) = DPDxyzCalculator(3, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:));
            %%% DdP/DX
            Htilde(2,1) = DdPDxyzCalculator(1, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat2_ECI(k,:));
            Htilde(2,2) = DdPDxyzCalculator(2, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat2_ECI(k,:));
            Htilde(2,3) = DdPDxyzCalculator(3, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat2_ECI(k,:));
            Htilde(2,4) = DdPDvxvyvzCalculator(1, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:));
            Htilde(2,5) = DdPDvxvyvzCalculator(2, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:));
            Htilde(2,6) = DdPDvxvyvzCalculator(3, refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:));
            %%% Prefit Residual
            yi = [Yi(1) - norm(refStatesSTM_CKF(k,1:3) - rStat2_ECI(k,:));...
                Yi(2) -  rangeRateCalculator(refStatesSTM_CKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat2_ECI(k,:))];
        
        %%% If measurement is from station 3
        elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) ~= 0
            %%% DP/DX
            Htilde(1,1) = DPDxyzCalculator(1, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:));
            Htilde(1,2) = DPDxyzCalculator(2, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:));
            Htilde(1,3) = DPDxyzCalculator(3, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:));
            %%% DdP/DX
            Htilde(2,1) = DdPDxyzCalculator(1, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat3_ECI(k,:));
            Htilde(2,2) = DdPDxyzCalculator(2, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat3_ECI(k,:));
            Htilde(2,3) = DdPDxyzCalculator(3, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat3_ECI(k,:));
            Htilde(2,4) = DdPDvxvyvzCalculator(1, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:));
            Htilde(2,5) = DdPDvxvyvzCalculator(2, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:));
            Htilde(2,6) = DdPDvxvyvzCalculator(3, refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:));
            %%% Prefit Residual
            yi = [Yi(1) - norm(refStatesSTM_CKF(k,1:3) - rStat3_ECI(k,:));...
                Yi(2) -  rangeRateCalculator(refStatesSTM_CKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_CKF(k,4:6), vStat3_ECI(k,:))];
            
        else
            warning('Trouble parsing noisy station measurements\n')
        end
        
        %%% DP/Dvx = DP/Dvy = DP/Dvz = 0
        Htilde(1,4) = 0;
        Htilde(1,5) = 0;
        Htilde(1,6) = 0;
        
        %%% Kalman Gain Matrix
        Ki = PBari*Htilde'*inv(Htilde*PBari*Htilde' + R);
        
        %%% Measurement Update
        xhat = xBari + Ki*(yi - Htilde*xBari);
        Pi = (eye(6) - Ki*Htilde)*PBari;
        
        %%% Storing Values
        xHats_CKF(k,:) = xhat;
        if Q_RIC_On == 1
            %%% Developing Q from Q_RIC
            r_temp = refStatesSTM_CKF(k,1:3) + xHats_CKF(k,1:3);
            v_temp = refStatesSTM_CKF(k,4:6) + xHats_CKF(k,4:6);
            Rbar = r_temp./norm(r_temp);
            h_temp = cross(r_temp,v_temp);
            Cbar = h_temp./norm(h_temp);
            Ibar = cross(Cbar,Rbar);
            B_RIC2ECI = [Rbar; Ibar; Cbar]';
            Q = B_RIC2ECI'*Q_RIC*B_RIC2ECI;
        end
        P_est_CKF(:,:,k) = Pi + gamma*Q*gamma';
        yis_C(:,k) = yi;
        
        %%% Post fit residual
        ei = yis_C(:,k) - Htilde*xHats_CKF(k,:)'; % 2x1
        temp = temp + ei'*R*ei; % 1x1
        eis_C(:,k) = ei; % 2x1

    elseif noisyRangeMeas(k) == 0 % If no measurement
        noMeasCount = noMeasCount + 1;
        
        %%% Propagate P and xHat with STM
        xHats_CKF(k,:) = refSTM_CKF*xHats_CKF(k-1,:)';
        
        if noMeasCount > processNoiseLimit % No process noise for big gaps
            P_est_CKF(:,:,k) = refSTM_CKF*P_est_CKF(:,:,k-1)*refSTM_CKF';
        else
            if Q_RIC_On == 1
                %%% Developing Q from Q_RIC
                r_temp = refStatesSTM_CKF(k,1:3) + xHats_CKF(k,1:3);
                v_temp = refStatesSTM_CKF(k,4:6) + xHats_CKF(k,4:6);
                Rbar = r_temp./norm(r_temp);
                h_temp = cross(r_temp,v_temp);
                Cbar = h_temp./norm(h_temp);
                Ibar = cross(Cbar,Rbar);
                B_RIC2ECI = [Rbar; Ibar; Cbar]';
                Q = B_RIC2ECI'*Q_RIC*B_RIC2ECI;
            end
            P_est_CKF(:,:,k) = refSTM_CKF*P_est_CKF(:,:,k-1)*refSTM_CKF' + gamma*Q*gamma';
        end
    end
    est_CKF(k,:) = refStatesSTM_CKF(k,1:6) + xHats_CKF(k,:);
end
clear r_temp v_temp h_temp
%%% Computing RMS values
m = 2*length(measTimes); % Total number of measurements
RMS_CKF = (temp/m)^(1/2);

%%% Computing estimate error
error_CKF = truthStateSTM(:,1:6) - est_CKF(:,:);

%%% CKF Plots
if CKFPlots == 1
    eis_C(eis_C == 0) = NaN;
    figure
    subplot(2,1,1); hold all; xlim([Times(1) Times(end)]);
    plot(Times,eis_C(1,:),'-mx')
    plot([Times(1) Times(end)],[3*pError 3*pError],'--r')
    plot([Times(1) Times(end)],[-3*pError -3*pError],'--r')
    legend('Post fit','3\sigma Error')
    PlotBoi2('','Range Residuals, km',14)
    subplot(2,1,2); hold all; xlim([Times(1) Times(end)]);
    plot(Times,eis_C(2,:),'-mx')
    plot([Times(1) Times(end)],[3*dpError 3*dpError],'--r')
    plot([Times(1) Times(end)],[-3*dpError -3*dpError],'--r')
    PlotBoi2('Times, sec','Range Rate Residuals, km/s',14)
    eis_C(isnan(eis_C)==1) = 0;
    
    figure
    subplot(3,2,1); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_CKF(:,1)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_CKF(1,1,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_CKF(1,1,:))),'--r')
    PlotBoi2('','X Position Error, km',fs)
    subplot(3,2,3); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_CKF(:,2)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_CKF(2,2,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_CKF(2,2,:))),'--r')
    PlotBoi2('','Y Position Error, km',fs)
    subplot(3,2,5); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_CKF(:,3)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_CKF(3,3,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_CKF(3,3,:))),'--r')
    PlotBoi2('Time, sec','Z Position Error, km',fs)
    subplot(3,2,2); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_CKF(:,4)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_CKF(4,4,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_CKF(4,4,:))),'--r')
    legend('Truth','3\sigma Covariance')
    PlotBoi2('','X Velocity Error, km/s',fs)
    subplot(3,2,4); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_CKF(:,5)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_CKF(5,5,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_CKF(5,5,:))),'--r')
    PlotBoi2('','Y Velocity Error, km/s',fs)
    subplot(3,2,6); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_CKF(:,6)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_CKF(6,6,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_CKF(6,6,:))),'--r')
    PlotBoi2('Time, sec','Z Velocity Error, km/s',fs)
    
end

fprintf('CKF Finished: %3.1f\n',toc)
end

% ========================================================================
%%% EKF (p212)
% ========================================================================
%%% Matrices to store results
refStatesSTM_EKF = zeros(length(Times),42);
xHats_EKF = zeros(length(Times),6);
P_est_EKF = zeros(6,6,length(Times));
yis_E = zeros(2,length(Times)); % Prefit residual
eis_E = zeros(2,length(Times)); % Posfit residual
est_EKF = zeros(length(Times),6); % for state estimates

P_est_EKF(:,:,1) = PBar0;
xHats_EKF(1,:) = dev;
Pi = P_est_EKF(:,:,1);
refStatesSTM_EKF(1,:) = refStatesSTM_CKF(1,:);
est_EKF(1,:) = refStatesSTM_EKF(1,1:6) + xHats_EKF(1,:);

refStatesSTM_EKF(1,:) = IC_ref';

xhat = dev;
PBari = PBar0;

if runEKF == 1
fprintf('EKF Started: %3.1f sec\n',toc)
measCounter = 0;
temp = 0; % For summing post fit residuals
k = 2;
fprintf('NOTE: EKF assumes Times(1) is a measurement time\n')
secondPass = 0;
fprintf('*****Later Note: Need to add gamma*Q*gamma in time update, not measurement update******\n')
for m = 2:length(measTimes) % Assumes no more states after last measurement .. Assumes starting on measurement
    %%% Grabbing previous state & stm
    pvsStateSTM = refStatesSTM_EKF(k-1,:);
    
    %%% Resetting previous stm to identity
    pvsStateSTM(7:end) = reshape(eye(6),36,1);
    
    %%% Integrate ref trajectory and STM since last measurement
    intTime = [measTimes(m-1):dt:measTimes(m)];
    [rT,rSSTM] = ode45(@Hw3_STMInt,intTime,pvsStateSTM,options,uE,RE,J2,0,0,0); % J3 = 0, drag = 0
    if (rT(end)-rT(1)) > dt % there are times w/out measurements ...rT(2:end-1)
        for mm = 2:length(rT)-1
            if Times(k) ~= rT(mm)
                warning('Time Indexing Broken')
                return
            end

            %%% Grab current state/stm
            refStatesSTM_EKF(k,:) = rSSTM(mm,:);

            %%% reshape new STM
            refSTM_EKF = reshape(refStatesSTM_EKF(k,7:end),6,6);
            if mm > 2 % At mm = 2, the STM already maps the previous state to current
                refSTM_EKF = refSTM_EKF * inv(reshape(refStatesSTM_EKF(k-1,7:end),6,6));
            end
            
            %%% Propagate P and xHat with STM
            xHats_EKF(k,:) = refSTM_EKF*xHats_EKF(k-1,:)';
            
            if (mm-1) > processNoiseLimit % No process noise for big gaps
                P_est_EKF(:,:,k) = refSTM_EKF*P_est_EKF(:,:,k-1)*refSTM_EKF';
            else
                if Q_RIC_On == 1
                    %%% Developing Q from Q_RIC
                    r_temp = refStatesSTM_EKF(k,1:3) + xHats_EKF(k,1:3);
                    v_temp = refStatesSTM_EKF(k,4:6) + xHats_EKF(k,4:6);
                    Rbar = r_temp./norm(r_temp);
                    h_temp = cross(r_temp,v_temp);
                    Cbar = h_temp./norm(h_temp);
                    Ibar = cross(Cbar,Rbar);
                    B_RIC2ECI = [Rbar; Ibar; Cbar]';
                    Q = B_RIC2ECI'*Q_RIC*B_RIC2ECI;
                end
                P_est_EKF(:,:,k) = refSTM_EKF*P_est_EKF(:,:,k-1)*refSTM_EKF' + gamma*Q*gamma';
            end
            
            %%% Store estimate
            refStatesSTM_EKF(k,1:6) = refStatesSTM_EKF(k,1:6) + xHats_EKF(k,:);
            est_EKF(k,:) = refStatesSTM_EKF(k,1:6);

            k = k+1;
        end
        measCounter = 0;
    end

    measCounter = measCounter + 1;

    %%% Grab current state/stm
    refStatesSTM_EKF(k,:) = rSSTM(end,:);

    %%% reshape new STM
    refSTM_EKF = reshape(refStatesSTM_EKF(k,7:end),6,6);
    if Times(k) == rT(end) && rT(end)-rT(1) > dt % If looking at the first first meas of non-first arc
        refSTM_EKF = refSTM_EKF*inv(reshape(refStatesSTM_EKF(k-1,7:end),6,6));
    end

    %%% Renaming Current Measurements
    Yi = [noisyRangeMeas(k); noisyRangeRateMeas(k)];

    %%% Time Update
    if measCounter > ekfLimit % If EKF
        xBari = zeros(6,1);
    else
        xBari = refSTM_EKF*xHats_EKF(k-1,:)'; % If CKF
    end

    PBari = refSTM_EKF*P_est_EKF(:,:,k-1)*refSTM_EKF';

    %%% Generating Htilde and Reading Observation
    Htilde = zeros(2,6); % Observation-State Matrix

    %%% If measurement is from station 1
    if range1_LEO_n(k) ~= 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) == 0
        %%% DP/DX
        Htilde(1,1) = DPDxyzCalculator(1, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:));
        Htilde(1,2) = DPDxyzCalculator(2, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:));
        Htilde(1,3) = DPDxyzCalculator(3, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:));
        %%% DdP/DX
        Htilde(2,1) = DdPDxyzCalculator(1, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat1_ECI(k,:));
        Htilde(2,2) = DdPDxyzCalculator(2, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat1_ECI(k,:));
        Htilde(2,3) = DdPDxyzCalculator(3, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat1_ECI(k,:));
        Htilde(2,4) = DdPDvxvyvzCalculator(1, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:));
        Htilde(2,5) = DdPDvxvyvzCalculator(2, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:));
        Htilde(2,6) = DdPDvxvyvzCalculator(3, refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:));
        %%% Prefit Residual
        yi = [Yi(1) - norm(refStatesSTM_EKF(k,1:3) - rStat1_ECI(k,:));...
            Yi(2) -  rangeRateCalculator(refStatesSTM_EKF(k,1:3), rStat1_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat1_ECI(k,:))];

    %%% If measurement is from station 2
    elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) ~= 0 && range3_LEO_n(k) == 0
        %%% DP/DX
        Htilde(1,1) = DPDxyzCalculator(1, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:));
        Htilde(1,2) = DPDxyzCalculator(2, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:));
        Htilde(1,3) = DPDxyzCalculator(3, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:));
        %%% DdP/DX
        Htilde(2,1) = DdPDxyzCalculator(1, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat2_ECI(k,:));
        Htilde(2,2) = DdPDxyzCalculator(2, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat2_ECI(k,:));
        Htilde(2,3) = DdPDxyzCalculator(3, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat2_ECI(k,:));
        Htilde(2,4) = DdPDvxvyvzCalculator(1, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:));
        Htilde(2,5) = DdPDvxvyvzCalculator(2, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:));
        Htilde(2,6) = DdPDvxvyvzCalculator(3, refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:));
        %%% Prefit Residual
        yi = [Yi(1) - norm(refStatesSTM_EKF(k,1:3) - rStat2_ECI(k,:));...
            Yi(2) -  rangeRateCalculator(refStatesSTM_EKF(k,1:3), rStat2_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat2_ECI(k,:))];

    %%% If measurement is from station 3
    elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) ~= 0
        %%% DP/DX
        Htilde(1,1) = DPDxyzCalculator(1, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:));
        Htilde(1,2) = DPDxyzCalculator(2, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:));
        Htilde(1,3) = DPDxyzCalculator(3, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:));
        %%% DdP/DX
        Htilde(2,1) = DdPDxyzCalculator(1, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat3_ECI(k,:));
        Htilde(2,2) = DdPDxyzCalculator(2, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat3_ECI(k,:));
        Htilde(2,3) = DdPDxyzCalculator(3, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat3_ECI(k,:));
        Htilde(2,4) = DdPDvxvyvzCalculator(1, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:));
        Htilde(2,5) = DdPDvxvyvzCalculator(2, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:));
        Htilde(2,6) = DdPDvxvyvzCalculator(3, refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:));
        %%% Prefit Residual
        yi = [Yi(1) - norm(refStatesSTM_EKF(k,1:3) - rStat3_ECI(k,:));...
            Yi(2) -  rangeRateCalculator(refStatesSTM_EKF(k,1:3), rStat3_ECI(k,:), refStatesSTM_EKF(k,4:6), vStat3_ECI(k,:))];
            
    else
        warning('Trouble parsing noisy station measurements\n')
    end
        
    %%% DP/Dvx = DP/Dvy = DP/Dvz = 0
    Htilde(1,4) = 0;
    Htilde(1,5) = 0;
    Htilde(1,6) = 0;

    %%% Kalman Gain Matrix
    Ki = PBari*Htilde'*inv(Htilde*PBari*Htilde' + R);

    %%% Measurement Update
    xhat = xBari + Ki*(yi - Htilde*xBari);
    Pi = (eye(6) - Ki*Htilde)*PBari;

    %%% Storing Values
    xHats_EKF(k,:) = xhat;
    if Q_RIC_On == 1
        %%% Developing Q from Q_RIC
        r_temp = refStatesSTM_EKF(k,1:3) + xHats_EKF(k,1:3);
        v_temp = refStatesSTM_EKF(k,4:6) + xHats_EKF(k,4:6);
        Rbar = r_temp./norm(r_temp);
        h_temp = cross(r_temp,v_temp);
        Cbar = h_temp./norm(h_temp);
        Ibar = cross(Cbar,Rbar);
        B_RIC2ECI = [Rbar; Ibar; Cbar]';
        Q = B_RIC2ECI'*Q_RIC*B_RIC2ECI;
    end
    P_est_EKF(:,:,k) = Pi + gamma*Q*gamma';
    yis_E(:,k) = yi;

    %%% Post fit residual
    ei = yis_E(:,k) - Htilde*xHats_EKF(k,:)'; % 2x1
    temp = temp + ei'*R*ei; % 1x1
    eis_E(:,k) = ei; % 2x1

    if measCounter > ekfLimit % If EKF
        refStatesSTM_EKF(k,1:6) = refStatesSTM_EKF(k,1:6) + xHats_EKF(k,:);
        est_EKF(k,:) = refStatesSTM_EKF(k,1:6);
    else % If CKF
        989
        est_EKF(k,:) = refStatesSTM_EKF(k,1:6) + xHats_EKF(k,:);
    end
%     if Times(k) == 38050
%         989
%         return
%     end

k = k+1;
end
if rT(end) < Times(end) % If there are more unmeasured times, propagate them
    %%% Grabbing previous state & stm
    pvsStateSTM = refStatesSTM_EKF(k-1,:);

    %%% Resetting previous stm to identity
    pvsStateSTM(7:end) = reshape(eye(6),36,1);

    %%% Integrate ref trajectory and STM since last measurement
    intTime = [Times(k-1):dt:Times(end)];
    [rT,rSSTM] = ode45(@Hw3_STMInt,intTime,pvsStateSTM,options,uE,RE,J2,0,0,0); % J3 = 0, drag off
    
    for mm = 2:length(rT)
        %%% Grab current state/stm
        refStatesSTM_EKF(k,:) = rSSTM(mm,:);

        %%% reshape new STM
        refSTM_EKF = reshape(refStatesSTM_EKF(k,7:end),6,6);
        if mm > 2 % At mm = 2, the STM already maps the previous state to current
            refSTM_EKF = refSTM_EKF * inv(reshape(refStatesSTM_EKF(k-1,7:end),6,6));
        end
        

        %%% Propagate P and xHat with STM
        xHats_EKF(k,:) = refSTM_EKF*xHats_EKF(k-1,:)';
        
        if (mm-1) > processNoiseLimit % No process noise for big gaps
            P_est_EKF(:,:,k) = refSTM_EKF*P_est_EKF(:,:,k-1)*refSTM_EKF';
        else
            if Q_RIC_On == 1
                %%% Developing Q from Q_RIC
                r_temp = refStatesSTM_EKF(k,1:3) + xHats_EKF(k,1:3);
                v_temp = refStatesSTM_EKF(k,4:6) + xHats_EKF(k,4:6);
                Rbar = r_temp./norm(r_temp);
                h_temp = cross(r_temp,v_temp);
                Cbar = h_temp./norm(h_temp);
                Ibar = cross(Cbar,Rbar);
                B_RIC2ECI = [Rbar; Ibar; Cbar]';
                Q = B_RIC2ECI'*Q_RIC*B_RIC2ECI;
            end
            P_est_EKF(:,:,k) = refSTM_EKF*P_est_EKF(:,:,k-1)*refSTM_EKF' + gamma*Q*gamma';
        end

        %%% Store estimate
        est_EKF(k,:) = refStatesSTM_EKF(k,1:6) + xHats_EKF(k,:);
        
        k = k+1;
    end
end
%%% Computing RMS values
m = 2*length(measTimes); % Total number of measurements
RMS_EKF = (temp/m)^(1/2);

%%% Computing estimate error
error_EKF = truthStateSTM(1:end,1:6) - est_EKF(1:end,:);

%%% EKF Plots
if EKFPlots == 1
    eis_E(eis_E == 0) = NaN;
    figure
    subplot(2,1,1); hold all; xlim([Times(1) Times(end)]);
    plot(Times,eis_E(1,:),'-mx')
    plot([Times(1) Times(end)],[3*pError 3*pError],'--r')
    plot([Times(1) Times(end)],[-3*pError -3*pError],'--r')
    legend('Post fit','3\sigma Error')
    PlotBoi2('','Range Residuals, km',14)
    subplot(2,1,2); hold all; xlim([Times(1) Times(end)]);
    plot(Times,eis_E(2,:),'-mx')
    plot([Times(1) Times(end)],[3*dpError 3*dpError],'--r')
    plot([Times(1) Times(end)],[-3*dpError -3*dpError],'--r')
    PlotBoi2('Times, sec','Range Rate Residuals, km/s',14)
    eis_E(isnan(eis_E)==1) = 0;
    
    figure
    subplot(3,2,1); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_EKF(:,1)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_EKF(1,1,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_EKF(1,1,:))),'--r')
    PlotBoi2('','X Position Error, km',fs)
    subplot(3,2,3); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_EKF(:,2)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_EKF(2,2,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_EKF(2,2,:))),'--r')
    PlotBoi2('','Y Position Error, km',fs)
    subplot(3,2,5); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_EKF(:,3)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_EKF(3,3,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_EKF(3,3,:))),'--r')
    PlotBoi2('Time, sec','Z Position Error, km',fs)
    subplot(3,2,2); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_EKF(:,4)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_EKF(4,4,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_EKF(4,4,:))),'--r')
    legend('Truth','3\sigma Covariance')
    PlotBoi2('','X Velocity Error, km/s',fs)
    subplot(3,2,4); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_EKF(:,5)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_EKF(5,5,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_EKF(5,5,:))),'--r')
    PlotBoi2('','Y Velocity Error, km/s',fs)
    subplot(3,2,6); hold all; xlim([Times(1) Times(end)]);
    plot(Times, error_EKF(:,6)) %%% difference between truth and Xref + xHat
    plot(Times,3*sqrt(squeeze(P_est_EKF(6,6,:))),'--r')
    plot(Times,-3*sqrt(squeeze(P_est_EKF(6,6,:))),'--r')
    PlotBoi2('Time, sec','Z Velocity Error, km/s',fs)
    
end

fprintf('EKF Finished: %3.1f\n',toc)
end

figure
subplot(3,2,1)
plot(Times,error_CKF(:,1)-error_EKF(:,1));
subplot(3,2,3)
plot(Times,error_CKF(:,2)-error_EKF(:,2));
subplot(3,2,5)
plot(Times,error_CKF(:,3)-error_EKF(:,3));
subplot(3,2,2)
plot(Times,error_CKF(:,4)-error_EKF(:,4));
subplot(3,2,4)
plot(Times,error_CKF(:,5)-error_EKF(:,5));
subplot(3,2,6)
plot(Times,error_CKF(:,6)-error_EKF(:,6));



% ------------------------------------------------------------------------
%%% RMS work (problem 1b)
% ------------------------------------------------------------------------

sumPos_CKF = 0;
sumVel_CKF = 0;
sumPos_EKF = 0;
sumVel_EKF = 0;
for k = 1:length(Times)
    sumPos_CKF = sumPos_CKF + norm(error_CKF(k,1:3));
    sumVel_CKF = sumVel_CKF + norm(error_CKF(k,4:6));
    
    sumPos_EKF = sumPos_EKF + norm(error_EKF(k,1:3));
    sumVel_EKF = sumVel_EKF + norm(error_EKF(k,4:6));
end

RMSPos_CKF = (sumPos_CKF/m)^(1/2)
RMSVel_CKF = (sumVel_CKF/m)^(1/2)

RMSPos_EKF = (sumPos_EKF/m)^(1/2)
RMSVel_EKF = (sumVel_EKF/m)^(1/2)

%%% Sigma values process noise (10^x)
sigs = [-15:1:-2];

%%% CKF & EKF post-fit residuals as function of sigs
ckf_RMS = [4.0320E-05, 4.0359E-05, 4.0302E-05, 4.0358E-05, 4.0568E-05, 3.9435E-05, 1.7499E-05, 3.1815E-06, 6.9995E-07, 6.5738E-07, 5.4358E-07, 2.8672E-07, 2.1591E-08, 2.5493E-10];
ekf_RMS = [1.4455E-04, 1.5108E-04, 1.4971E-04, 1.5121E-04, 1.5352E-04, 1.3436E-04, 9.6325E-05, 1.7069E-04, 4.0040E-01, 1.1355E+00, 2.3288E+00, 1.3905E+00, 1.0352E+00, 6.3550E-01];

%%% 3D RMS for position and velocity
ckfPos_RMS = [1.1065, 1.1186, 1.112, 1.1054, 1.1354, 1.0926, 0.8558, 0.7635, 0.7133, 0.7469, 1.3005, 4.7374, 8.3757, 9.4959];
ckfVel_RMS = [0.0341, 0.0346, 0.0344, 0.0341, 0.0355, 0.0347, 0.0278, 0.025, 0.0236, 0.0245, 0.0427, 0.1536, 0.2745, 0.3082];

if rmsPlotsOn == 1
    figure
    plot(sigs,ckf_RMS,'-rx','markersize',15,'linewidth',2)
    PlotBoi2('\sigma (10^x) [km/s^2]','Post Fit Residuals RMS (CKF)',16)
    
    figure
    plot(sigs,ekf_RMS,'-rx','markersize',15,'linewidth',2)
    PlotBoi2('\sigma (10^x) [km/s^2]','Post Fit Residuals RMS (EKF)',16)
    
    figure
    plot(sigs,ckfPos_RMS,'-rx','markersize',15,'linewidth',2)
    PlotBoi2('\sigma (10^x) [km/s^2]','3D RMS - Position (CKF)',16)
    
    figure
    plot(sigs,ckfVel_RMS,'-rx','markersize',15,'linewidth',2)
    PlotBoi2('\sigma (10^x) [km/s^2]','3D RMS - Velocity (CKF)',16)
    
    figure
    plot(sigs,ckf_RMS+ckfPos_RMS+ckfVel_RMS,'-rx','markersize',15,'linewidth',2)
    PlotBoi2('\sigma (10^x) [km/s^2]','Post Fit Residuals RMS + Pos/Vel RMS (CKF)',16)
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------
% %% Numerical Integrating Truth Trajectory
% ------------------------------------------------------------------------
function [ dY ] = u_J2_J3_OrbitIntegrator(t,Y,u,RE,J2,J3)
    dY = zeros(6,1);

    %%% Unpack the state vector (ECI)
    x = Y(1);
    y = Y(2);
    z = Y(3);
    dy = Y(4:6); % Satellite Velocity, km/s

    %%% Output the derivative of fthe state
    dY(1:3) = dy; % km/s
    % -EQM(4)
    dY(4) = (3*J2*RE^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*x)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
    % -EQM(5)
    dY(5) = (3*J2*RE^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*y)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
    % -EQM(6)
    dY(6) = (9*J2*RE^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)) - (3*J3*RE^3*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*z)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*RE^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J3*RE^3*u*z^2)/(x^2 + y^2 + z^2)^(7/2) - (35*J3*RE^3*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2));
end

% ------------------------------------------------------------------------
%%% Numerical Integrating Reference Orbit and STM
% ------------------------------------------------------------------------
function [ dY ] = Hw3_STMInt(t,Y,u,RE,J2,J3,dragOn,aDrag_RIC)
    %%% Size dY to fit state and all of reshaped (n^2,1) STM
    dY = zeros(6+6^2,1);
    
    %%% Unpack state
    x = Y(1);
    y = Y(2);
    z = Y(3);
    dx = Y(4);
    dy = Y(5);
    dz = Y(6);

    %%% Reshape (n^2,1) stm to (n,n)
    stm = reshape(Y(7:end),6,6);

    %%% Build A matrix and evaluate at current state
    A = zeros(6,6);
    A(1:3,4:6) = eye(3,3);
    % Asym(4,1)
    A(4,1) = -(u/(x^2 + y^2 + z^2)^(3/2) - (3*u*x^2)/(x^2 + y^2 + z^2)^(5/2) - (3*J2*RE^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J3*RE^3*u*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*x^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (35*J3*RE^3*u*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*x^2*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*x^2*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*x^2*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(4,2)
    A(4,2) = -((15*J2*RE^2*u*x*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*y)/(x^2 + y^2 + z^2)^(5/2) + (105*J3*RE^3*u*x*y*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*x*y*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*x*y*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(4,3)
    A(4,3) = -((45*J2*RE^2*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*x)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*x*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*x*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(5,1)
    A(5,1) = -((15*J2*RE^2*u*x*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*y)/(x^2 + y^2 + z^2)^(5/2) + (105*J3*RE^3*u*x*y*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*x*y*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*x*y*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(5,2)
    A(5,2) = -(u/(x^2 + y^2 + z^2)^(3/2) - (3*u*y^2)/(x^2 + y^2 + z^2)^(5/2) - (3*J2*RE^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J3*RE^3*u*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*y^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (35*J3*RE^3*u*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*y^2*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*y^2*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*y^2*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(5,3)
    A(5,3) = -((45*J2*RE^2*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*y*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*y*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*y*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(6,1)
    A(6,1) = -((45*J2*RE^2*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*x)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*x*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*x*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(6,2)
    A(6,2) = -((45*J2*RE^2*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*y*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*y*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*y*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
    % Asym(6,3)
    A(6,3) = -(u/(x^2 + y^2 + z^2)^(3/2) - (3*u*z^2)/(x^2 + y^2 + z^2)^(5/2) - (9*J2*RE^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (75*J3*RE^3*u*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (45*J2*RE^2*u*z^2)/(x^2 + y^2 + z^2)^(7/2) - (105*J2*RE^2*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2)) + (175*J3*RE^3*u*z^3)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*z^5)/(2*(x^2 + y^2 + z^2)^(11/2)));

    %%% Calculate new STM
    stm_dot = A*stm;

    %%% Creating X-dot
    % Spacecraft velocities
    dY(1:3) = [dx; dy; dz];
    
    %%% Creating drag term
    if dragOn == 1
        %%% Developing rotation matrix
        r_temp = [x, y, z];
        v_temp = [dx, dy, dz];
        Rbar = r_temp./norm(r_temp);
        h_temp = cross(r_temp,v_temp);
        Cbar = h_temp./norm(h_temp);
        Ibar = cross(Cbar,Rbar);
        B_RIC2ECI = [Rbar; Ibar; Cbar]';
        
        aDrag = B_RIC2ECI*aDrag_RIC;
    else
        aDrag = zeros(3,1);
    end
    
    
    %%% Using u, J2, and J3 terms as dynamics
    % -EQM(4)
    dY(4) = (3*J2*RE^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*x)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) - aDrag(1);
    % -EQM(5)
    dY(5) = (3*J2*RE^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*y)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) - aDrag(2);
    % -EQM(6)
    dY(6) = (9*J2*RE^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)) - (3*J3*RE^3*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*z)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*RE^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J3*RE^3*u*z^2)/(x^2 + y^2 + z^2)^(7/2) - (35*J3*RE^3*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2)) - aDrag(3);

    % Filling in reshaped (6^2,1) STM to state
    dY(7:end) = reshape(stm_dot,36,1);
end

% ------------------------------------------------------------------------
%%% DPDX (x, y, z only) Calculator
% ------------------------------------------------------------------------
function [ DPDxyz ] = DPDxyzCalculator(n, Rs, Rstat)
    DPDxyz = (Rs(n) - Rstat(n)) / norm(Rs - Rstat);
end

% ------------------------------------------------------------------------
%%% DdPDX (x, y, z only) Calculator
% ------------------------------------------------------------------------
function [DdPDXxyz] = DdPDxyzCalculator(n, Rs, Rstat, Vs, Vstat)
    p = norm(Rs-Rstat);
    DdPDXxyz = (Vs(n) - Vstat(n))/p - (dot((Rs-Rstat),(Vs-Vstat))*(Rs(n)-Rstat(n)))/(p^3);
end

% ------------------------------------------------------------------------
%%% DdPDX (vx, vy, vz only) Calculator
% ------------------------------------------------------------------------
function [DdPDvxvyvz] = DdPDvxvyvzCalculator(n, Rs, Rstat)
    DdPDvxvyvz = (Rs(n) - Rstat(n)) / norm(Rs - Rstat);
end

% ------------------------------------------------------------------------
%%% Range Rate Calculator
% ------------------------------------------------------------------------
function [rangeRate] = rangeRateCalculator(Rs, Rstat, Vs, Vstat)
    rangeRate = dot((Rs-Rstat),(Vs-Vstat))/norm(Rs-Rstat);
end




