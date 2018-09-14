clear
clc
close all
addpath('../../bin')
addpath('../../bin/FilterPlots')
addpath('../../bin/FilterFuncs')
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
processNoiseLimit = 250; % # of no-measurements before turning off process noise
Q_RIC_On          = 0; % Form Q_RIC then Q (Q_ECI)?
dragOn            = 0; % Include drag in truth dynamics

%%% Batch Options
runBatch   = 0;
BatchPlots = 0;

%%% CKF Options
runCKF   = 1;
CKFPlots = 1;

%%% SRIF Options
runSRIF = 1;

%%% RMS Options
rmsPlotsOn = 0;

%%% Smoothing
smoothingOn = 0;

%%% RMS stuff
RMSstuffOn = 0;

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

diff(-U_J2J3,x)
diff(-U_J2J3,y)
diff(-U_J2J3,z)
989
return
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
%%% Creating matrix with relevant station parameters for each time (ECI)
% ------------------------------------------------------------------------
statPrms_ECI = zeros(length(Times), 6);
for k = 1:length(Times)
    %%% If measurement is from station 1
    if range1_LEO_n(k) ~= 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) == 0
        statPrms_ECI(k,:) = [rStat1_ECI(k,:), vStat1_ECI(k,:)];
    %%% If measurement is from station 2
    elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) ~= 0 && range3_LEO_n(k) == 0
        statPrms_ECI(k,:) = [rStat2_ECI(k,:), vStat2_ECI(k,:)];
    %%% If measurement is from station 3
    elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) ~= 0
        statPrms_ECI(k,:) = [rStat3_ECI(k,:), vStat3_ECI(k,:)];
    end
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
%%% Batch (2a) (p196)
% ========================================================================
if runBatch == 1
%%% Preparing for State and STM propagation
XBar0_B = XBar0; 
IC_Batch = IC_ref;

%%% For storing deviation estimates
xHats_Batch = zeros(length(Times),6);

refSTMs_Batch = zeros(6,6,length(Times));
refStates_Batch = zeros(length(Times),6);
%%% Run Batch
fprintf('Batch Started: %3.1f sec\n',toc)
for bl = 1:10 % Outer Batch Loop
if bl > 1
    IC_Batch(1:6) = IC_Batch(1:6) + xHat0_Batch;
end
%%% Generating Reference State (STM)
[rT,refStatesSTM_Batch] = ode45(@Hw3_STMInt,Times,IC_Batch,options,uE,RE,J2,0,0); % J3 = 0

for k = 1:length(rT)
    %%% Reference States
    refStates_Batch(k,:) = refStatesSTM_Batch(k,1:6);
    
    %%% STMs
    refSTMs_Batch(:,:,k) = reshape(refStatesSTM_Batch(k,7:end),6,6);
end

%%% Start doing Batch things
lamda = zeros(6,6);
N = zeros(6,1);
l = 1; % Index for counting measurements
yis_B = zeros(2,length(Times)); % Prefit residual
Hs = zeros(2,6,length(Times));

if bl > 1 % If this is not the first loop of the outer batch
    XBar0_B = XBar0_B - xHat0_Batch;
    N = lamda*XBar0_B;
end
    
for k = 1:length(Times) % Inner Batch Loop
    
    if noisyRangeMeas(k) ~= 0 % If we have a measurement at this time step
        l = l + 1; % counting measurement times
        
        %%% Renaming Current Measurements
        Yi = [noisyRangeMeas(k); noisyRangeRateMeas(k)];
        
        %%% Calculating Htilde
        [ Htilde ] = calculateHtilde(refStates_Batch(k,1:6), statPrms_ECI(k,:));
        
        %%% Calculating pre-fit residuals
        [ yi ] = calculatePreFitResiduals(refStates_Batch(k,1:6), statPrms_ECI(k,:), Yi);
        
        %%% Processesing Observation
        H = Htilde*refSTMs_Batch(:,:,k); % 2x6
        lamda = lamda + H'*inv(R)*H; % 6x6
        N = N + H'*inv(R)*yi; % 6x1
        Hs(:,:,k) = H; % 2x6x1
        yis_B(:,k) = yi; % 2x1
    end
end

xHat0_Batch = linsolve(lamda,N);
P0 = inv(lamda);

%%% Calculate Post-Fit Residuals
[ eis_B ] = calcPostFitResiduals(length(Times), yis_B, Hs, xHat0_Batch);

%%% Calculate RMS value
m = 2*l; % Total number of measurements
[ RMS_Batch_bl ] = calculateRMS(length(Times), eis_B, R, m);
RMS_Batch(bl) = RMS_Batch_bl;


if bl > 1
    if  abs(RMS_Batch(bl)-RMS_Batch(bl-1))/RMS_Batch(bl-1) < 1e-8
        break % break if RMS %change tolerance met
    end
end

end

P_est_Batch = zeros(6,6,length(Times));
P_est_Batch(:,:,1) = P0;
for k = 2:length(Times)
    refSTM_Batch = refSTMs_Batch(:,:,k) * inv(refSTMs_Batch(:,:,k-1));
    P_est_Batch(:,:,k) = refSTM_Batch*P_est_Batch(:,:,k-1)*refSTM_Batch';
end

est_Batch = refStates_Batch; % for state estimates
%%% Computing estimate error
error_Batch = truthStateSTM(:,1:6) - est_Batch;

%%% Batch Plots of Residuals and Errors
fs = 14; % fontsize
if BatchPlots == 1
    PlotPostFitResiduals(Times, eis_B, pError, dpError, fs)
    PlotEstimateError(Times, error_Batch, P_est_Batch,fs)
end
fprintf('Batch Finished: %3.1f sec\n',toc)

end

% ========================================================================
%%% CKF (p203)
% ========================================================================
%%% Matrices to store results
refStatesSTM_CKF = zeros(length(Times),42);
refSTMs_CKF = zeros(6,6,length(Times));
xHats_CKF = zeros(length(Times),6);
P_est_CKF = zeros(6,6,length(Times));
P_bar_CKF = zeros(6,6,length(Times));
yis_C = zeros(2,length(Times)); % Prefit residual
eis_C = zeros(2,length(Times)); % Posfit residual
est_CKF = zeros(length(Times),6); % for state estimates
Htildes = zeros(2,6,length(Times));

P_est_CKF(:,:,1) = PBar0;
P_bar_CKF(:,:,1) = PBar0;
xHats_CKF(1,:) = dev;
Pi = PBar0;

refStatesSTM_CKF(1,:) = IC_ref';
est_CKF(1,:) = refStatesSTM_CKF(1,1:6) + xHats_CKF(1,:);
refSTMs_CKF(:,:,1) = reshape(refStatesSTM_CKF(1,7:end),6,6);

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
    refSTMs_CKF(:,:,k) = refSTM_CKF;
    refSTM_CKF =  refSTM_CKF * inv(reshape(refStatesSTM_CKF(k-1,7:end),6,6));
        
    if noisyRangeMeas(k) ~= 0 % If there is a measurement
        if ismember(Times(k),measTimes) ~= 1
            warning('CKF Time index error\n')
        end
        noMeasCount = 0;
        %%% Renaming Current Measurements
        Yi = [noisyRangeMeas(k); noisyRangeRateMeas(k)];
        
        %%% Time Update
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
        xBari = refSTM_CKF*xHats_CKF(k-1,:)';
        PBari = refSTM_CKF*P_est_CKF(:,:,k-1)*refSTM_CKF' + gamma*Q*gamma';
        P_bar_CKF(:,:,k) = PBari;
        
        %%% Calculating Htilde
        [ Htilde ] = calculateHtilde(refStatesSTM_CKF(k,1:6), statPrms_ECI(k,:));
        Htildes(:,:,k) = Htilde;
        
        %%% Calculating pre-fit residuals
        [ yi ] = calculatePreFitResiduals(refStatesSTM_CKF(k,1:6), statPrms_ECI(k,:), Yi);
        
        %%% Kalman Gain Matrix
        Ki = PBari*Htilde'*inv(Htilde*PBari*Htilde' + R);
        
        %%% Measurement Update
        xhat = xBari + Ki*(yi - Htilde*xBari);
        Pi = (eye(6) - Ki*Htilde)*PBari;
        
        xHats_CKF(k,:) = xhat;
        P_est_CKF(:,:,k) = Pi;
        yis_C(:,k) = yi;

    elseif noisyRangeMeas(k) == 0 % If no measurement
        noMeasCount = noMeasCount + 1;
        
        %%% Propagate P and xHat with STM
        xHats_CKF(k,:) = refSTM_CKF*xHats_CKF(k-1,:)';
        
        if noMeasCount > processNoiseLimit % No process noise for big gaps
            P_est_CKF(:,:,k) = refSTM_CKF*P_est_CKF(:,:,k-1)*refSTM_CKF';
            P_bar_CKF(:,:,k) = P_est_CKF(:,:,k);
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
            P_bar_CKF(:,:,k) = P_est_CKF(:,:,k);
        end
    end
    est_CKF(k,:) = refStatesSTM_CKF(k,1:6) + xHats_CKF(k,:);
end
clear r_temp v_temp h_temp

%%% Calculate Post-Fit Residuals
[ eis_C ] = calcPostFitResiduals(length(Times), yis_C, Htildes, xHats_CKF');

%%% Calculate RMS value
m = 2*length(measTimes); % Total number of measurements
[ RMS_CKF ] = calculateRMS(length(Times), eis_C, R, m);

%%% Computing estimate error
error_CKF = truthStateSTM(:,1:6) - est_CKF(:,:);

%%% CKF Plots
if CKFPlots == 1
    PlotPostFitResiduals(Times, eis_C, pError, dpError, fs)
    PlotEstimateError(Times, error_CKF, P_est_CKF,fs)
end

fprintf('CKF Finished: %3.1f\n',toc)
end

% ========================================================================
%%% Smoothing (p 246)
% ========================================================================
if smoothingOn == 1
fprintf('Smoothing Started: %3.1f sec\n',toc)
%%% Preparing
estSmooth_CKF = zeros(size(est_CKF));
estSmooth_CKF(end,:) = est_CKF(end,:);
xHats_lk = zeros(size(xHats_CKF));
Ps_lk = zeros(size(P_est_CKF));

%%% CKF Smoothing
xHats_lk(end,:) = xHats_CKF(end,:)';
Ps_lk(:,:,end) = P_est_CKF(:,:,end);
for l = length(Times):-1:2
    k = l-1;
    
    phi_tkp1_tk = refSTMs_CKF(:,:,k+1) * inv(refSTMs_CKF(:,:,k)); % Phi(tk+1, tk)

    Sk = P_est_CKF(:,:,k) * phi_tkp1_tk' * inv(P_bar_CKF(:,:,k+1));
    xHats_lk(k,:) = (xHats_CKF(k,:)' + Sk*(xHats_lk(k+1,:)' - phi_tkp1_tk*xHats_CKF(k,:)'))'; 
    Ps_lk(:,:,k) = P_est_CKF(:,:,k) + Sk*(Ps_lk(:,:,k+1) - P_bar_CKF(:,:,k+1))*Sk';
    
    
%     phi_tkp1_tk = refSTMs_CKF(:,:,k+1) * inv(refSTMs_CKF(:,:,k)); % Phi(tk+1, tk)
% 
%     Sk = P_est_CKF(:,:,k) * phi_tkp1_tk' * inv(P_bar_CKF(:,:,k+1));
%     xHats_lk(k,:) = (xHats_CKF(k,:)' + Sk*(xHats_lk(k+1,:)' - phi_tkp1_tk*xHats_CKF(k,:)'))'; 
%     Ps_lk(:,:,k) = P_est_CKF(:,:,k) + Sk*(Ps_lk(:,:,k+1) - P_bar_CKF(:,:,k+1))*Sk';
    


%     if processNoiseOn == 1 && noisyRangeMeas(k) ~= 0 % process noise (at tk)
%         Sk = P_est_CKF(:,:,k) * phi_tkp1_tk' * inv(P_bar_CKF(:,:,k+1));
%         xHats_lk(k,:) = (xHats_CKF(k,:)' + Sk*(xHats_lk(k+1,:)' - phi_tkp1_tk*xHats_CKF(k,:)'))'; 
%         Ps_lk(:,:,k) = P_est_CKF(:,:,k) + Sk*(Ps_lk(:,:,k+1) - P_bar_CKF(:,:,k+1))*Sk';
%     else % No process noise
% %         phi_tk_tl = refSTMs_CKF(:,:,k) * inv(refSTMs_CKF(:,:,l)); % Phi(tk, tl)
%         Sk = inv(phi_tkp1_tk);
%         xHats_lk(k,:) = Sk * xHats_lk(k+1,:)';
%         Ps_lk(:,:,k) = Sk * Ps_lk(:,:,k+1) * Sk';
%     end
    
    estSmooth_CKF(k,:) = refStatesSTM_CKF(k,1:6) + xHats_lk(k,:);
end


errorSmooth_CKF = truthStateSTM(:,1:6) - estSmooth_CKF;
PlotEstimateError(Times, errorSmooth_CKF, Ps_lk,fs)


fprintf('Smoothing Finished: %3.1f sec\n',toc)
end

%%% Specifically for problem 1c
% error_SmoothingBatch = error_Batch - errorSmooth_CKF;
% figure;
% subplot(3,2,1);
% plot(Times,error_SmoothingBatch(:,1))
% PlotBoi2('','X Position Error, km',fs)
% xlim([Times(1) Times(end)])
% subplot(3,2,3);
% plot(Times,error_SmoothingBatch(:,2))
% PlotBoi2('','Y Position Error, km',fs)
% xlim([Times(1) Times(end)])
% subplot(3,2,5);
% plot(Times,error_SmoothingBatch(:,3))
% PlotBoi2('Times, sec','Z Position Error, km',fs)
% xlim([Times(1) Times(end)])
% subplot(3,2,2);
% plot(Times,error_SmoothingBatch(:,4))
% PlotBoi2('','X Velocity Error, km',fs)
% xlim([Times(1) Times(end)])
% subplot(3,2,4);
% plot(Times,error_SmoothingBatch(:,5))
% PlotBoi2('','Y Velocity Error, km',fs)
% xlim([Times(1) Times(end)])
% subplot(3,2,6);
% plot(Times,error_SmoothingBatch(:,6))
% PlotBoi2('Times, sec','Z Velocity Error, km',fs)
% xlim([Times(1) Times(end)])

% ------------------------------------------------------------------------
%%% RMS work (problem 1d)
% ------------------------------------------------------------------------
if RMSstuffOn == 1
sumXPos_CKF = 0;
sumXVel_CKF = 0;
sumXPos_SCKF = 0;
sumXVel_SCKF = 0;

sumYPos_CKF = 0;
sumYVel_CKF = 0;
sumYPos_SCKF = 0;
sumYVel_SCKF = 0;

sumZPos_CKF = 0;
sumZVel_CKF = 0;
sumZPos_SCKF = 0;
sumZVel_SCKF = 0;
for k = 1:length(Times)
    sumXPos_CKF = sumXPos_CKF + norm(error_CKF(k,1));
    sumXVel_CKF = sumXVel_CKF + norm(error_CKF(k,4));
    sumXPos_SCKF = sumXPos_SCKF + norm(errorSmooth_CKF(k,1));
    sumXVel_SCKF = sumXVel_SCKF + norm(errorSmooth_CKF(k,4));
    
    sumYPos_CKF = sumYPos_CKF + norm(error_CKF(k,2));
    sumYVel_CKF = sumYVel_CKF + norm(error_CKF(k,5));
    sumYPos_SCKF = sumYPos_SCKF + norm(errorSmooth_CKF(k,2));
    sumYVel_SCKF = sumYVel_SCKF + norm(errorSmooth_CKF(k,5));
    
    sumZPos_CKF = sumZPos_CKF + norm(error_CKF(k,3));
    sumZVel_CKF = sumZVel_CKF + norm(error_CKF(k,6));
    sumZPos_SCKF = sumZPos_SCKF + norm(errorSmooth_CKF(k,3));
    sumZVel_SCKF = sumZVel_SCKF + norm(errorSmooth_CKF(k,6));
end

RMSXPos_CKF = (sumXPos_CKF/m)^(1/2);
RMSYPos_CKF = (sumYPos_CKF/m)^(1/2);
RMSZPos_CKF = (sumZPos_CKF/m)^(1/2);
RMSXVel_CKF = (sumXVel_CKF/m)^(1/2);
RMSYVel_CKF = (sumYVel_CKF/m)^(1/2);
RMSZVel_CKF = (sumZVel_CKF/m)^(1/2);

RMSXPos_SCKF = (sumXPos_SCKF/m)^(1/2);
RMSYPos_SCKF = (sumYPos_SCKF/m)^(1/2);
RMSZPos_SCKF = (sumZPos_SCKF/m)^(1/2);
RMSXVel_SCKF = (sumXVel_SCKF/m)^(1/2);
RMSYVel_SCKF = (sumYVel_SCKF/m)^(1/2);
RMSZVel_SCKF = (sumZVel_SCKF/m)^(1/2);

RMS_Comp_CKF = [RMSXPos_CKF; RMSYPos_CKF; RMSZPos_CKF; RMSXVel_CKF; RMSYVel_CKF; RMSZVel_CKF];
RMS_Comp_SCKF = [RMSXPos_SCKF; RMSYPos_SCKF; RMSZPos_SCKF; RMSXVel_SCKF; RMSYVel_SCKF; RMSZVel_SCKF];

for k = 1:6
    fprintf('-----k\n')
    pc = ((RMS_Comp_SCKF(k) - RMS_Comp_CKF(k))/RMS_Comp_CKF(k)) * 100;
end
    
end
% ------------------------------------------------------------------------
%%% SRIF (pg 510)
% ------------------------------------------------------------------------
if runSRIF == 1
    fprintf('SRIF Started: %3.1f sec\n',toc)
%%% Matrices to store results
yks_SRIF = zeros(2,length(Times)); % Prefit residual
Htildes_SRIF = zeros(2,6,length(Times));
bs_SRIF = zeros(6,length(Times));
Rs_SRIF = zeros(6,6,length(Times));
Ps_SRIF = zeros(6,6,length(Times));
xHats_SRIF = zeros(6,length(Times));
est_SRIF = zeros(length(Times),6); % for state estimates

%%% Integrate ref trajectory and STM
[rT,refStatesSTM_SRIF] = ode45(@Hw3_STMInt,Times,IC_ref,options,uE,RE,J2,0,0,0); % J3 = 0, drag off

%%% Setting initial values
Ps_SRIF(:,:,1) = PBar0;
xHats_SRIF(:,1) = dev;
Rs_SRIF(:,:,1) = chol(Ps_SRIF(:,:,1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not entirely sure what's going on here
% Initial R matrix is zeros
est_SRIF(1,:) = refStatesSTM_SRIF(1,1:6) + xHats_SRIF(:,1)';

% Yk = [noisyRangeMeas(1); noisyRangeRateMeas(1)];
% yk = calculatePreFitResiduals(refStatesSTM_SRIF(1,1:6), statPrms_ECI(1,:), Yk);
% yks_SRIF(:,1) = yk;
fprintf('Getting error when entering 2nd arc (k=268, T=38040)\n')
fprintf('line 434... Rs_SRIF(k-1) is zeros for first data on new arc\n')
for k = 2:length(Times)
    %%% Grabbing current stm ... phi(tk, tk-1)
    stm = reshape(refStatesSTM_SRIF(k,7:end),6,6) * inv(reshape(refStatesSTM_SRIF(k-1,7:end),6,6));

    if noisyRangeMeas(k) ~= 0 % If there is a measurement
        %%%%% Time update
        x_bar = stm*xHats_SRIF(:,k-1);
        R_bar = Rs_SRIF(:,:,k-1) * inv(stm);
        [R_bar] = houseHolderTrans(R_bar);
        b_bar = R_bar * x_bar;

        %%%%% Compute yi, Hitilde
        %%% Current Measurements
        Yk = [noisyRangeMeas(k); noisyRangeRateMeas(k)];
        %%% Calculating Htilde
        [ Htilde ] = calculateHtilde(refStatesSTM_SRIF(k,1:6), statPrms_ECI(k,:));
        %%% Calculating pre-fit residuals
        [ yk ] = calculatePreFitResiduals(refStatesSTM_SRIF(k,1:6), statPrms_ECI(k,:), Yk);

        %%%%% Whitening measurements
        Vk = chol(R);
        yk = Vk\yk; % same as inv(Vk) * yk
        Htilde = Vk\Htilde; % same as inv(Vk) * Htilde
        yks_SRIF (:,k) = yk;
        Htildes_SRIF(:,:,k) = Htilde;
        
        %%%%% Householder Transformation (pg 316-323)
        T = [R_bar, b_bar;...
             Htilde, yk];
        [T] = houseHolderTrans(T);
        n_temp = size(T,2)-1;
        m_temp = size(T,1)-n_temp;
        Rk = T(1:n_temp,1:n_temp);
        bk = T(1:n_temp,end);
        ek = T(n_temp+1:end,end);
        %%% Storing some values
        bs_SRIF(:,k) = bk;
        Rs_SRIF(:,:,k) = Rk;

        %%%%% Measurement Update
        xHats_SRIF(:,k) = Rk\bk;

    elseif noisyRangeMeas(k) == 0 % If no measurement
        %%% Propagate P and xHat with STM
        xHats_SRIF(:,k) = stm*xHats_SRIF(:,k-1);
        R_bar = Rs_SRIF(:,:,k-1) * inv(stm);
        [R_bar] = houseHolderTrans(R_bar);
        Rs_SRIF(:,:,k) = R_bar;
    end
    
    %%%%% Updating Estimate
    est_SRIF(k,:) = refStatesSTM_SRIF(k,1:6) + xHats_SRIF(:,k)';
    
end

%%% Computing Covariances
for k = 2:length(Times)
    Ps_SRIF(:,:,k) = inv(Rs_SRIF(:,:,k)) * inv(Rs_SRIF(:,:,k))';
end

%%% Computing estimate error
error_SRIF = truthStateSTM(:,1:6) - est_SRIF(:,:);

%%% Calculate Post-Fit Residuals
[ eks_SRIF ] = calcPostFitResiduals(length(Times), yks_SRIF, Htildes_SRIF, xHats_SRIF);

%%% Un-Whiten
for k = 1:length(Times)
    eks_SRIF(:,k) = Vk*eks_SRIF(:,k);
end

PlotPostFitResiduals(Times, eks_SRIF, pError, dpError, fs)
PlotEstimateError(Times, error_SRIF, Ps_SRIF,fs)

end % end SRIF section





































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