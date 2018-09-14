clear
clc
% close all
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

%%% SRIF Options
runSRIF = 1;

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
        yk = Vk\yk;
        Htilde = Vk\Htilde;
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


PlotPostFitResiduals(Times, eks_SRIF, pError, dpError, fs)
PlotEstimateError(Times, error_SRIF, Ps_SRIF,fs)

end % end SRIF section


































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

