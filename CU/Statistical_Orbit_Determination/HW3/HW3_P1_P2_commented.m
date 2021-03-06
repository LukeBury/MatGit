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
devOn = 1; % Turn on deviation from ref trajectory
J3On = 0; % Turn on J3 for truth trajectory

%%% Filter on?
runCKF = 0;
runEKF = 0;

%%% Filter Results
BatchPlots = 1;
CKFPlots = 0;
EKFPlots = 0;

fprintf('Turn on J3 and deviation!\n')

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
%%% Setting time period
t0 = 0; % (sec)
dt = 10; % (sec)
tf = 24*3600; % (sec)
time = t0:dt:tf; % (sec)

%%% Setting integrator options
tol = 1E-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% %%% Propagating LEO
% [Times,truthStateSTM] = ode45(@u_J2_J3_OrbitIntegrator,time,X0_Sat,options,uE,RE,J2,J3);
%%% Preparing for State and STM propagation
stm0 = eye(6);
IC_truth = [X0_Sat; reshape(stm0,36,1)];

%%% Propagating the State (STM)
[Times,truthStateSTM] = ode45(@Hw3_STMInt,time,IC_truth,options,uE,RE,J2,J3);

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

pError = 0.001;
dpError = 0.000001;
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
%%% Processing the State Transition Matrix
% ------------------------------------------------------------------------
% %%% Preparing for State and STM propagation
% stm0 = eye(6);
% IC_truth = [X0_Sat; reshape(stm0,36,1)];
% 
% %%% Propagating the State (STM)
% [Times,truthStateSTM] = ode45(@Hw3_STMInt,time,IC_truth,options,uE,RE,J2,J3);

%%% Reshaping ref states and STMs
States_truth = zeros(length(Times),6);
STMs_truth = zeros(6,6,length(Times));
for k = 1:length(Times)
    %%% Reference States
    States_truth(k,:) = truthStateSTM(k,1:6);
    
    %%% STMs
    STMs_truth(:,:,k) = reshape(truthStateSTM(k,7:end),6,6);
end

% ========================================================================
%%% Batch (2a) (p196)
% ========================================================================
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
measTimes = [0];
for k = 1:length(Times)
   if noisyRangeMeas(k) ~= 0
       measTimes = [measTimes; Times(k)];
   end
end

%%% Setting up Initial Conditions

PBar0 = blkdiag(eye(3)*(pError^2),eye(3)*(dpError^2));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Square the errors
if devOn == 1
    dev = [.003; .001; .002; .000003; .000001;.000002];
else
    dev = zeros(6,1);
end
XBar0_ref = States_truth(1,:)' + dev;

R = [pError^2 0; 0 dpError^2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Reference Trajectory %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preparing for State and STM propagation
stm0 = eye(6);
IC_ref = [XBar0_ref; reshape(stm0,36,1)];
XBar0_B = XBar0_ref; 

%%% Propagating the State (STM)
[refTimes,refStateSTM] = ode45(@Hw3_STMInt,time,IC_ref,options,uE,RE,J2,0); % J3 = 0

%%% Reshaping ref states and STMs
refStates_Batch = zeros(length(refTimes),6);
refSTMs_Batch = zeros(6,6,length(refTimes));
for k = 1:length(refTimes)
    %%% Reference States
    refStates_Batch(k,:) = refStateSTM(k,1:6);
    
    %%% STMs
    refSTMs_Batch(:,:,k) = reshape(refStateSTM(k,7:end),6,6);
end

%%% Initiallizing matrices to store values 
% xHats_Batch = zeros(length(Times),6);


%%% Run Batch
fprintf('Batch Started: %3.1f sec\n',toc)
for bl = 1:10 % Outer Batch Loop
lamda = zeros(6,6);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Testing
N = zeros(6,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Testing
l = 1; % Index for counting measurements
yis_B = zeros(2,length(Times)); % Prefit residual
Hs = zeros(2,6,length(Times));

if bl > 1 % If this is not the first loop of the outer batch
    XBar0_B = XBar0_B - xHat0;
    N = lamda*XBar0_B;
end
    
for k = 1:length(Times) % Inner Batch Loop
    
    if noisyRangeMeas(k) ~= 0 % If we have a measurement at this time step
        l = l + 1; % counting measurement times
        
        %%% Renaming Current Measurements
        Yi = [noisyRangeMeas(k); noisyRangeRateMeas(k)];

        %%% Generating Htilde and Reading Observation
        Htilde = zeros(2,6); % Observation-State Matrix
        
        %%% If measurement is from station 1
        if range1_LEO_n(k) ~= 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) == 0
            %%% DP/DX
            Htilde(1,1) = DPDxyzCalculator(1, refStates_Batch(k,1:3), rStat1_ECI(k,:));
            Htilde(1,2) = DPDxyzCalculator(2, refStates_Batch(k,1:3), rStat1_ECI(k,:));
            Htilde(1,3) = DPDxyzCalculator(3, refStates_Batch(k,1:3), rStat1_ECI(k,:));
            %%% DdP/DX
            Htilde(2,1) = DdPDxyzCalculator(1, refStates_Batch(k,1:3), rStat1_ECI(k,:), refStates_Batch(k,4:6), vStat1_ECI(k,:));
            Htilde(2,2) = DdPDxyzCalculator(2, refStates_Batch(k,1:3), rStat1_ECI(k,:), refStates_Batch(k,4:6), vStat1_ECI(k,:));
            Htilde(2,3) = DdPDxyzCalculator(3, refStates_Batch(k,1:3), rStat1_ECI(k,:), refStates_Batch(k,4:6), vStat1_ECI(k,:));
            Htilde(2,4) = DdPDvxvyvzCalculator(1, refStates_Batch(k,1:3), rStat1_ECI(k,:));
            Htilde(2,5) = DdPDvxvyvzCalculator(2, refStates_Batch(k,1:3), rStat1_ECI(k,:));
            Htilde(2,6) = DdPDvxvyvzCalculator(3, refStates_Batch(k,1:3), rStat1_ECI(k,:));
            %%% Prefit Residual
            yi = [Yi(1) - norm(refStates_Batch(k,1:3) - rStat1_ECI(k,:));...
                Yi(2) -  rangeRateCalculator(refStates_Batch(k,1:3), rStat1_ECI(k,:), refStates_Batch(k,4:6), vStat1_ECI(k,:))];
            
        %%% If measurement is from station 2
        elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) ~= 0 && range3_LEO_n(k) == 0
            %%% DP/DX
            Htilde(1,1) = DPDxyzCalculator(1, refStates_Batch(k,1:3), rStat2_ECI(k,:));
            Htilde(1,2) = DPDxyzCalculator(2, refStates_Batch(k,1:3), rStat2_ECI(k,:));
            Htilde(1,3) = DPDxyzCalculator(3, refStates_Batch(k,1:3), rStat2_ECI(k,:));
            %%% DdP/DX
            Htilde(2,1) = DdPDxyzCalculator(1, refStates_Batch(k,1:3), rStat2_ECI(k,:), refStates_Batch(k,4:6), vStat2_ECI(k,:));
            Htilde(2,2) = DdPDxyzCalculator(2, refStates_Batch(k,1:3), rStat2_ECI(k,:), refStates_Batch(k,4:6), vStat2_ECI(k,:));
            Htilde(2,3) = DdPDxyzCalculator(3, refStates_Batch(k,1:3), rStat2_ECI(k,:), refStates_Batch(k,4:6), vStat2_ECI(k,:));
            Htilde(2,4) = DdPDvxvyvzCalculator(1, refStates_Batch(k,1:3), rStat2_ECI(k,:));
            Htilde(2,5) = DdPDvxvyvzCalculator(2, refStates_Batch(k,1:3), rStat2_ECI(k,:));
            Htilde(2,6) = DdPDvxvyvzCalculator(3, refStates_Batch(k,1:3), rStat2_ECI(k,:));
            %%% Prefit Residual
            yi = [Yi(1) - norm(refStates_Batch(k,1:3) - rStat2_ECI(k,:));...
                Yi(2) -  rangeRateCalculator(refStates_Batch(k,1:3), rStat2_ECI(k,:), refStates_Batch(k,4:6), vStat2_ECI(k,:))];
        
        %%% If measurement is from station 3
        elseif range1_LEO_n(k) == 0 && range2_LEO_n(k) == 0 && range3_LEO_n(k) ~= 0
            %%% DP/DX
            Htilde(1,1) = DPDxyzCalculator(1, refStates_Batch(k,1:3), rStat3_ECI(k,:));
            Htilde(1,2) = DPDxyzCalculator(2, refStates_Batch(k,1:3), rStat3_ECI(k,:));
            Htilde(1,3) = DPDxyzCalculator(3, refStates_Batch(k,1:3), rStat3_ECI(k,:));
            %%% DdP/DX
            Htilde(2,1) = DdPDxyzCalculator(1, refStates_Batch(k,1:3), rStat3_ECI(k,:), refStates_Batch(k,4:6), vStat3_ECI(k,:));
            Htilde(2,2) = DdPDxyzCalculator(2, refStates_Batch(k,1:3), rStat3_ECI(k,:), refStates_Batch(k,4:6), vStat3_ECI(k,:));
            Htilde(2,3) = DdPDxyzCalculator(3, refStates_Batch(k,1:3), rStat3_ECI(k,:), refStates_Batch(k,4:6), vStat3_ECI(k,:));
            Htilde(2,4) = DdPDvxvyvzCalculator(1, refStates_Batch(k,1:3), rStat3_ECI(k,:));
            Htilde(2,5) = DdPDvxvyvzCalculator(2, refStates_Batch(k,1:3), rStat3_ECI(k,:));
            Htilde(2,6) = DdPDvxvyvzCalculator(3, refStates_Batch(k,1:3), rStat3_ECI(k,:));
            %%% Prefit Residual
            yi = [Yi(1) - norm(refStates_Batch(k,1:3) - rStat3_ECI(k,:));...
                Yi(2) -  rangeRateCalculator(refStates_Batch(k,1:3), rStat3_ECI(k,:), refStates_Batch(k,4:6), vStat3_ECI(k,:))];

        else
            warning('Trouble parsing noisy station measurements\n')
        end

        %%% DP/Dvx = DP/Dvy = DP/Dvz = 0
        Htilde(1,4) = 0;
        Htilde(1,5) = 0;
        Htilde(1,6) = 0;

        %%% Processesing Observation
        H = Htilde*refSTMs_Batch(:,:,k); % 2x6
        lamda = lamda + H'*inv(R)*H; % 6x6
        N = N + H'*inv(R)*yi; % 6x1
        Hs(:,:,k) = H; % 2x6x1
        yis_B(:,k) = yi; % 2x1
    end    
end

xHat0 = linsolve(lamda,N);
P0 = inv(lamda);
temp = 0;
eis_B = zeros(2,length(Times));
for k = 1:length(Times)
    if noisyRangeMeas(k) ~= 0 % If there was a measurement at this time
        ei = yis_B(:,k) - Hs(:,:,k)*xHat0; % 2x1
        temp = temp + ei'*R*ei; % 1x1
        eis_B(:,k) = ei; % 2x1
    end
end
m = 2*l; % Total number of measurements
RMS(bl) = (temp/m)^(1/2);

refStates_Batch(1,:) = refStates_Batch(1,:) + xHat0';

if bl > 1
    if  abs(RMS(bl)-RMS(bl-1))/RMS(bl-1) < 1e-8
        break % break if RMS %change tolerance met
    end
end

end

P_est_Batch = zeros(6,6,length(Times));
P_est_Batch(:,:,1) = P0;
for k = 2:length(Times)
    P_est_Batch(:,:,k) = refSTMs_Batch(:,:,k)*P_est_Batch(:,:,1)*refSTMs_Batch(:,:,k)';
end

%%% Batch Plots of Residuals and Errors
fs = 14; % fontsize
if BatchPlots == 1
    eis_B(eis_B == 0) = NaN;
    figure
    subplot(2,1,1); hold all
    plot(Times,eis_B(1,:),'-mx')
    plot([Times(1) Times(end)],[3*pError 3*pError],'--r')
    plot([Times(1) Times(end)],[-3*pError -3*pError],'--r')
    legend('Post fit','3\sigma Error')
    PlotBoi2('','Range Residuals, km',14)
    subplot(2,1,2); hold all
    plot(Times,eis_B(2,:),'-mx')
    plot([Times(1) Times(end)],[3*dpError 3*dpError],'--r')
    plot([Times(1) Times(end)],[-3*dpError -3*dpError],'--r')
    PlotBoi2('Times, sec','Range Rate Residuals, km/s',14)
    eis_B(isnan(eis_B)==1) = 0;
    
    figure
    subplot(3,2,1); hold all
    plot(Times(2:end), truthStateSTM(2:end,1) - (refStates_Batch(2:end,1))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_Batch(1,1,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_Batch(1,1,2:end))),'--r')
    PlotBoi2('','X Position Error, km',fs)
    subplot(3,2,3); hold all
    plot(Times(2:end), truthStateSTM(2:end,2) - (refStates_Batch(2:end,2))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_Batch(2,2,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_Batch(2,2,2:end))),'--r')
    PlotBoi2('','Y Position Error, km',fs)
    subplot(3,2,5); hold all
    plot(Times(2:end), truthStateSTM(2:end,3) - (refStates_Batch(2:end,3))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_Batch(3,3,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_Batch(3,3,2:end))),'--r')
    PlotBoi2('Time, sec','Z Position Error, km',fs)
    subplot(3,2,2); hold all
    plot(Times(2:end), truthStateSTM(2:end,4) - (refStates_Batch(2:end,4))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_Batch(4,4,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_Batch(4,4,2:end))),'--r')
    legend('Truth','3\sigma Covariance')
    PlotBoi2('','X Velocity Error, km/s',fs)
    subplot(3,2,4); hold all
    plot(Times(2:end), truthStateSTM(2:end,5) - (refStates_Batch(2:end,5))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_Batch(5,5,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_Batch(5,5,2:end))),'--r')
    PlotBoi2('','Y Velocity Error, km/s',fs)
    subplot(3,2,6); hold all
    plot(Times(2:end), truthStateSTM(2:end,6) - (refStates_Batch(2:end,6))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_Batch(6,6,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_Batch(6,6,2:end))),'--r')
    PlotBoi2('Time, sec','Z Velocity Error, km/s',fs)
end
fprintf('Batch Finished: %3.1f sec\n',toc)

% ========================================================================
%%% CKF (2b) (p203)
% ========================================================================
%%% Matrices to store results
refStatesSTM_CKF = zeros(length(Times),42);
xHats_CKF = zeros(length(Times),6);
P_est_CKF = zeros(6,6,length(Times));
yis_C = zeros(2,length(Times)); % Prefit residual
eis_C = zeros(2,length(Times)); % Posfit residual

PBar0 = blkdiag(eye(3)*pError^2,eye(3)*dpError^2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eye(6) / errors squared??
xBar0 = IC_ref(1:6);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% should this be full state?

P_est_CKF(:,:,1) = PBar0;
xHats_CKF(1,:) = dev;
Pi = PBar0;

refStatesSTM_CKF(1,:) = IC_ref';

xhat = dev;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% should this just be deviation?
PBari = PBar0;
R = [pError^2 0; 0 dpError^2];

if runCKF == 1
fprintf('CKF Started: %3.1f sec\n',toc)
%%% Integrate ref trajectory and STM
[rT,refStatesSTM_CKF] = ode45(@Hw3_STMInt,time,IC_ref,options,uE,RE,J2,0); % J3 = 0


l = 1; % counter for ref states

arcMeasurements = 0;
for k = 2:length(Times)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Times minus 1?
        
    %%% reshape new STM
    refSTM_CKF = reshape(refStatesSTM_CKF(k,7:end),6,6);
        
        
    if noisyRangeMeas(k) ~= 0 % If there is a measurement
        arcMeasurements = 1;
        if ismember(Times(k),measTimes) ~= 1
            warning('CKF Time index error\n')
        end
        l = l + 1; % Indexing measurement number
        
        %%% Renaming Current Measurements
        Yi = [noisyRangeMeas(k); noisyRangeRateMeas(k)];
        
        %%% Time Update
        xBari = refSTM_CKF*xhat;
        PBari = refSTM_CKF*Pi*refSTM_CKF';
        
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
        P_est_CKF(:,:,k) = Pi;
        yis_C(:,k) = yi;
        
        %%% Post fit residual
        temp = 0;
        ei = yis_C(:,k) - Htilde*xHats_CKF(k,:)'; % 2x1
        temp = temp + ei'*R*ei; % 1x1
        eis_C(:,k) = ei; % 2x1
        
    elseif noisyRangeMeas(k) == 0 % If no measurement

        %%% Propagate P and xHat with STM
        P_est_CKF(:,:,k) = refSTM_CKF*P_est_CKF(:,:,1)*refSTM_CKF';
        xHats_CKF(k,:) = refSTM_CKF*xHats_CKF(1,:)';
        
    end
    
end

%%% CKF Plots
if CKFPlots == 1
    eis_C(eis_C == 0) = NaN;
    figure
    subplot(2,1,1); hold all
    plot(Times,eis_C(1,:),'-mx')
    plot([Times(1) Times(end)],[3*pError 3*pError],'--r')
    plot([Times(1) Times(end)],[-3*pError -3*pError],'--r')
    legend('Post fit','3\sigma Error')
    PlotBoi2('','Range Residuals, km',14)
    subplot(2,1,2); hold all
    plot(Times,eis_C(2,:),'-mx')
    plot([Times(1) Times(end)],[3*dpError 3*dpError],'--r')
    plot([Times(1) Times(end)],[-3*dpError -3*dpError],'--r')
    PlotBoi2('Times, sec','Range Rate Residuals, km/s',14)
    eis_C(isnan(eis_C)==1) = 0;
    
    figure
    subplot(3,2,1); hold all
    plot(Times(2:end), truthStateSTM(2:end,1) - (refStatesSTM_CKF(2:end,1) + xHats_CKF(2:end,1))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_CKF(1,1,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_CKF(1,1,2:end))),'--r')
    PlotBoi2('','X Position Error, km',fs)
    subplot(3,2,3); hold all
    plot(Times(2:end), truthStateSTM(2:end,2) - (refStatesSTM_CKF(2:end,2) + xHats_CKF(2:end,2))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_CKF(2,2,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_CKF(2,2,2:end))),'--r')
    PlotBoi2('','Y Position Error, km',fs)
    subplot(3,2,5); hold all
    plot(Times(2:end), truthStateSTM(2:end,3) - (refStatesSTM_CKF(2:end,3) + xHats_CKF(2:end,3))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_CKF(3,3,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_CKF(3,3,2:end))),'--r')
    PlotBoi2('Time, sec','Z Position Error, km',fs)
    subplot(3,2,2); hold all
    plot(Times(2:end), truthStateSTM(2:end,4) - (refStatesSTM_CKF(2:end,4) + xHats_CKF(2:end,4))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_CKF(4,4,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_CKF(4,4,2:end))),'--r')
    legend('Truth','3\sigma Covariance')
    PlotBoi2('','X Velocity Error, km/s',fs)
    subplot(3,2,4); hold all
    plot(Times(2:end), truthStateSTM(2:end,5) - (refStatesSTM_CKF(2:end,5) + xHats_CKF(2:end,5))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_CKF(5,5,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_CKF(5,5,2:end))),'--r')
    PlotBoi2('','Y Velocity Error, km/s',fs)
    subplot(3,2,6); hold all
    plot(Times(2:end), truthStateSTM(2:end,6) - (refStatesSTM_CKF(2:end,6) + xHats_CKF(2:end,6))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_CKF(6,6,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_CKF(6,6,2:end))),'--r')
    PlotBoi2('Time, sec','Z Velocity Error, km/s',fs)
    
%     xlim([35370 35680]) % First arc
end

fprintf('CKF Finished: %3.1f\n',toc)
end

% ========================================================================
%%% EKF (2c) (p212)
% ========================================================================
%%% Matrices to store results
refStatesSTM_EKF = zeros(length(Times),42);
xHats_EKF = zeros(length(Times),6);
P_est_EKF = zeros(6,6,length(Times));
yis_E = zeros(2,length(Times)); % Prefit residual
eis_E = zeros(2,length(Times)); % Posfit residual

PBar0 = blkdiag(eye(3)*pError^2,eye(3)*dpError^2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eye(6) / errors squared??
xBar0 = IC_ref(1:6);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% should this be full state?

P_est_EKF(:,:,1) = PBar0;
xHats_EKF(1,:) = dev;
Pi = PBar0;

refStatesSTM_EKF(1,:) = IC_ref';

xhat = dev;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% should this just be deviation?
PBari = PBar0;
R = [pError^2 0; 0 dpError^2];

if runEKF == 1
l = 1; % counter for ref states
fprintf('EKF Started: %3.1f sec\n',toc)
arcMeasurements = 0; % Switch for when an arc begins
arcCounter = 0; % Counter to get 20 into arc as CKF
for k = 2:length(Times)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Times minus 1?
    %%% Integrate ref trajectory and STM since last measurement
        intTime = [Times(k-1), Times(k)];
        pvsStateSTM = refStatesSTM_EKF(k-1,:);
        [rT,rSSTM] = ode45(@Hw3_STMInt,intTime,pvsStateSTM,options,uE,RE,J2,0); % J3 = 0
        refStatesSTM_EKF(k,:) = rSSTM(end,:); % Store value from Times(k+1)
        
    %%% reshape new STM
    refSTM_EKF = reshape(refStatesSTM_EKF(k,7:end),6,6);
        
        
    if noisyRangeMeas(k) ~= 0 % If there is a measurement
        arcMeasurements = 1;
        arcCounter = arcCounter + 1;
        
        if ismember(Times(k),measTimes) ~= 1
            warning('EKF Time index error\n')
        end
        l = l + 1; % Indexing measurement number
        
        %%% Renaming Current Measurements
        Yi = [noisyRangeMeas(k); noisyRangeRateMeas(k)];
                
        %%% Time Update
%         xBari = refSTM_EKF*xhat;
        PBari = refSTM_EKF*Pi*refSTM_EKF';
        
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
        xhat = Ki*yi;
        Pi = (eye(6) - Ki*Htilde)*PBari;
        
        %%% Storing Values
        xHats_EKF(k,:) = xhat;
        P_est_EKF(:,:,k) = Pi;
        yis_E(:,k) = yi;
        
        %%% Post fit residual
        temp = 0;
        ei = yis_E(:,k) - Htilde*xHats_EKF(k,:)'; % 2x1
        temp = temp + ei'*R*ei; % 1x1
        eis_E(:,k) = ei; % 2x1
        
    elseif noisyRangeMeas(k) == 0 % If no measurement
        arcMeasurements = 0;
        arcCounter = 0;
        %%% Propagate P and xHat with STM
        P_est_EKF(:,:,k) = refSTM_EKF*P_est_EKF(:,:,1)*refSTM_EKF';
        xHats_EKF(k,:) = refSTM_EKF*xHats_EKF(1,:)';

%         P_est_EKF(:,:,k) = refSTM_EKF*P_est_EKF(:,:,k-1)*refSTM_EKF';
%         xHats_EKF(k,:) = refSTM_EKF*xHats_EKF(k-1,:)';

%         P_est_EKF(:,:,k) = P_est_EKF(:,:,k-1);
%         xHats_EKF(k,:) = xHats_EKF(k-1,:);
        
    end
    
    if arcCounter > 20 % If we've had at least 20 measurements on this arc
        refStatesSTM_EKF(k,1:6) = refStatesSTM_EKF(k,1:6) + xHats_EKF(k,:);
    end
    
end

%%% EKF Plots
if EKFPlots == 1
    eis_E(eis_E == 0) = NaN;
    figure
    subplot(2,1,1); hold all
    plot(Times,eis_E(1,:),'-mx')
    plot([Times(1) Times(end)],[3*pError 3*pError],'--r')
    plot([Times(1) Times(end)],[-3*pError -3*pError],'--r')
    legend('Post fit','3\sigma Error')
    PlotBoi2('','Range Residuals, km',14)
    subplot(2,1,2); hold all
    plot(Times,eis_E(2,:),'-mx')
    plot([Times(1) Times(end)],[3*dpError 3*dpError],'--r')
    plot([Times(1) Times(end)],[-3*dpError -3*dpError],'--r')
    PlotBoi2('Times, sec','Range Rate Residuals, km/s',14)
    eis_E(isnan(eis_E)==1) = 0;
    
    figure
    subplot(3,2,1); hold all
    plot(Times(2:end), truthStateSTM(2:end,1) - (refStatesSTM_EKF(2:end,1) + xHats_EKF(2:end,1))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_EKF(1,1,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_EKF(1,1,2:end))),'--r')
    PlotBoi2('','X Position Error, km',fs)
    subplot(3,2,3); hold all
    plot(Times(2:end), truthStateSTM(2:end,2) - (refStatesSTM_EKF(2:end,2) + xHats_EKF(2:end,2))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_EKF(2,2,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_EKF(2,2,2:end))),'--r')
    PlotBoi2('','Y Position Error, km',fs)
    subplot(3,2,5); hold all
    plot(Times(2:end), truthStateSTM(2:end,3) - (refStatesSTM_EKF(2:end,3) + xHats_EKF(2:end,3))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_EKF(3,3,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_EKF(3,3,2:end))),'--r')
    PlotBoi2('Time, sec','Z Position Error, km',fs)
    subplot(3,2,2); hold all
    plot(Times(2:end), truthStateSTM(2:end,4) - (refStatesSTM_EKF(2:end,4) + xHats_EKF(2:end,4))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_EKF(4,4,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_EKF(4,4,2:end))),'--r')
    legend('Truth','3\sigma Covariance')
    PlotBoi2('','X Velocity Error, km/s',fs)
    subplot(3,2,4); hold all
    plot(Times(2:end), truthStateSTM(2:end,5) - (refStatesSTM_EKF(2:end,5) + xHats_EKF(2:end,5))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_EKF(5,5,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_EKF(5,5,2:end))),'--r')
    PlotBoi2('','Y Velocity Error, km/s',fs)
    subplot(3,2,6); hold all
    plot(Times(2:end), truthStateSTM(2:end,6) - (refStatesSTM_EKF(2:end,6) + xHats_EKF(2:end,6))) %%% difference between truth and Xref + xHat
    plot(Times(2:end),3*sqrt(squeeze(P_est_EKF(6,6,2:end))),'--r')
    plot(Times(2:end),-3*sqrt(squeeze(P_est_EKF(6,6,2:end))),'--r')
    PlotBoi2('Time, sec','Z Velocity Error, km/s',fs)
    
%     xlim([35370 35680]) % First arc
end

fprintf('EKF Finished: %3.1f\n',toc)
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
function [ dY ] = Hw3_STMInt(t,Y,u,RE,J2,J3)
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

    %%% Using u, J2, and J3 terms as dynamics
    % -EQM(4)
    dY(4) = (3*J2*RE^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*x)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
    % -EQM(5)
    dY(5) = (3*J2*RE^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*y)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
    % -EQM(6)
    dY(6) = (9*J2*RE^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)) - (3*J3*RE^3*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*z)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*RE^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J3*RE^3*u*z^2)/(x^2 + y^2 + z^2)^(7/2) - (35*J3*RE^3*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2));

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




