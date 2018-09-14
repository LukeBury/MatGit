clear
clc
close all
addpath('../../bin')
addpath('../../Research/Europa_Hopper_Analysis/ProjectBin')
addpath('../../bin/FilterPlots')
addpath('../../bin/FilterFuncs')
addpath('../../bin/BPlaneFuncs')
tic

% ------------------------------------------------------------------------
%%% Plot/Run Options (0 = no/off, 1 = yes/on)
% ------------------------------------------------------------------------
%%% Reiterate ref trajectory and SRIF??
rerunRefAndSRIF = 0;

%%% Recalculate Station Paramters?
recalculateStationParameters = 0;

%%% Smooth SRIF data?
smoothingOn = 1;

%%% Movie
movieOn = 0;

%%% Setting time period
t0 = 0; % (sec)
dt = 60; % (sec)
tf = 2.5e7;

% %%% Setting time vector
% time = t0:dt:tf; % (sec)

% ------------------------------------------------------------------------
%%% Determining Equations of Motion
% ------------------------------------------------------------------------
syms x y z dx dy dz xs ys zs uE uS Cr Pphi Am c AU pb1 pb2 pb3 dpb1 dpb2 dpb3

r = [x; y; z]; % km
rs = [xs; ys; zs]; % km
Rs = rs-r;

a2B = -uE*r/(norm(r)^3);
a3B = -uS*((r-rs)/(norm(Rs)^3) + rs/(norm(rs)^3));

%%% SRP equations of motion (acceleration)
aSRP = (-Cr*Pphi*Am*(AU^2)/((norm(Rs)^3)*c))*Rs;

%%% Defining Equations of Motion, State Variables, and A Matrix
EQM = [dx; dy; dz; a2B(1) + a3B(1) + aSRP(1); a2B(2) + a3B(2) + aSRP(2); a2B(3) + a3B(3) + aSRP(3); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
state = [x; y; z; dx; dy; dz; Cr; uE; uS; AU; pb1; pb2; pb3; dpb1; dpb2; dpb3];
Asym = jacobian(EQM, state);

% syms xS yS zS dxS dyS dzS pb dpb
% dr = [dx; dy; dz];
% rStat = [xS; yS; zS];
% vStat = [dxS; dyS; dzS];
% p = norm(r - rStat) + pb;
% dp = dot((r-rStat),(dr-vStat))/p + dpb;
% 
% state_mod = [x; y; z; dx; dy; dz; Cr; Am; uE; uS; AU; pb; pb; pb; dpb; dpb; dpb];
% jacobian([p;dp],state_mod)

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
rE = 6378.1363; % km
rE_SOI = 925000; % km
uE = 3.98600432896939e5; % km^3/s^2
uS = 132712440017.987; % km^3/s^2
wE = [0, 0, 7.29211585275553e-5]; % rad/s
AU =  149597870.7; % km
Am = 0.01e-6; % km^2/kg 
Pphi = 1357; % kg/s^3
c = 299792.458; % km/s
JD0 =  2456296.25; % days
Cr = 1.0;
pb1 = 0;
pb2 = 0;
pb3 = 0;
dpb1 = 0;
dpb2 = 0;
dpb3 = 0;

% ------------------------------------------------------------------------
%%% Measurements
% ------------------------------------------------------------------------
%%% Measurement noise / uncertainty
pError = 0.005; % km
dpError = 0.0000005; % km/s
R = [pError^2 0; 0 dpError^2];

%%% Loading Measurements
% Time since Epoch, DSS34 Range (km), DSS65 Range, DSS13 Range, DSS34 Range-Rate (km/sec), DSS65 Range-Rate, DSS13 Range-Rate 
dataSet = csvread('Project1b_Obs.txt',2,0);
% dataSet = csvread('Project1b_Obs.txt',2,0);

%%% Checking for One Station per Time
for k = 1:size(dataSet,1)
    sum = 0;
    if isnan(dataSet(k,2)) == 1
        sum = sum+1;
    end
    if isnan(dataSet(k,3)) == 1
        sum = sum+1;
    end
    if isnan(dataSet(k,4)) == 1
        sum = sum+1;
    end
    if sum ~= 2
        warning('Some time has multiple measurements (dataA)')
    end
end

%%% Collapsing data matrix
data = zeros(size(dataSet,1),4);
for k = 1:size(dataSet,1)
    if isnan(dataSet(k,2)) == 0
        data(k,:) = [1,dataSet(k,1),dataSet(k,2),dataSet(k,5)];
    elseif isnan(dataSet(k,3)) == 0
        data(k,:) = [2,dataSet(k,1),dataSet(k,3),dataSet(k,6)];
    elseif isnan(dataSet(k,4)) == 0
        data(k,:) = [3,dataSet(k,1),dataSet(k,4),dataSet(k,7)];
    end
end

%%% Building Time Vector
% %%% Setting time vector
time = [t0:dt:18814168-28]; % (sec)
time = [time, 18814168:dt:tf-12]; 

% ------------------------------------------------------------------------
%%% Station Positions and Velocities
% ------------------------------------------------------------------------
if recalculateStationParameters == 1
    %%% Initial lat/lon/alts of Stations
    latLonAlt = [-35.398333, 148.981944, 0.691750;
        40.427222, -4.250556, 0.834539;
        35.247164, 243.205, 1.07114904]; % deg | deg | km

    %%% Calculating Intial Positions and Velocities of Stations in ECI
    statStates = zeros(6,3,length(time));
    for k = 1:3
        %%% Acquiring ECEF/ECI station position
        [r] = latlon2surfECEF(latLonAlt(k,1), latLonAlt(k,2), latLonAlt(k,3) + rE); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correct way of handling altitude?

        %%% Storing initial station positions and velocities
        statStates(1:3,k,1) = r'; % km
        statStates(4:6,k,1) = cross(wE',r'); % km/s
    end; clear r

    %%% Propagating Station States in ECI
    for k = 2:length(time)
        %%% Calculating current rotation angle
        th = wE(3)*time(k); % rad 
        for s = 1:3
            statStates(1:3,s,k) = R3(statStates(1:3,s,1),th); % km
            statStates(4:6,s,k) = R3(statStates(4:6,s,1),th); % km/s
        end
    end
    
    %%% Saving Station Parameters
    save('statStates_B.mat','statStates');
else
    %%% Loading Station Parameters
    load('statStates_B.mat')
end

% ------------------------------------------------------------------------
%%% Setting Initial Satellite State                               
% ------------------------------------------------------------------------
r0_Sat = [-274096770.76544; -92859266.4499061; -40199493.6677441]; % km
v0_Sat = [32.6704564599943; -8.93838913761049; -3.87881914050316]; % km/s

X0 = [r0_Sat; v0_Sat; Cr; uE; uS; AU; pb1; pb2; pb3; dpb1; dpb2; dpb3]; % (km, km/s)
stm0 = eye(16);
IC_ref = [X0; reshape(stm0,256,1)];

% ------------------------------------------------------------------------
%%% Preparing for Filter
% ------------------------------------------------------------------------
r_sig = 100; % km
v_sig = 0.1; % km/s
Cr_sig = 0.1;
uE_sig = 10; % km^3/s^2
uS_sig = 500; % km^3/s^2
AU_sig = 10; % km
pb_sig = 100; % km
dpb_sig = 50; % km/s


PBar0 = blkdiag(eye(3)*(r_sig^2),eye(3)*(v_sig^2),Cr_sig^2,uE_sig^2,uS_sig^2,AU_sig^2,eye(3)*pb_sig^2,eye(3)*dpb_sig^2);
Vk = chol(R);

% ------------------------------------------------------------------------
%%% Propagating States
% ------------------------------------------------------------------------
%%% Setting integrator options
tol = 1E-13;
options = odeset('RelTol',tol,'AbsTol',tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rerunRefAndSRIF == 1
    iterN = 5;
for iter = 1:iterN
    covExplosion = 0;
    fprintf('SRIF Process Started: %3.1f sec\n',toc)
    if iter == 1
        time = [t0:dt:50*86400];
    elseif iter == 2
        time = [t0:dt:100*86400];
    elseif iter == 3
        time = [t0:dt:150*86400];
    elseif iter == 4
        time = [t0:dt:(18814168-28)];
    elseif iter == 5
        time = [t0:dt:18814168-28]; % (sec)
        time = [time, 18814168:dt:data(end,2)]; 
%     elseif iter == iterN
%         time = [t0:dt:18814168-28]; % (sec)
%         time = [time, 18814168:dt:tf-12];
    end
    
%%% Propagating the Ref State (Traj & STM)
[Times,refStatesSTMs] = ode45(@proj1Integrator_B,time,IC_ref,options,Pphi,c,JD0,Am);

%%% Calculating Sun Positions
rSuns = zeros(length(Times),3);
for k = 1:length(Times)
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,3,'EME2000'); % km, km/s, km^3/s^2
    rSuns(k,:) = -r;
end

% ------------------------------------------------------------------------
%%% SRIF
% ------------------------------------------------------------------------
%%% Matrices to store results
yks_SRIF = zeros(2,length(Times)); % Prefit residual
Htildes_SRIF = zeros(2,16,length(Times));
bs_SRIF = zeros(16,length(Times));
Rs_SRIF = zeros(16,16,length(Times));
Ps_SRIF = zeros(16,16,length(Times));
xHats_SRIF = zeros(16,length(Times));
est_SRIF = zeros(length(Times),16); % for state estimates

%%% Setting initial values
Ps_SRIF(:,:,1) = PBar0;
Rs_SRIF(:,:,1) = chol(Ps_SRIF(:,:,1)); 
est_SRIF(1,:) = refStatesSTMs(1,1:16) + xHats_SRIF(:,1)';
Yk = [data(1,3); data(1,4)]; % First measurement set
yk = calculatePreFitResiduals(refStatesSTMs(1,1:6), statStates(:,data(1,1),1)', Yk); 
yks_SRIF(:,1) = Vk\yk;

% kEnd = 290000;
kEnd = length(Times);
for k = 2:kEnd
    warning('off','MATLAB:nearlySingularMatrix')
    %%% Grabbing current stm ... phi(tk, tk-1)
    stm = reshape(refStatesSTMs(k,17:end),16,16) * inv(reshape(refStatesSTMs(k-1,17:end),16,16));

    if any(data(:,2)==Times(k)) == 1 % If there is a measurement
        %%% Find row of measurement
        m = find(data(:,2)==Times(k));
        
        %%%%%%% Time update
        if Times(k) > (18814168-100) && covExplosion == 0
            R_bar = (blkdiag(Rs_SRIF(1:6,1:6,1)./10000, Rs_SRIF(7:16,7:16,k-1)) * inv(reshape(refStatesSTMs(k,17:end),16,16)));
            covExplosion = 1;
        else
            R_bar = Rs_SRIF(:,:,k-1) * inv(stm);
        end
        x_bar = stm*xHats_SRIF(:,k-1);
        [R_bar] = houseHolderTrans(R_bar);
        b_bar = R_bar * x_bar;

        %%%%%%% Compute yi, Hitilde
        %%% Current Measurements
        Yk = [data(m,3); data(m,4)];
        %%% Calculating Htilde
        stationNum = data(m,1);
        refRangeBias = refStatesSTMs(k,11+stationNum);
        [ Htilde ] = calculateHtilde_B(refStatesSTMs(k,1:6), statStates(:,stationNum,k)', stationNum, refRangeBias);
        %%% Calculating pre-fit residuals
        [ yk ] = calculatePreFitResiduals(refStatesSTMs(k,1:6), statStates(:,data(m,1),k)', Yk);

        %%%%%%% Whitening measurements
        yk = Vk\yk; % same as inv(Vk) * yk
        Htilde = Vk\Htilde; % same as inv(Vk) * Htilde
        yks_SRIF (:,k) = yk;
        Htildes_SRIF(:,:,k) = Htilde;
        
        %%%%%%% Householder Transformation (pg 316-323)
        T = [R_bar, b_bar;...
             Htilde, yk];
        [T] = houseHolderTrans(T);
        n_temp = size(T,2)-1;
        Rk = T(1:n_temp,1:n_temp);
        bk = T(1:n_temp,end);
        ek = T(n_temp+1:end,end);
        %%% Storing some values
        bs_SRIF(:,k) = bk;
        Rs_SRIF(:,:,k) = Rk;

        %%%%%%% Measurement Update
        xHats_SRIF(:,k) = Rk\bk;
        
    else % If no measurement
        %%% Propagate P and xHat with STM
        xHats_SRIF(:,k) = stm*xHats_SRIF(:,k-1);
        R_bar = Rs_SRIF(:,:,k-1) * inv(stm);
        [R_bar] = houseHolderTrans(R_bar);
        Rs_SRIF(:,:,k) = R_bar;
    end
    
    %%%%%%% Updating Estimate
    est_SRIF(k,:) = refStatesSTMs(k,1:16) + xHats_SRIF(:,k)';
    
end

%%% For iterating, recalculate a priori using estimate at time of last measurement
% mEnd = find(Times==data(end,2)); % Last measurement time index
x0_iter = inv(reshape(refStatesSTMs(end,17:end),16,16)) * xHats_SRIF(:,end); % Mapping deviation estimate back to time 0
X0 = refStatesSTMs(1,1:16) + x0_iter';

%%% For iterating, reset a priori
IC_ref = [X0'; reshape(stm0,256,1)];


end % end iter

%%% Saving Ref State for Loading Purposes
save('RefAndSRIFResults_B.mat','kEnd','Times','refStatesSTMs','rSuns','yks_SRIF','Htildes_SRIF','bs_SRIF','Ps_SRIF','xHats_SRIF','est_SRIF','Rs_SRIF','Vk');

else % Skip ref and SRIF and just reload them
    load('RefAndSRIFResults_B.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerunRefAndSRIF == 1
    warning('off','MATLAB:nearlySingularMatrix')
    %%% Computing Covariances
    for k = 2:kEnd 
        Ps_SRIF(:,:,k) = inv(Rs_SRIF(:,:,k)) * inv(Rs_SRIF(:,:,k))';
    end

    %%% Calculate Post-Fit Residuals
    [ eks_SRIF ] = calcPostFitResiduals(length(Times), yks_SRIF, Htildes_SRIF, xHats_SRIF);

    %%% Un-Whiten
    for k = 1:kEnd 
        eks_SRIF(:,k) = Vk*eks_SRIF(:,k);
    end
    
    %%% Save Data
    save('PsAndEs_B.mat','Ps_SRIF','eks_SRIF')
else
    
    load('PsAndEs_B.mat')
end

%%% Storing Final Meas. Time for Plots, Baby
mEnd = find(Times==data(end,2));

%%% Plotting 3-Sig Covariances (to end of measurements)
lw = 1.5;
figure
subplot(3,2,1); hold all; ylim([-.1 .1])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(1,1,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(1,1,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','X Error, km',14)
subplot(3,2,3); hold all; ylim([-10 10])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(2,2,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(2,2,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','Y Error, km',14)
subplot(3,2,5); hold all; ylim([-60 60])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(3,3,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(3,3,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('Times, days','Z Error, km',14)
subplot(3,2,2); hold all; ylim([-2e-5 2e-5])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(4,4,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(4,4,1:mEnd))),'--r','linewidth',lw)
legend('3\sigma Covariance')
PlotBoi2('','dX Error, km',14)
subplot(3,2,4); hold all; ylim([-2e-5 2e-5])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(5,5,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(5,5,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','dY Error, km',14)
subplot(3,2,6); hold all; ylim([-2e-5 2e-5])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(6,6,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(6,6,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','dZ Error, km',14)

figure
subplot(5,2,1); hold all; ylim([-2*3*sqrt(Ps_SRIF(7,7,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(7,7,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(7,7,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(7,7,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','Cr Error',14)
subplot(5,2,2); hold all; ylim([-2*3*sqrt(Ps_SRIF(8,8,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(8,8,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(8,8,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(8,8,1:mEnd))),'--r','linewidth',lw)
legend('3\sigma Covariance')
PlotBoi2('','\muEarth Error',14)
subplot(5,2,3); hold all; ylim([-2*3*sqrt(Ps_SRIF(9,9,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(9,9,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(9,9,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(9,9,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','\muSun Error',14)
subplot(5,2,4); hold all; ylim([-2*3*sqrt(Ps_SRIF(10,10,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(10,10,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(10,10,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(10,10,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','AU Error',14)
subplot(5,2,5); hold all; ylim([-2*3*sqrt(Ps_SRIF(11,11,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(11,11,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(11,11,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(11,11,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','\rho Bias 1 Error',14)
subplot(5,2,7); hold all; ylim([-2*3*sqrt(Ps_SRIF(12,12,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(12,12,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(12,12,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(12,12,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','\rho Bias 2 Error',14)
subplot(5,2,9); hold all; ylim([-2*3*sqrt(Ps_SRIF(13,13,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(13,13,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(13,13,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(13,13,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','\rho Bias 3 Error',14)
subplot(5,2,6); hold all; ylim([-2*3*sqrt(Ps_SRIF(14,14,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(14,14,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(14,14,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(14,14,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','\rhoRate Bias 1 Error',14)
subplot(5,2,8); hold all; ylim([-2*3*sqrt(Ps_SRIF(15,15,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(15,15,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(15,15,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(15,15,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('','\rhoRate Bias 2 Error',14)
subplot(5,2,10); hold all; ylim([-2*3*sqrt(Ps_SRIF(16,16,mEnd/2+.5)) 2*3*sqrt(Ps_SRIF(16,16,mEnd/2+.5))])
plot(Times(1:mEnd)./86400,3*sqrt(squeeze(Ps_SRIF(16,16,1:mEnd))),'--r','linewidth',lw)
plot(Times(1:mEnd)./86400,-3*sqrt(squeeze(Ps_SRIF(16,16,1:mEnd))),'--r','linewidth',lw)
PlotBoi2('Times, days','\rhoRate Bias 3 Error',14)

%%% Plotting Pre-Fit Residuals
PlotPreFitResiduals(Times(1:mEnd), yks_SRIF(:,1:mEnd), pError, dpError, 16) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not full time

%%% Plotting Post-Fit Residuals
PlotPostFitResiduals(Times(1:mEnd), eks_SRIF(:,1:mEnd), pError, dpError, 16) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not full time


% ------------------------------------------------------------------------
%%% B-Plane
% ------------------------------------------------------------------------
figure(22); hold all
figure(33); hold all
%%% Earth Equator 
th = 0:.01:2*pi;
t = rE * cos(th);
r = rE * sin(th);
plot(t, r,'b','linewidth',1.5);
for dk = 1:4
%%% Integrate day 50/100/150/200 day estimate
d50 = 86400*50*dk/60+1; % Index of day 50
%%% Setting Initial Condition
IC_50 = [est_SRIF(d50,:)';reshape(eye(16),256,1)];
%%% Integrating from 50/100/150/200 days to B-Plane Crossing
[Time50,est50] = ode45(@proj1Integrator_B,[Times(d50):dt:tf],IC_50,options,Pphi,c,JD0,Am);

%%% Finding 3SOI Time
for k = 1:length(Time50)
    if norm(est50(k,1:3)) < 3*rE_SOI
        % Store Velocity as vInfinity
        vInf = est_SRIF(k,4:6); % km/s
        
        % Store the index at which vInfinity is determined
        vInfk = k;
        
        %%% Running BPlane Function, but only to get LTOF from vInf time
        [ BdT, BdR, DCM_vInf, LTOF ] = BPlaneCalculator(est50(vInfk,1:6), uE, vInf);
        break
    end
end

%%% Finding B-Plane Cross time
for k = vInfk:length(Time50) % From Time of vInf Calculation Onward..
    if Time50(k) > (Time50(vInfk) + LTOF) % If Time of Crossing is Found
        crossk = k; % Time index when B-Plane is crossed
        
        %%% Running BPlane Function to get BdT, BdR, and DCM_Cross
        [ BdT, BdR, DCM_Cross, LTOF ] = BPlaneCalculator(est50(crossk,1:6), uE, vInf);
        break
    end
end

%%% Rotating Estimated State at Crossing into B-Plane
rCross_B = DCM_Cross*est50(crossk,1:3)';

%%% Propagating Covariance from 50/100/150/200 days to B-Plane Cross Time
PCross = reshape(est50(crossk,17:end),16,16) * Ps_SRIF(:,:,d50) * reshape(est50(crossk,17:end),16,16)';

%%% Rotating Estimated Covariance at Crossing into B-Plane
PCross_B = DCM_Cross * PCross(1:3,1:3) * DCM_Cross';

%%% Plotting Error Ellipse
figure(22)
error_ellipse(PCross_B(2:3,2:3),'conf',0.9973)
PlotBoi2('T-axis, km','R-axis, km',16)

%%% Plotting Crossing Position
figure(33)
plot(rCross_B(2),rCross_B(3),'x','markersize',10)
PlotBoi2('T-axis, km','R-axis, km',16)
axis equal
end
figure(22); legend('50','100','150','200')
figure(33); legend('Earth','50','100','150','200')











fprintf('Done: %3.1f sec\n',toc)
% ------------------------------------------------------------------------
%%% Movie of Motion
% ------------------------------------------------------------------------
if movieOn == 1
figure
for i = 1:length(Times)
    if rem(i,1000) == 0
        clf
        hold all
        trackWidth = 2;
        lineWidth = 1.5;
        %%% Plotting Current Location Marker
        plot3(refStatesSTMs(i,1),refStatesSTMs(i,2),refStatesSTMs(i,3),'rX','linewidth',trackWidth,'markersize',10);
        
        %%% Plotting Past Trajectory
        plot3(refStatesSTMs(1:i-1,1),refStatesSTMs(1:i-1,2),refStatesSTMs(1:i-1,3),'-m')

        %%% Earth Equator 
        th = 0:.01:2*pi;
        x = 1000*rE * cos(th);
        y = 1000*rE * sin(th);
        plot(x, y,'b','linewidth',lineWidth);

        %%% Sun Equator 
        th = 0:.01:2*pi;
        x = 3000*rE * cos(th) + rSuns(i,1);
        y = 3000*rE * sin(th) + rSuns(i,2);
        plot(x, y,'b','linewidth',lineWidth);
        
        PlotBoi3('X, km','Y, km','Z, km',16)
        axis equal
        view(-20,60)

        drawnow limitrate
    end
end
end

