clear
clc
close all
addpath('../../bin')
addpath('../../bin/FilterPlots')
addpath('../../bin/FilterFuncs')
tic

% ------------------------------------------------------------------------
%%% Determining Equations of Motion and Measurement Model
% ------------------------------------------------------------------------
syms xsE ysE zsE dxsE dysE dzsE xsJ ysJ zsJ u1 u2 real
syms xEJ yEJ zEJ dxEJ dyEJ dzEJ real
syms xJSun yJSun zJSun dxJSun dyJSun dzJSun real
syms xEaSun yEaSun zEaSun dxEaSun dyEaSun dzEaSun real
syms xstEa ystEa zstEa dxstEa dystEa dzstEa real

rsE = [xsE; ysE; zsE]; % Sat position in Europa-Inertial, km
rsJ = [xsJ; ysJ; zsJ]; % Sat position in Jupiter-Inertial, km
rEJ = rsJ-rsE; % Europa position in Jupiter-Inertial

a2B = -u2*rsE/(norm(rsE)^3); % km/s^2
a3B = -u1*(rsJ/(norm(rsJ)^3) - (rEJ)/(norm(rEJ)^3)); % km/s^2

%%% Defining Equations of Motion, State Variables, and A Matrix
EQM = [dxsE; dysE; dzsE; a2B(1) + a3B(1); a2B(2) + a3B(2); a2B(3) + a3B(3)];
state = [xsE; ysE; zsE; dxsE; dysE; dzsE];
Asym = jacobian(EQM, state);


%%% Measurement Model
vsE = [dxsE; dysE; dzsE];
rEJ = [xEJ; yEJ; zEJ]; vEJ = [dxEJ; dyEJ; dzEJ];
rJSun = [xJSun;yJSun;zJSun]; vJSun = [dxJSun;dyJSun;dzJSun];
rEaSun = [xEaSun;yEaSun;zEaSun]; vEaSun = [dxEaSun;dyEaSun;dzEaSun];
rstEa = [xstEa;ystEa;zstEa]; vstEa = [dxstEa;dystEa;dzstEa];

% ECI Range and Range-Rate from Europa-Cenetered-Inertial States
range = sqrt(dot((rsE + rEJ + rJSun - rEaSun - rstEa),(rsE + rEJ + rJSun - rEaSun - rstEa)));
rangeRate = dot((rsE + rEJ + rJSun - rEaSun - rstEa),(vsE + vEJ + vJSun - vEaSun - vstEa))/range;

measFunction = [range; rangeRate];
HtildeEqns = jacobian(measFunction,state);

clear

% ------------------------------------------------------------------------
%%% Plot/Run Options (0 = no/off, 1 = yes/on)
% ------------------------------------------------------------------------
measOn        = 1; % use measurements?
devOn         = 1;
linearityTest = 0;

runSRIF        = 0;
iterMax_SRIF   = 1;

runCKF             = 1;
iterMax_CKF        = 3;
processNoiseOn_CKF = 1;

plotAltitudes = 1;
% ------------------------------------------------------------------------
%%% System Parameters / Givens
% ------------------------------------------------------------------------
AU = 149597870; % km
d2r = pi/180;
day = 86400; % sec

%%% Europa
aEur = 671100; % semimajor axis
wEur = [0,0,2.047200349303344e-05]; % rad/s

%%% Gravitational Parameters
u1 = 126672520; % Jupiter, km^3 / s^2
u2 = 3203.413216; % Europa, km^3 / s^2

%%% Body Radii
rad1= 71492; % radius of Jupiter
rad2 = 1560.8; % radius of Europa

%%% Initial Julian Date
[JD0] = JD(2025,3,5,0,0);

%%% Earth
rE = 6378.1363; % km
wDay = [0,0,2*pi/day]; % rad/s

%%% Europa period (days)
pEur = 3.551181; 
%%% Setting normalized time vector
ti = 0;
dt = 60;
tf = 1*day;
time = (ti:dt:tf);

% ------------------------------------------------------------------------
%%% Measurements
% ------------------------------------------------------------------------
%%% Measurement noise / uncertainty
pError = (5e-3)*2; % km
dpError = 5e-7; % km/s
R = [pError^2 0; 0 dpError^2];

%%% Loading Measurements
% Time since Epoch, DSS34 Range (km), DSS65 Range, DSS13 Range, DSS34 Range-Rate (km/sec), DSS65 Range-Rate, DSS13 Range-Rate 
dataSet = csvread('finalProjectData_ECI2.txt');

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

%%% Loading Truth Trajectory
load('states_EurCI_truth.mat');
% ------------------------------------------------------------------------
%%% Station Positions and Velocities
% ------------------------------------------------------------------------
%%% Initial lat/lon/alts of Stations
latLonAlt = [-35.398333, 148.981944, 0.691750;
    40.427222, -4.250556, 0.834539;
    35.247164, 243.205, 1.07114904]; % deg | deg | km

%%% Calculating Intial Positions and Velocities of Stations in ECI
stations_ECI = zeros(6,3,length(time));
for k = 1:3
    %%% Acquiring ECEF/ECI station position
    [rStat_ECEF] = latlon2surfECEF(latLonAlt(k,1), latLonAlt(k,2), latLonAlt(k,3) + rE); % km

    %%% Storing initial station positions and velocities
    stations_ECI(1:3,k,1) = rStat_ECEF'; % km
    stations_ECI(4:6,k,1) = cross(wDay',rStat_ECEF'); % km/s
end; clear rStat_ECEF

%%% Propagating Station States in ECI
for k = 2:length(time)
    %%% Calculating current rotation angle
    th = wDay(3)*time(k); % rad 
    for s = 1:3
        stations_ECI(1:3,s,k) = R3(stations_ECI(1:3,s,1),th); % km
        stations_ECI(4:6,s,k) = R3(stations_ECI(4:6,s,1),th); % km/s
    end
end

% ------------------------------------------------------------------------
%%% Setting Initial Satellite State                               
% ------------------------------------------------------------------------
load('OrbitalInfo_nn.mat')
[rH0_EurCI, vH0_EurCI] = OE2ECI(a, e, i, raan, w, ta, u2);
rH0_EurCI = rH0_EurCI';
vH0_EurCI = vH0_EurCI';

%%% Deviation from truth
devR = [0.001; 0.007; 0.001];
devV = [1e-7; 1e-7; 1e-6];

stm0 = eye(6);
if devOn == 1
    X0_EurCI = [rH0_EurCI' + devR; vH0_EurCI' + devV];
else
    X0_EurCI = [rH0_EurCI'; vH0_EurCI'];
end
IC_BCR_ref = [X0_EurCI; reshape(stm0,36,1)];

% ------------------------------------------------------------------------
%%% Preparing for Filter
% ------------------------------------------------------------------------
r_sig = 2; % km
v_sig = 0.001; % km/s
 
PBar0 = blkdiag(eye(3)*(r_sig^2),eye(3)*(v_sig^2));
Vk = chol(R);

%%% Setting integrator options
tol = 1E-13;
options = odeset('Events',@impact_FP_nn,'RelTol',tol,'AbsTol',tol);

% ------------------------------------------------------------------------
%%% Testing Linearity
% ------------------------------------------------------------------------
%%% Doing this part regardless so I've got a "Times" to work with before SRIF
IC_BCR_ref_t = IC_BCR_ref;
[Times,x1] = ode45(@integrator_FP_nn,time,IC_BCR_ref_t,options,rad2,u1,u2,aEur,wEur);

if linearityTest == 1
IC_BCR_ref_t(1:6) = IC_BCR_ref_t(1:6) + [devR; devV];
[Times,x2] = ode45(@integrator_FP_nn,time,IC_BCR_ref_t,options,rad2,u1,u2,aEur,wEur);

dx0 = x2(1,1:6) - x1(1,1:6);

x2Hat = zeros(6,length(Times));
x2Hat(:,1) = dx0';
for k = 2:length(Times)
    % Phi(tk, tk-1)
    stm = reshape(x1(k,7:end),6,6) * inv(reshape(x1(k-1,7:end),6,6));
    
    x2Hat(:,k) = stm * (x2Hat(:,k-1));
end

x22 = x1(:,1:6) + x2Hat';

plot6StateError(Times,x2(:,1:6),x22)
subplot(3,2,1); title('x2 - x22')

plot6StateError(Times,x2(:,1:6),x1(:,1:6))
subplot(3,2,1); title('x2 - x1')

plot6StateError(Times,x1(:,1:6),x22)
subplot(3,2,1); title('x1 - x22')

figure; hold all
plot(Times,x1(:,1))
plot(Times,x2(:,1))
plot(Times,x22(:,1))
legend('x1','x2','x22')


end

% ------------------------------------------------------------------------
%%% Calculating other necessary states
% ------------------------------------------------------------------------
states_EJ = zeros(length(Times),6); 
states_EaSun = zeros(length(Times),6);
states_JSun = zeros(length(Times),6);

%%% Barycenter Inertial State (BCI)
for k = 1:length(Times)
    %%% Rotation Angle
    th = Times(k) * wEur(3); % rad
    
    %%% Europa in JCI
    states_EJ(k,1:3) = R3([aEur, 0, 0],th); % km
    states_EJ(k,4:6) = R3(cross(wEur,states_EJ(1,1:3)),th); % km/s
	
    %%% Jupiter relative to Sun Inertial
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,5,'EME2000'); % km, km/s, km^3/s^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    states_JSun(k,1:6) = [r',v'];
    
    %%% Earth relative to Sun Inertial
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,3,'EME2000'); % km, km/s, km^3/s^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    states_EaSun(k,1:6) = [r',v'];
end

% ========================================================================
%%% SRIF
% ========================================================================
if runSRIF == 1
    
for iter = 1:iterMax_SRIF
%%% Propagating the Ref State (Traj & STM)
[Times,refStatesSTMs_EurCI] = ode45(@integrator_FP_nn,time,IC_BCR_ref,options,rad2,u1,u2,aEur,wEur);

% ------------------------------------------------------------------------
%%% Preparing for SRIF
% ------------------------------------------------------------------------
%%% Matrices to store results
yks_SRIF = zeros(2,length(Times)); % Prefit residual
Htildes_SRIF = zeros(2,6,length(Times));
bs_SRIF = zeros(6,length(Times));
Rs_SRIF = zeros(6,6,length(Times));
Ps_SRIF = zeros(6,6,length(Times));
xHats_SRIF = zeros(6,length(Times));
est_SRIF = zeros(length(Times),6); % for state estimates

%%% Setting initial values
Ps_SRIF(:,:,1) = PBar0;
Rs_SRIF(:,:,1) = chol(Ps_SRIF(:,:,1)); 
est_SRIF(1,:) = refStatesSTMs_EurCI(1,1:6) + xHats_SRIF(:,1)';
if data(1,2) == 0
    Yk = [data(1,3); data(1,4)]; % First measurement set
    yk = calculatePreFitResiduals_EurCI(refStatesSTMs_EurCI(1,1:6), states_EJ(1,1:6), states_JSun(1,1:6), states_EaSun(1,1:6), stations_ECI(:,data(1,1),k)', Yk); 
    yks_SRIF(:,1) = Vk\yk;
end

kEnd = length(Times);
for k = 2:kEnd
    %%% Grabbing current stm ... phi(tk, tk-1)
    stm = reshape(refStatesSTMs_EurCI(k,7:end),6,6) * inv(reshape(refStatesSTMs_EurCI(k-1,7:end),6,6));

    if any(data(:,2)==Times(k)) == 1 && measOn == 1 % If there is a measurement (if the data time is in my time vector)
        %%% Find row of measurement
        m = find(data(:,2)==Times(k));
        
        %%%%%%% Time update
        x_bar = stm*xHats_SRIF(:,k-1);
        R_bar = Rs_SRIF(:,:,k-1) * inv(stm);
        [R_bar] = houseHolderTrans(R_bar);
        b_bar = R_bar * x_bar;

        %%%%%%% Compute yi, Hitilde
        %%% Current Measurements
        Yk = [data(m,3); data(m,4)];
        %%% Calculating Htilde
        [ Htilde ] = calculateHtilde_EurCI(refStatesSTMs_EurCI(k,1:6), states_EJ(k,1:6), states_JSun(k,1:6), states_EaSun(k,1:6), stations_ECI(:,data(m,1),k)');
        %%% Calculating pre-fit residuals
        [ yk ] = calculatePreFitResiduals_EurCI(refStatesSTMs_EurCI(k,1:6), states_EJ(k,1:6), states_JSun(k,1:6), states_EaSun(k,1:6), stations_ECI(:,data(m,1),k)', Yk);

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
    est_SRIF(k,:) = refStatesSTMs_EurCI(k,1:6) + xHats_SRIF(:,k)';
    
end

if iter < iterMax_SRIF
    %%%%% Prepare for next iteration
    %%% Propagating back to t0 w/ stm
    dev0_iter = inv(reshape(refStatesSTMs_EurCI(end,7:end),6,6)) * xHats_SRIF(:,end); % Mapping deviation estimate back to time 0
    X0_iter = refStatesSTMs_EurCI(1,1:6) + dev0_iter';
    
    IC_BCR_ref = [X0_iter'; reshape(stm0,36,1)];
    

end
end % iterations

%%% Computing Covariances
for k = 2:kEnd
    Ps_SRIF(:,:,k) = inv(Rs_SRIF(:,:,k)) * inv(Rs_SRIF(:,:,k))';
end

%%% Calculate Post-Fit Residuals
[ eks_SRIF ] = calcPostFitResiduals(length(Times), yks_SRIF, Htildes_SRIF, xHats_SRIF);

%%% Un-Whiten
for k = 1:kEnd
    eks_SRIF(:,k) = Vk*eks_SRIF(:,k);
    yks_SRIF(:,k) = Vk*yks_SRIF(:,k);
end
warning('Unwhitening both residuals... \n')

%%% Compute Altitudes
europaAltitudes_SRIF = zeros(length(Times),1); % km
for k = 1:length(Times)
    europaAltitudes_SRIF(k) = norm(est_SRIF(k,1:3)) - rad2; % km
end

%%% Calculating Error
error = est_SRIF - states_EurCI_truth;

%%% Plotting 3-Sig Covariances (to end of measurements)
lw = 1;
if round(Times(end)) ~= Times(end) % If there was a collision, don't plot last error b/c it's ugly
    pEnd = length(Times)-1;
else
    pEnd = length(Times);
end
figure
subplot(3,2,1); hold all; ylim([-.06 .06])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_SRIF(1,1,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_SRIF(1,1,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error(1:pEnd,1),'k')
PlotBoi2('','X Error, km',14)
subplot(3,2,3); hold all;ylim([-.08 .08])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_SRIF(2,2,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_SRIF(2,2,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error(1:pEnd,2),'k')
PlotBoi2('','Y Error, km',14)
subplot(3,2,5); hold all;ylim([-.1 .1])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_SRIF(3,3,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_SRIF(3,3,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error(1:pEnd,3),'k')
PlotBoi2('Time, hr','Z Error, km',14)
subplot(3,2,2); hold all;ylim([-5e-5 5e-5])
p1 = plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_SRIF(4,4,1:pEnd))),'--r','linewidth',lw);
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_SRIF(4,4,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error(1:pEnd,4),'k')
legend([p1],'3\sigma Covariance')
PlotBoi2('','dX Error, km/s',14)
subplot(3,2,4); hold all;ylim([-5e-5 5e-5])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_SRIF(5,5,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_SRIF(5,5,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error(1:pEnd,5),'k')
PlotBoi2('','dY Error, km/s',14)
subplot(3,2,6); hold all;ylim([-10e-5 10e-5])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_SRIF(6,6,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_SRIF(6,6,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error(1:pEnd,6),'k')
PlotBoi2('Time, hr','dZ Error, km/s',14)

%%% Plotting Pre-Fit Residuals
PlotPreFitResiduals(Times(1:pEnd), yks_SRIF(:,1:pEnd), data(:,1:2), pError, dpError, 16)

%%% Plotting Post-Fit Residuals
PlotPostFitResiduals(Times(1:pEnd), eks_SRIF(:,1:pEnd), data(:,1:2), pError, dpError, 16)

end % runSRIF

% ========================================================================
%%% CKF
% ========================================================================
if runCKF == 1
    
% ------------------------------------------------------------------------
%%% Process Noise
% ------------------------------------------------------------------------
pnSig = 10^-6;
%%% Process noise / uncertainty
if processNoiseOn_CKF == 1
    Q = [(pnSig)^2 0 0;
        0 (pnSig)^2 0;
        0 0 (pnSig)^2];

    gamma = dt*[eye(3)*(dt/2); eye(3)];
else
    Q = zeros(3,3);
    gamma = zeros(6,3);
end

%%% Resetting normalized time vector
ti = 0;
dt = 60;
tf = 1*day;
time = (ti:dt:tf);

%%% Setting Initial Conditions
IC_BCR_ref = [X0_EurCI; reshape(stm0,36,1)];

for iter = 1:iterMax_CKF
    
%%% Propagating the Ref State (Traj & STM)
[Times,refStatesSTMs_CKF_EurCI] = ode45(@integrator_FP_nn,time,IC_BCR_ref,options,rad2,u1,u2,aEur,wEur);

% ------------------------------------------------------------------------
%%% Preparing for CKF
% ------------------------------------------------------------------------
%%% Matrices to store results
xHats_CKF = zeros(6,length(Times));
Ps_CKF = zeros(6,6,length(Times));
yks_CKF = zeros(2,length(Times)); % Prefit residual
est_CKF_EurCI = zeros(length(Times),6); % for state estimates
Htildes_CKF = zeros(2,6,length(Times));

%%% Storing Initial Values
Ps_CKF(:,:,1) = PBar0;
Pi = PBar0;
est_CKF_EurCI(1,:) = refStatesSTMs_CKF_EurCI(1,1:6) + xHats_CKF(:,1)';

if data(1,2) == 0 % If first measurement is at t = 0;
    Yk = [data(1,3); data(1,4)]; % First measurement set
    yk = calculatePreFitResiduals_EurCI(refStatesSTMs_CKF_EurCI(1,1:6), states_EJ(1,1:6), states_JSun(1,1:6), states_EaSun(1,1:6), stations_ECI(:,data(1,1),1)', Yk);
    yks_CKF(:,1) = yk;
end

% ------------------------------------------------------------------------
%%% Running CKF
% ------------------------------------------------------------------------
for k = 2:length(Times)
    %%% Grabbing current stm ... phi(tk, tk-1)
    stm = reshape(refStatesSTMs_CKF_EurCI(k,7:end),6,6) * inv(reshape(refStatesSTMs_CKF_EurCI(k-1,7:end),6,6));
        
    if any(data(:,2)==Times(k)) == 1 && measOn == 1 % If there is a measurement (if the data time is in my time vector)
        %%% Find row of measurement
        m = find(data(:,2)==Times(k));
        
        %%% Time Update
        xBari = stm*xHats_CKF(:,k-1);
        PBari = stm*Ps_CKF(:,:,k-1)*stm';
        
        %%%%%%% Compute yi, Hitilde
        %%% Current Measurements
        Yk = [data(m,3); data(m,4)];
        %%% Calculating Htilde
        [ Htilde ] = calculateHtilde_EurCI(refStatesSTMs_CKF_EurCI(k,1:6), states_EJ(k,1:6), states_JSun(k,1:6), states_EaSun(k,1:6), stations_ECI(:,data(m,1),k)');
        Htildes_CKF(:,:,k) = Htilde;
        
        %%% Calculating pre-fit residuals
        [ yk ] = calculatePreFitResiduals_EurCI(refStatesSTMs_CKF_EurCI(k,1:6), states_EJ(k,1:6), states_JSun(k,1:6), states_EaSun(k,1:6), stations_ECI(:,data(m,1),k)', Yk);
        yks_CKF(:,k) = yk;
        
        %%% Kalman Gain Matrix
        Ki = PBari*Htilde'*inv(Htilde*PBari*Htilde' + R);
        
        %%% Measurement Update
        xhat = xBari + Ki*(yk - Htilde*xBari);
        Pi = (eye(6) - Ki*Htilde)*PBari;
        
        %%% Storing Values
        xHats_CKF(:,k) = xhat;
        Ps_CKF(:,:,k) = Pi + gamma*Q*gamma';
        
    else % If no measurement        
        %%% Propagate P and xHat with STM
        xHats_CKF(:,k) = stm*xHats_CKF(:,k-1);
        Ps_CKF(:,:,k) = stm*Ps_CKF(:,:,k-1)*stm' + gamma*Q*gamma';
    end
    
    %%% Update estimate
    est_CKF_EurCI(k,:) = refStatesSTMs_CKF_EurCI(k,1:6) + xHats_CKF(:,k)';
end

if iter < iterMax_CKF
    %%%%% Prepare for next iteration
    %%% Propagating back to t0 w/ stm
    dev0_iter = inv(reshape(refStatesSTMs_CKF_EurCI(end,7:end),6,6)) * xHats_CKF(:,end); % Mapping deviation estimate back to time 0
    X0_iter = refStatesSTMs_CKF_EurCI(1,1:6) + dev0_iter';
    
    IC_BCR_ref = [X0_iter'; reshape(stm0,36,1)];
end

end % iter
%%% Calculate Post-Fit Residuals
[ eks_CKF ] = calcPostFitResiduals(length(Times), yks_CKF, Htildes_CKF, xHats_CKF);


%%% Computing estimate error
error_CKF =  est_CKF_EurCI - states_EurCI_truth;

%%% Compute Altitudes
europaAltitudes_CKF = zeros(length(Times),1); % km
maxCovs = zeros(length(Times),1);
maxSigs = zeros(length(Times),1);
for k = 1:length(Times)
    europaAltitudes_CKF(k) = norm(est_CKF_EurCI(k,1:3)) - rad2; % km
%     maxCovs(k) = norm([Ps_CKF(1,1,1097), Ps_CKF(2,2,1097), Ps_CKF(3,3,1097)]);
    maxCovs(k) = norm([Ps_CKF(1,1,k), Ps_CKF(2,2,k), Ps_CKF(3,3,k)]);
    maxSigs(k) = 3*sqrt(maxCovs(k));
end

% figure; hold all
% plot3(est_CKF_EurCI(:,4),est_CKF_EurCI(:,5),est_CKF_EurCI(:,6),'m','linewidth',2)
% plot3(states_EurCI_truth(:,4),states_EurCI_truth(:,5),states_EurCI_truth(:,6),'b','linewidth',2)
% PlotBoi3('X','Y','Z',16)
% title('Estimated EurCI Velocity (CKF)')
% axis equal
% figure; hold all
% plot3(est_CKF_EurCI(:,1),est_CKF_EurCI(:,2),est_CKF_EurCI(:,3),'m','linewidth',2)
% plot3(states_EurCI_truth(:,1),states_EurCI_truth(:,2),states_EurCI_truth(:,3),'b','linewidth',2)
% PlotBoi3('X','Y','Z',16)
% title('Estimated EurCI Position (CKF)')
% axis equal
%%% Plotting 3-Sig Covariances (to end of measurements)
lw = 1;
if round(Times(end)) ~= Times(end) % If there was a collision, don't plot last error b/c it's ugly
    pEnd = length(Times)-1;
else
    pEnd = length(Times);
end
figure
subplot(3,2,1); hold all; %ylim([-.06 .06])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_CKF(1,1,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_CKF(1,1,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error_CKF(1:pEnd,1),'k')
PlotBoi2('','X Error, km',14)
subplot(3,2,3); hold all; %ylim([-.08 .08])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_CKF(2,2,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_CKF(2,2,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error_CKF(1:pEnd,2),'k')
PlotBoi2('','Y Error, km',14)
subplot(3,2,5); hold all; %ylim([-.1 .1])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_CKF(3,3,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_CKF(3,3,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error_CKF(1:pEnd,3),'k')
PlotBoi2('Time, hr','Z Error, km',14)
subplot(3,2,2); hold all; %ylim([-5e-5 5e-5])
p1 = plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_CKF(4,4,1:pEnd))),'--r','linewidth',lw);
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_CKF(4,4,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error_CKF(1:pEnd,4),'k')
legend([p1],'3\sigma Covariance')
PlotBoi2('','dX Error, km/s',14)
subplot(3,2,4); hold all; %ylim([-5e-5 5e-5])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_CKF(5,5,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_CKF(5,5,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error_CKF(1:pEnd,5),'k')
PlotBoi2('','dY Error, km/s',14)
subplot(3,2,6); hold all; %ylim([-10e-5 10e-5])
plot(Times(1:pEnd)./3600,3*sqrt(squeeze(Ps_CKF(6,6,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,-3*sqrt(squeeze(Ps_CKF(6,6,1:pEnd))),'--r','linewidth',lw)
plot(Times(1:pEnd)./3600,error_CKF(1:pEnd,6),'k')
PlotBoi2('Time, hr','dZ Error, km/s',14)

%%% Plotting Pre-Fit Residuals
PlotPreFitResiduals(Times(1:pEnd), yks_CKF(:,1:pEnd), data(:,1:2), pError, dpError, 16)

%%% Plotting Post-Fit Residuals
PlotPostFitResiduals(Times(1:pEnd), eks_CKF(:,1:pEnd), data(:,1:2), pError, dpError, 16)


% %%% RMS Error Analysis
% sumPos_CKF = 0;
% sumVel_CKF = 0;
% m = size(data,1)*2;
% for k = 1:length(Times)
%     sumPos_CKF = sumPos_CKF + norm(error_CKF(k,1:3));
%     sumVel_CKF = sumVel_CKF + norm(error_CKF(k,4:6));
% end
% 
% RMSPos_CKF = (sumPos_CKF/m)^(1/2)
% RMSVel_CKF = (sumVel_CKF/m)^(1/2)

end % runCKF

if plotAltitudes == 1
    figure; hold all
    plot(Times./3600,europaAltitudes_CKF,'m','linewidth',1.5)
    
    p1 = plot(Times./3600,europaAltitudes_CKF - maxSigs, '--r');
    PlotBoi2('Times,hr','Europa Altitude, km',16)
    legend([p1], '3\sigma Lower Bound')
end
% ------------------------------------------------------------------------
%%% Comparing SRIF and CKF results
% ------------------------------------------------------------------------
% figure
% subplot(3,2,1); hold all; %ylim([-.06 .06])
% plot(Times(6:end)./3600,est_CKF_EurCI(6:end,1)-est_SRIF(6:end,1))
% PlotBoi2('','X Error, km',14)
% subplot(3,2,3); hold all; %ylim([-.08 .08])
% plot(Times(6:end)./3600,est_CKF_EurCI(6:end,2)-est_SRIF(6:end,2))
% PlotBoi2('','Y Error, km',14)
% subplot(3,2,5); hold all; %ylim([-.1 .1])
% plot(Times(6:end)./3600,est_CKF_EurCI(6:end,3)-est_SRIF(6:end,3))
% PlotBoi2('Time, hr','Z Error, km',14)
% subplot(3,2,2); hold all; %ylim([-5e-5 5e-5])
% plot(Times(6:end)./3600,est_CKF_EurCI(6:end,4)-est_SRIF(6:end,4))
% PlotBoi2('','dX Error, km/s',14)
% subplot(3,2,4); hold all; %ylim([-5e-5 5e-5])
% plot(Times(6:end)./3600,est_CKF_EurCI(6:end,5)-est_SRIF(6:end,5))
% PlotBoi2('','dY Error, km/s',14)
% subplot(3,2,6); hold all; %ylim([-10e-5 10e-5])
% plot(Times(6:end)./3600,est_CKF_EurCI(6:end,6)-est_SRIF(6:end,6))
% PlotBoi2('Time, hr','dZ Error, km/s',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Testing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Theoretical Ranges and Range Rates
% ranges = zeros(length(Times),3);
% rangeRates = zeros(length(Times),3);
% for k = 1:length(Times)
%     %%% Find row of measurement
%     m = find(data(:,2)==Times(k));
%     if m ~= 0
%         ranges(k,1) = sqrt(dot(refStatesSTMS_EurCI(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI(1:3,data(m,1),k)',refStatesSTMS_EurCI(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI(1:3,data(m,1),k)'));
%         rangeRates(k,1) = dot(refStatesSTMS_EurCI(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI(1:3,data(m,1),k)',refStatesSTMS_EurCI(k,4:6) + europa_JCI(k,4:6) + jupiter_HC(k,4:6) - earth_HC(k,4:6) - stations_ECI(4:6,data(m,1),k)')/ranges(k,1);
%     else
%         ranges(k,:) = NaN;
%         rangeRates(k,:) = NaN;
%     end
% end
% rangeData = zeros(length(Times),6,4);
% rangeData(:,:,1) = refStatesSTMS_EurCI(:,1:6);
% rangeData(:,:,2) = europa_JCI(:,1:6);
% rangeData(:,:,3) = jupiter_HC(:,1:6);
% rangeData(:,:,4) = earth_HC(:,1:6);
% 
% plot6StateError(Times,rangeData(:,:,1), rangeData_truth(:,:,1))
% PlotBoi2('Times','EurCI Error',16)
% plot6StateError(Times,rangeData(:,:,2), rangeData_truth(:,:,2))
% PlotBoi2('Times','europa_JCI Error',16)
% plot6StateError(Times,rangeData(:,:,3), rangeData_truth(:,:,3))
% PlotBoi2('Times','jupiter_HC Error',16)
% plot6StateError(Times,rangeData(:,:,4), rangeData_truth(:,:,4))
% PlotBoi2('Times','earth_HC Error',16)
% plot6StateError(Times,squeeze(stations_ECI(:,1,:))', squeeze(stations_ECI_truth(:,1,:))')
% PlotBoi2('Times','Station 1 Error',16)
% plot6StateError(Times,squeeze(stations_ECI(:,2,:))', squeeze(stations_ECI_truth(:,2,:))')
% PlotBoi2('Times','Station 2 Error',16)
% plot6StateError(Times,squeeze(stations_ECI(:,3,:))', squeeze(stations_ECI_truth(:,3,:))')
% PlotBoi2('Times','Station 3 Error',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^^^^Testing^^^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


