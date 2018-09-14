clear
clc
% close all
addpath('../../bin')
addpath('../../bin/FilterPlots')
addpath('../../bin/FilterFuncs')
tic

% ------------------------------------------------------------------------
%%% Determining Equations of Motion and Measurement Model
% ------------------------------------------------------------------------
syms xsE ysE zsE dxsE dysE dzsE xsJ ysJ zsJ uEur uJ wEur real
syms xsER ysER zsER dxsER dysER dzsER real
syms xEJ yEJ zEJ dxEJ dyEJ dzEJ real
syms xJSun yJSun zJSun dxJSun dyJSun dzJSun real
syms xEaSun yEaSun zEaSun dxEaSun dyEaSun dzEaSun real
syms xstEa ystEa zstEa dxstEa dystEa dzstEa real

wEur = [0;0;wEur];
rsE = [xsE; ysE; zsE]; % Sat position in Europa-Inertial, km
rsJ = [xsJ; ysJ; zsJ]; % Sat position in Jupiter-Inertial, km
rEJ = rsJ-rsE; % Europa position in Jupiter-Inertial

a2B = -uEur*rsE/(norm(rsE)^3); % km/s^2
a3B = -uJ*(rsJ/(norm(rsJ)^3) - (rEJ)/(norm(rEJ)^3)); % km/s^2

%%% Defining Equations of Motion, State Variables, and A Matrix
EQM = [dxsE; dysE; dzsE; a2B(1) + a3B(1); a2B(2) + a3B(2); a2B(3) + a3B(3)];
state = [xsE; ysE; zsE; dxsE; dysE; dzsE];
Asym = jacobian(EQM, state);

%%% Measurement Model
rEJ = [xEJ; yEJ; zEJ]; vEJ = [dxEJ; dyEJ; dzEJ];
vsE = [dxsE; dysE; dzsE] + cross(wEur,rsE);

rJSun = [xJSun;yJSun;zJSun]; vJSun = [dxJSun;dyJSun;dzJSun];
rEaSun = [xEaSun;yEaSun;zEaSun]; vEaSun = [dxEaSun;dyEaSun;dzEaSun];
rstEa = [xstEa;ystEa;zstEa]; vstEa = [dxstEa;dystEa;dzstEa];

% ECI Range and Range-Rate from Europa-Cenetered-Inertial States
range = sqrt(dot((rsE + rEJ + rJSun - rEaSun - rstEa),(rsE + rEJ + rJSun - rEaSun - rstEa)));
rangeRate = dot((rsE + rEJ + rJSun - rEaSun - rstEa),(vsE + vEJ + vJSun - vEaSun - vstEa))/range;

measFunction = [range; rangeRate];
HtildeEqns = jacobian(measFunction,state);

clear
fprintf('This\n')% how do the htilde plots look when dev = 0?
% how do they look for Proj1?
% ------------------------------------------------------------------------
%%% Plot/Run Options (0 = no/off, 1 = yes/on)
% ------------------------------------------------------------------------
devOn = 1;
measOn = 1; % use measurements?
linearityTest = 1;

plotEst_EurCI                   = 1;
plotEst_EurCEF                  = 0;
plotAltitude                    = 0;
plotResidualsAndCovariance_SRIF = 1;

runCKF = 1;
plotResidualsAndCovariance_CKF = 1;
% ------------------------------------------------------------------------
%%% System Parameters / Givens
% ------------------------------------------------------------------------
AU = 149597870; % km
d2r = pi/180;
day = 86400; % sec

%%% Europa
aEur = 671100; % semimajor axis (distance-normalizing value)
wEur = [0,0,2.047200349303344e-05]; % rad/s
w_n = [0, 0, wEur(3)/wEur(3)]; % angular velocity

%%% Normalizing Factors (real / xN = normalized)
wNorm = 1/wEur(3); % Angular Velocity (rad/s)
rNorm = aEur; % Position (km)
vNorm = rNorm/wNorm; % Velocity (km/s)
tNorm = wNorm;

%%% Gravitational Parameters
uJ = 126672520; % Jupiter, km^3 / s^2
uEur = 3203.413216; % Europa, km^3 / s^2
u = uEur/uJ;

%%% Body Radii
rad1_n = 71492/rNorm; % radius of Jupiter (km normalized)
rad2_n = 1560.8/rNorm; % radius of Europa (km normalized)
rEur = rad2_n.*rNorm;

%%% Rotating frame cooridinates
rB1_BCR_n = [-u, 0, 0];
rB2_BCR_n = [1-u, 0, 0];

%%% Initial Julian Date
[JD0] = JD(2025,3,5,0,0);
fprintf('***JD0 from buffinton paper\n') 

%%% Jupiter
rJ = rad1_n*rNorm; % km

%%% Earth
rE = 6378.1363; % km
wDay = [0,0,2*pi/86400]; % rad/s
fprintf('***Planetary Prms from Valado\n')
fprintf('***Assuming coplanar\n')

%%% Europa period (days)
pEur = 3.551181; 
%%% Setting normalized time vector
ti = 0;
dt = 60;
tf = 1*day;
time = (ti:dt:tf)./tNorm;

% ------------------------------------------------------------------------
%%% Measurements
% ------------------------------------------------------------------------
%%% Loading Truth Traj
load('truthTraj_EurCI.mat');

%%% Loading Measurements
% Time since Epoch, DSS34 Range (km), DSS65 Range, DSS13 Range, DSS34 Range-Rate (km/sec), DSS65 Range-Rate, DSS13 Range-Rate 
dataSet = csvread('finalProjectData_ECI.txt');

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

% ------------------------------------------------------------------------
%%% Station Positions and Velocities
% ------------------------------------------------------------------------
Times = [ti:dt:tf]; % sec

%%% Initial lat/lon/alts of Stations
latLonAlt = [-35.398333, 148.981944, 0.691750;
    40.427222, -4.250556, 0.834539;
    35.247164, 243.205, 1.07114904]; % deg | deg | km

%%% Calculating Intial Positions and Velocities of Stations in ECI
stations_ECI = zeros(6,3,length(Times));
for k = 1:3
    %%% Acquiring ECEF/ECI station position
    [rStat_ECEF] = latlon2surfECEF(latLonAlt(k,1), latLonAlt(k,2), latLonAlt(k,3) + rE); % km

    %%% Storing initial station positions and velocities
    stations_ECI(1:3,k,1) = rStat_ECEF'; % km
    stations_ECI(4:6,k,1) = cross(wDay',rStat_ECEF'); % km/s
end; clear rStat_ECEF

%%% Propagating Station States in ECI
for k = 2:length(Times)
    %%% Calculating current rotation angle
    th = wDay(3)*Times(k); % rad 
    for s = 1:3
        stations_ECI(1:3,s,k) = R3(stations_ECI(1:3,s,1),th); % km
        stations_ECI(4:6,s,k) = R3(stations_ECI(4:6,s,1),th); % km/s
    end
end

% ------------------------------------------------------------------------
%%% Setting Initial Satellite State                               
% ------------------------------------------------------------------------
%%% Starting w/ same state as data creation (truth)
load('OrbitalInfo.mat')
[rH0_ECR_n, vH0_BCR_n] = OE2ECI(a, e, i, raan, w, ta, uEur);
pSat = 2*pi*sqrt((a^3)/uEur); % sec
rH0_ECR_n = (rH0_ECR_n./rNorm)';
v0_BCR_n = (vH0_BCR_n./vNorm)';
% Position wrt Barycenter (normalized)
r0_BCR_n = rH0_ECR_n + rB2_BCR_n;
clear a e i raan w ta rH0_ECR_n vH0_BCR_n

%%% Deviation from truth
devR = [0.001; 0.001; 0.001]./rNorm;
devV = [1e-7; 1e-7; 1e-7]./vNorm;

stm0 = eye(6);
if devOn == 1
    X0_BCR_n = [r0_BCR_n' + devR; v0_BCR_n' + devV];
else
    X0_BCR_n = [r0_BCR_n'; v0_BCR_n'];
end
IC_BCR_n_ref = [X0_BCR_n; reshape(stm0,36,1)];

% ------------------------------------------------------------------------
%%% Preparing for Filter
% ------------------------------------------------------------------------
%%% Measurement noise / uncertainty
pError = (5e-3)*2; % km
fprintf('***Note: pError *2\n')
dpError = 5e-7; % km/s
R = [pError^2 0; 0 dpError^2]/100;

r_sig = 2; % km
v_sig = 0.001; % km/s
 
PBar0 = blkdiag(eye(3)*(r_sig^2),eye(3)*(v_sig^2));
Vk = chol(R);

%%% Setting integrator options
tol = 1E-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% ------------------------------------------------------------------------
%%% Testing Linearity
% ------------------------------------------------------------------------
if linearityTest == 1
IC_BCR_n_ref_t = IC_BCR_n_ref;
[Times_n,x1] = ode45(@integrator_FP,time,IC_BCR_n_ref_t,options,u,rB1_BCR_n,rB2_BCR_n,rad2_n,tNorm,rNorm,vNorm,JD0,uJ,uEur);

x1(:,1:3) = x1(:,1:3).*rNorm;
x1(:,4:6) = x1(:,4:6).*vNorm;

IC_BCR_n_ref_t(1:6) = IC_BCR_n_ref_t(1:6) + [devR; devV];
[Times_n,x2] = ode45(@integrator_FP,time,IC_BCR_n_ref_t,options,u,rB1_BCR_n,rB2_BCR_n,rad2_n,tNorm,rNorm,vNorm,JD0,uJ,uEur);

x2(:,1:3) = x2(:,1:3).*rNorm;
x2(:,4:6) = x2(:,4:6).*vNorm;

dx0 = x2(1,1:6) - x1(1,1:6);

x2Hat = zeros(6,length(Times_n));
x2Hat(:,1) = dx0';
for k = 2:length(Times_n)
    % Phi(tk, tk-1)
    stm = reshape(x1(k,7:end),6,6) * inv(reshape(x1(k-1,7:end),6,6));
    
    x2Hat(:,k) = stm * (x2Hat(:,k-1));
end

x22 = x1(:,1:6) + x2Hat';

plot6StateError(Times_n.*tNorm/3600,x2(:,1:6),x22)
subplot(3,2,1); title('Linearity Test')
plot6StateError(Times_n.*tNorm/3600,x2(:,1:6),x1(:,1:6))
subplot(3,2,1); title('Linearity Test 2')
plot6StateError(Times_n.*tNorm/3600,x1(:,1:6),x22)
subplot(3,2,1); title('Linearity Test 3')

% plot6StateError(Times_n.*tNorm/3600,x2(:,1:6) - x1(:,1:6),x2Hat(1:6,:)')
% subplot(3,2,1); title('2nd Linearity Test')

figure; hold all
plot3(x2(:,1),x2(:,2),x2(:,3),'r','linewidth',2)
plot3(x22(:,1),x22(:,2),x22(:,3),'b','linewidth',2) % plot this stuff w/ i = 0
plot3(x1(:,1),x1(:,2),x1(:,3),'m','linewidth',2) % plot this stuff w/ i = 0
legend('Prop r','STM r'); PlotBoi3('X','Y','Z',16); axis equal

figure; hold all
plot3(x2(:,4),x2(:,5),x2(:,6),'r','linewidth',2)
plot3(x22(:,4),x22(:,5),x22(:,6),'b','linewidth',2) % plot this stuff w/ i = 0
plot3(x1(:,4),x1(:,5),x1(:,6),'m','linewidth',2) % plot this stuff w/ i = 0
legend('Prop v','\Phi v'); PlotBoi3('X','Y','Z',16); axis equal
end

989
return
% ------------------------------------------------------------------------
%%% Reference Trajectory
% ------------------------------------------------------------------------
%%% Creating Earth and Jupiter Heliocentric Inertial States
bc_JCI = zeros(length(Times),6);
europa_BCI = zeros(length(Times),6);
europa_JCI = zeros(length(Times),6); 
earth_HC = zeros(length(Times),6);
jupiter_HC = zeros(length(Times),6);
for k = 1:length(Times)
    th = wEur(3)*Times(k); % rad 
    
    %%% Europa in BCI
    europa_BCI(k,1:3) = R3([1-u; 0; 0],th).*rNorm; % km
    europa_BCI(k,4:6) = cross(wEur,europa_BCI(k,1:3)); % km/s
    
    %%% Barycenter in Jupiter-Centered Inertial (JCI)
    bc_JCI(k,1:3) = R3([u; 0; 0],th).*rNorm; % km
    bc_JCI(k,4:6) = cross(wEur,bc_JCI(k,1:3)); % km/s
    
    %%% Europa in JCI
    europa_JCI(k,1:3) = bc_JCI(k,1:3) + europa_BCI(k,1:3); % km
    europa_JCI(k,4:6) = cross(wEur,europa_JCI(k,1:3)); % km/s
    
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,3,'EME2000'); % km, km/s, km^3/s^2
    earth_HC(k,1:6) = [r',v'];
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,5,'EME2000'); % km, km/s, km^3/s^2
    jupiter_HC(k,1:6) = [r',v'];
end

iterMax_SRIF = 1;
for iter = 1:iterMax_SRIF
%     if iter < iterMax
%         time = (ti:dt:tf/24)./tNorm;
%     elseif iter == iterMax
%         time = (ti:dt:tf)./tNorm;
%     end
%%% Propagating the Ref State (Traj & STM)
[Times_n,refStatesSTMs_BCR_n] = ode45(@integrarefStatesSTMs_BCR_ntor_FP,time,IC_BCR_n_ref,options,u,rB1_BCR_n,rB2_BCR_n,rad2_n,tNorm,rNorm,vNorm,JD0,uJ,uEur);

%%% Un-Normalizing Time
Times = Times_n.*tNorm; % sec
if tf == Times(end)
    Times = [ti:dt:tf];
else
    fprintf('*** COLLISION: T = %5.2f sec\n',Times(end))
    return
end

% ------------------------------------------------------------------------
%%% Converting Reference Trajectory from BCR_n to EurCI and Getting Necessary
%%% Vectors
% ------------------------------------------------------------------------
refStates_BCI = zeros(length(Times),6);
states_EurCR = zeros(length(Times),6);
states_EurCI2 = zeros(length(Times),6);

%%% Barycenter Inertial State (BCI)
for k = 1:length(Times)
    %%% Rotation Angle
    th = Times(k) * wEur(3); % rad
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing EurCI
    %%% States in EurCR
    states_EurCR(k,1:3) = (refStatesSTMs_BCR_n(k,1:3) - rB2_BCR_n).*rNorm; % km/s
    states_EurCR(k,4:6) = refStatesSTMs_BCR_n(k,4:6).*vNorm; % km/s
    
    %%% EurCI 2!!!
    states_EurCI2(k,1:3) = R3(states_EurCR(k,1:3),th);
    states_EurCI2(k,4:6) = R3(states_EurCR(k,4:6),th) + cross(wEur,states_EurCI2(k,1:3));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% States in BCI
    refStates_BCI(k,1:3) = R3(refStatesSTMs_BCR_n(k,1:3),th).*rNorm; % km
    refStates_BCI(k,4:6) = R3(refStatesSTMs_BCR_n(k,4:6),th).*vNorm + cross(wEur,refStates_BCI(k,1:3)); % km/s
    
end

%%% Creating EurCI States
refStatesSTMS_EurCI = zeros(length(Times),6+6^2);
for k = 1:length(Times)
    refStatesSTMS_EurCI(k,1:3) = refStates_BCI(k,1:3) - europa_BCI(k,1:3); % km
    refStatesSTMS_EurCI(k,4:6) = refStates_BCI(k,4:6) - europa_BCI(k,4:6); % km/s
end

%%% Adding STM (calculated from ECI)
refStatesSTMS_EurCI(:,7:end) = refStatesSTMs_BCR_n(:,7:end);

% ------------------------------------------------------------------------
%%% SRIF
% ------------------------------------------------------------------------
%%% Matrices to store results
yks_SRIF = zeros(2,length(Times)); % Prefit residual
Htildes_SRIF = zeros(2,6,length(Times));
bs_SRIF = zeros(6,length(Times));
Rs_SRIF = zeros(6,6,length(Times));
Ps_SRIF = zeros(6,6,length(Times));
xHats_SRIF = zeros(6,length(Times));
est_SRIF_EurCI = zeros(length(Times),6); % for state estimates

%%% Setting initial values
Ps_SRIF(:,:,1) = PBar0;
Rs_SRIF(:,:,1) = chol(Ps_SRIF(:,:,1)); 
est_SRIF_EurCI(1,:) = refStatesSTMS_EurCI(1,1:6) + xHats_SRIF(:,1)';

if data(1,2) == 0 % If first measurement is at t = 0;
    Yk = [data(1,3); data(1,4)]; % First measurement set
    yk = calculatePreFitResiduals_EurCI(refStatesSTMS_EurCI(1,1:6), europa_JCI(1,1:6), jupiter_HC(1,1:6), earth_HC(1,1:6), stations_ECI(:,data(1,1),k)', Yk); 
    yks_SRIF(:,1) = Vk\yk;
end

kEnd = length(Times);
for k = 2:kEnd
    %%% Grabbing current stm ... phi(tk, tk-1)
    stm = reshape(refStatesSTMS_EurCI(k,7:end),6,6) * inv(reshape(refStatesSTMS_EurCI(k-1,7:end),6,6));

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
        tempEurCI = zeros(1,6);
        tempEurCI(1:3) = refStatesSTMS_EurCI(k,1:3);
        tempEurCI(4:6) = refStatesSTMS_EurCI(k,4:6) - cross(wEur,refStatesSTMS_EurCI(k,1:3));
        [ Htilde ] = calculateHtilde_EurCI(tempEurCI(1,1:6), europa_JCI(k,1:6), jupiter_HC(k,1:6), earth_HC(k,1:6), stations_ECI(:,data(m,1),k)', wEur(3));
%         [ Htilde ] = calculateHtilde_EurCI(refStatesSTMS_EurCI(k,1:6), europa_JCI(k,1:6), jupiter_HC(k,1:6), earth_HC(k,1:6), stations_ECI(:,data(m,1),k)', wEur(3));
        
        %%% Calculating pre-fit residuals
        [ yk ] = calculatePreFitResiduals_EurCI(refStatesSTMS_EurCI(k,1:6), europa_JCI(k,1:6), jupiter_HC(k,1:6), earth_HC(k,1:6), stations_ECI(:,data(m,1),k)', Yk);

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
    est_SRIF_EurCI(k,:) = refStatesSTMS_EurCI(k,1:6) + xHats_SRIF(:,k)';
    
end

if iter < iterMax_SRIF
%     figure
%     plot3(est_SRIF_EurCI(:,1),est_SRIF_EurCI(:,2),est_SRIF_EurCI(:,3))
%     plot6StateError(Times./3600,est_SRIF_EurCI(:,1:6),States_EurCI_truth(1:length(Times),1:6))
    
    %%%%% Prepare for next iteration
    %%% Propagating back to t0 w/ stm
    dev0_iter = inv(reshape(refStatesSTMS_EurCI(end,7:end),6,6)) * xHats_SRIF(:,end); % Mapping deviation estimate back to time 0
    X0_iter = refStatesSTMS_EurCI(1,1:6) + dev0_iter';
    
    %%% EurCI States
    r0_EurCI_iter = X0_iter(1:3);
    v0_EurCI_iter = X0_iter(4:6);
    
    %%% BCR States (alligned w/ BCI at t0)
    r0_BCR_n_iter = ((europa_BCI(1,1:3) + r0_EurCI_iter)./rNorm)';
%     v0_BCR_n_iter = ((europa_BCI(1,4:6) + v0_EurCI_iter)./vNorm)';
    v0_BCR_n_iter = ((v0_EurCI_iter - cross(wEur,r0_EurCI_iter)./vNorm))';
    
    X0_BCR_n_iter = [r0_BCR_n_iter; v0_BCR_n_iter];
    IC_BCR_n_ref = [X0_BCR_n_iter; reshape(stm0,36,1)];
    
end
end % iterations

plot6StateError(Times./3600, est_SRIF_EurCI(:,1:6), States_EurCI_truth(:,1:6))
title('est, truth (SRIF)')
plot6StateError(Times./3600, refStatesSTMS_EurCI(:,1:6), States_EurCI_truth(:,1:6))
title('ref, truth (SRIF)')
% ------------------------------------------------------------------------
%%% Post-SRIF
% ------------------------------------------------------------------------
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

%%% Plotting Estimated EurCI States
if plotEst_EurCI == 1
    figure; hold all
    plot3(est_SRIF_EurCI(:,1),est_SRIF_EurCI(:,2),est_SRIF_EurCI(:,3),'m','linewidth',2)
    bodySurface3(rad2_n.*rNorm, [0,0,0], [0,1,1])
    PlotBoi3('X','Y','Z',16)
    title('Estimated EurCI')
    axis equal
end

%%% Plotting Estimated EurCEF States
est_SRIF_EurCEF = zeros(size(est_SRIF_EurCI));
for k = 1:length(Times)
    th = wEur(3)*Times(k);
    est_SRIF_EurCEF(k,1:3) = R3(est_SRIF_EurCI(k,1:3),-th);
    est_SRIF_EurCEF(k,4:6) = R3(est_SRIF_EurCI(k,4:6),-th);
end
if plotEst_EurCEF == 1
    figure; hold all
    plot3(est_SRIF_EurCEF(:,1),est_SRIF_EurCEF(:,2),est_SRIF_EurCEF(:,3),'m','linewidth',2)
    bodySurface3(rad2_n.*rNorm, [0,0,0], [0,1,1])
    PlotBoi3('X','Y','Z',16)
    title('Estimated EurCEF')
    axis equal
end

%%% Plotting Altitude
europaAltitudes = zeros(length(Times),1);
for k = 1:length(Times)
    europaAltitudes(k) = norm(est_SRIF_EurCI(k,1:3)) - rad2_n*rNorm;
end
if plotAltitude == 1
    figure
    plot(Times./3600,europaAltitudes,'linewidth',1.5)
    PlotBoi2('Times, hr','Europa Altitude, km',16)
end

%%% Creating Error from Truth
error = est_SRIF_EurCI - States_EurCI_truth;

%%% Plotting 3-Sig Covariances (to end of measurements)
if plotResidualsAndCovariance_SRIF == 1
    lw = 1.5;
    figure
    subplot(3,2,1); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_SRIF(1,1,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_SRIF(1,1,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error(:,1),'k')
    PlotBoi2('','X Error, km',14)
    subplot(3,2,3); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_SRIF(2,2,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_SRIF(2,2,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error(:,2),'k')
    PlotBoi2('','Y Error, km',14)
    subplot(3,2,5); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_SRIF(3,3,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_SRIF(3,3,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error(:,3),'k')
    PlotBoi2('Time, hr','Z Error, km',14)
    subplot(3,2,2); hold all;
    p1 = plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_SRIF(4,4,1:end))),'--r','linewidth',lw);
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_SRIF(4,4,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error(:,4),'k')
    legend([p1],'3\sigma Covariance')
    PlotBoi2('','dX Error, km/s',14)
    subplot(3,2,4); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_SRIF(5,5,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_SRIF(5,5,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error(:,5),'k')
    PlotBoi2('','dY Error, km/s',14)
    subplot(3,2,6); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_SRIF(6,6,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_SRIF(6,6,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error(:,6),'k')
    PlotBoi2('Time, hr','dZ Error, km/s',14)

    %%% Plotting Pre-Fit Residuals
    PlotPreFitResiduals(Times(1:end), yks_SRIF(:,1:end), data(:,1:2), pError, dpError, 16)

    %%% Plotting Post-Fit Residuals
    PlotPostFitResiduals(Times(1:end), eks_SRIF(:,1:end), data(:,1:2), pError, dpError, 16)
end
figure; hold all
plot3(est_SRIF_EurCI(:,4),est_SRIF_EurCI(:,5),est_SRIF_EurCI(:,6),'m','linewidth',2)
plot3(States_EurCI_truth(:,4),States_EurCI_truth(:,5),States_EurCI_truth(:,6),'b','linewidth',2)
PlotBoi3('X','Y','Z',16)
title('Estimated EurCI Velocity (SRIF)')
axis equal
figure; hold all
plot3(est_SRIF_EurCI(:,1),est_SRIF_EurCI(:,2),est_SRIF_EurCI(:,3),'m','linewidth',2)
plot3(States_EurCI_truth(:,1),States_EurCI_truth(:,2),States_EurCI_truth(:,3),'b','linewidth',2)
PlotBoi3('X','Y','Z',16)
title('Estimated EurCI Position (SRIF)')
axis equal








% ========================================================================
%%% CKF (p203)
% ========================================================================
% 989
% return
if runCKF == 1
%%% Resetting normalized time vector
ti = 0;
dt = 60;
tf = 1*day;
time = (ti:dt:tf)./tNorm;

%%% Setting Initial Conditions
IC_BCR_n_ref = [X0_BCR_n; reshape(stm0,36,1)];

iterMax_CKF = 1;
for iter = 1:iterMax_CKF
%     if iter < iterMax_CKF
%         time = (ti:dt:tf/24)./tNorm;
%     elseif iter == iterMax_CKF
%         time = (ti:dt:tf)./tNorm;
%     end

%%% Propagating the Ref State (Traj & STM)
[Times_n,refStatesSTMs_CKF_BCR_n] = ode45(@integrator_FP,time,IC_BCR_n_ref,options,u,rB1_BCR_n,rB2_BCR_n,rad2_n,tNorm,rNorm,vNorm,JD0,uJ,uEur);

%%% Un-Normalizing Time
Times = Times_n.*tNorm; % sec
if tf == Times(end)
    Times = [ti:dt:tf];
else
    warning('Weird collision or something\n')
end

% ------------------------------------------------------------------------
%%% Converting BCR_n to EurCI
% ------------------------------------------------------------------------
refStates_CKF_BCI = zeros(length(Times),6);
refStatesSTMS_CKF_EurCI = zeros(length(Times),6+6^2);

%%% Barycenter Inertial State (BCI)
for k = 1:length(Times)
    %%% Rotation Angle
    th = Times(k) * wEur(3); % rad
    
    %%% States in BCI
    refStates_CKF_BCI(k,1:3) = R3(refStatesSTMs_CKF_BCR_n(k,1:3),th).*rNorm; % km
    refStates_CKF_BCI(k,4:6) = R3(refStatesSTMs_CKF_BCR_n(k,4:6),th).*vNorm + cross(wEur,refStates_CKF_BCI(k,1:3)); % km/s

end

%%% Creating EurCI States
for k = 1:length(Times)
    refStatesSTMS_CKF_EurCI(k,1:3) = refStates_CKF_BCI(k,1:3) - europa_BCI(k,1:3); % km
    refStatesSTMS_CKF_EurCI(k,4:6) = refStates_CKF_BCI(k,4:6) - europa_BCI(k,4:6); % km/s
end

%%% Adding STM (calculated from ECI)
refStatesSTMS_CKF_EurCI(:,7:end) = refStatesSTMs_CKF_BCR_n(:,7:end);

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
est_CKF_EurCI(1,:) = refStatesSTMS_CKF_EurCI(1,1:6) + xHats_CKF(:,1)';

if data(1,2) == 0 % If first measurement is at t = 0;
    Yk = [data(1,3); data(1,4)]; % First measurement set
    yk = calculatePreFitResiduals_EurCI(refStatesSTMS_CKF_EurCI(1,1:6), europa_JCI(1,1:6), jupiter_HC(1,1:6), earth_HC(1,1:6), stations_ECI(:,data(1,1),k)', Yk); 
    yks_CKF(:,1) = yk;
end

% ------------------------------------------------------------------------
%%% Running CKF
% ------------------------------------------------------------------------
for k = 2:length(Times)
    %%% Grabbing current stm ... phi(tk, tk-1)
    stm = reshape(refStatesSTMS_CKF_EurCI(k,7:end),6,6) * inv(reshape(refStatesSTMS_CKF_EurCI(k-1,7:end),6,6));
        
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
        tempEurCI = zeros(1,6);
        tempEurCI(1:3) = refStatesSTMS_CKF_EurCI(k,1:3);
        tempEurCI(4:6) = refStatesSTMS_CKF_EurCI(k,4:6) - cross(wEur,refStatesSTMS_CKF_EurCI(k,1:3));
        [ Htilde ] = calculateHtilde_EurCI(tempEurCI(1:6), europa_JCI(k,1:6), jupiter_HC(k,1:6), earth_HC(k,1:6), stations_ECI(:,data(m,1),k)', wEur(3));
%         [ Htilde ] = calculateHtilde_EurCI(refStatesSTMS_CKF_EurCI(k,1:6), europa_JCI(k,1:6), jupiter_HC(k,1:6), earth_HC(k,1:6), stations_ECI(:,data(m,1),k)', wEur(3));
        Htildes_CKF(:,:,k) = Htilde;
        
        %%% Calculating pre-fit residuals
        [ yk ] = calculatePreFitResiduals_EurCI(refStatesSTMS_CKF_EurCI(k,1:6), europa_JCI(k,1:6), jupiter_HC(k,1:6), earth_HC(k,1:6), stations_ECI(:,data(m,1),k)', Yk, wEur(3));
        yks_CKF(:,k) = yk;
        
        %%% Kalman Gain Matrix
        Ki = PBari*Htilde'*inv(Htilde*PBari*Htilde' + R);
        
        %%% Measurement Update
        xhat = xBari + Ki*(yk - Htilde*xBari);
        Pi = (eye(6) - Ki*Htilde)*PBari;
        
        %%% Storing Values
        xHats_CKF(:,k) = xhat;
        Ps_CKF(:,:,k) = Pi;
        
    else % If no measurement        
        %%% Propagate P and xHat with STM
        xHats_CKF(:,k) = stm*xHats_CKF(:,k-1);
        Ps_CKF(:,:,k) = stm*Ps_CKF(:,:,k-1)*stm';
    end
    
    %%% Update estimate
    est_CKF_EurCI(k,:) = refStatesSTMS_CKF_EurCI(k,1:6) + xHats_CKF(:,k)';
end

%%% Calculate Post-Fit Residuals
[ eks_CKF ] = calcPostFitResiduals(length(Times), yks_CKF, Htildes_CKF, xHats_CKF);


%%% Computing estimate error
error_CKF =  est_CKF_EurCI(:,:) - States_EurCI_truth(:,1:6);

if iter < iterMax_CKF
    d
end

end % iter

figure; hold all
plot3(est_CKF_EurCI(:,4),est_CKF_EurCI(:,5),est_CKF_EurCI(:,6),'m','linewidth',2)
plot3(States_EurCI_truth(:,4),States_EurCI_truth(:,5),States_EurCI_truth(:,6),'b','linewidth',2)
PlotBoi3('X','Y','Z',16)
title('Estimated EurCI Velocity (CKF)')
axis equal
figure; hold all
plot3(est_CKF_EurCI(:,1),est_CKF_EurCI(:,2),est_CKF_EurCI(:,3),'m','linewidth',2)
plot3(States_EurCI_truth(:,1),States_EurCI_truth(:,2),States_EurCI_truth(:,3),'b','linewidth',2)
PlotBoi3('X','Y','Z',16)
title('Estimated EurCI Position (CKF)')
axis equal

plot6StateError(Times./3600, est_CKF_EurCI(:,1:6), States_EurCI_truth(:,1:6))
title('est, truth (CKF)')
plot6StateError(Times./3600, refStatesSTMS_CKF_EurCI(:,1:6), States_EurCI_truth(:,1:6))
title('ref, truth (CKF)')

%%% Plotting 3-Sig Covariances (to end of measurements)
if plotResidualsAndCovariance_SRIF == 1
    lw = 1.5;
    figure
    subplot(3,2,1); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_CKF(1,1,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_CKF(1,1,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error_CKF(:,1),'k')
    PlotBoi2('','X Error, km',14)
    subplot(3,2,3); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_CKF(2,2,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_CKF(2,2,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error_CKF(:,2),'k')
    PlotBoi2('','Y Error, km',14)
    subplot(3,2,5); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_CKF(3,3,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_CKF(3,3,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error_CKF(:,3),'k')
    PlotBoi2('Time, hr','Z Error, km',14)
    subplot(3,2,2); hold all;
    p1 = plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_CKF(4,4,1:end))),'--r','linewidth',lw);
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_CKF(4,4,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error_CKF(:,4),'k')
    legend([p1],'3\sigma Covariance')
    PlotBoi2('','dX Error, km/s',14)
    subplot(3,2,4); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_CKF(5,5,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_CKF(5,5,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error_CKF(:,5),'k')
    PlotBoi2('','dY Error, km/s',14)
    subplot(3,2,6); hold all;
    plot(Times(1:end)./3600,3*sqrt(squeeze(Ps_CKF(6,6,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,-3*sqrt(squeeze(Ps_CKF(6,6,1:end))),'--r','linewidth',lw)
    plot(Times(1:end)./3600,error_CKF(:,6),'k')
    PlotBoi2('Time, hr','dZ Error, km/s',14)

    %%% Plotting Pre-Fit Residuals
    PlotPreFitResiduals(Times(1:end), yks_CKF(:,1:end), data(:,1:2), pError, dpError, 16)

    %%% Plotting Post-Fit Residuals
    PlotPostFitResiduals(Times(1:end), eks_CKF(:,1:end), data(:,1:2), pError, dpError, 16)
end

end % runCKF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Testing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% Comparing Range/RangeRate components calculated here vs those
% %%%%% caluculated in CreateData_FP.m
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
Htildes_CKF(Htildes_CKF == 0) = NaN
figure
subplot(6,2,1)
plot(Times./3600,squeeze(((Htildes_CKF(1,1,:)-Htildes_CKF(1,1,2))/Htildes_CKF(1,1,2))).*100)
subplot(6,2,3)
plot(Times./3600,squeeze(((Htildes_CKF(1,2,:)-Htildes_CKF(1,2,2))/Htildes_CKF(1,2,2))).*100)
subplot(6,2,5)
plot(Times./3600,squeeze(((Htildes_CKF(1,3,:)-Htildes_CKF(1,3,2))/Htildes_CKF(1,3,2))).*100)
subplot(6,2,7)
plot(Times./3600,squeeze(((Htildes_CKF(1,4,:)-Htildes_CKF(1,4,2))/Htildes_CKF(1,4,2))).*100)
subplot(6,2,9)
plot(Times./3600,squeeze(((Htildes_CKF(1,5,:)-Htildes_CKF(1,5,2))/Htildes_CKF(1,5,2))).*100)
subplot(6,2,11)
plot(Times./3600,squeeze(((Htildes_CKF(1,6,:)-Htildes_CKF(1,6,2))/Htildes_CKF(1,6,2))).*100)
subplot(6,2,2)
plot(Times./3600,squeeze(((Htildes_CKF(2,1,:)-Htildes_CKF(2,1,2))/Htildes_CKF(2,1,2))).*100)
subplot(6,2,4)
plot(Times./3600,squeeze(((Htildes_CKF(2,2,:)-Htildes_CKF(2,2,2))/Htildes_CKF(2,2,2))).*100)
subplot(6,2,6)
plot(Times./3600,squeeze(((Htildes_CKF(2,3,:)-Htildes_CKF(2,3,2))/Htildes_CKF(2,3,2))).*100)
subplot(6,2,8)
plot(Times./3600,squeeze(((Htildes_CKF(2,4,:)-Htildes_CKF(2,4,2))/Htildes_CKF(2,4,2))).*100)
subplot(6,2,10)
plot(Times./3600,squeeze(((Htildes_CKF(2,5,:)-Htildes_CKF(2,5,2))/Htildes_CKF(2,5,2))).*100)
subplot(6,2,12)
plot(Times./3600,squeeze(((Htildes_CKF(2,6,:)-Htildes_CKF(2,6,2))/Htildes_CKF(2,6,2))).*100)
