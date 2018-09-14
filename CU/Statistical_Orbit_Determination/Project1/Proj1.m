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
% dt = 1000; % (sec)
% tf = 86400*300;

% ------------------------------------------------------------------------
%%% Determining Equations of Motion
% ------------------------------------------------------------------------
syms x y z dx dy dz xs ys zs uE uS Cr Pphi Am c AU

r = [x; y; z]; % km
rs = [xs; ys; zs]; % km
Rs = rs-r;

a2B = -uE*r/(norm(r)^3);
a3B = -uS*((r-rs)/(norm(Rs)^3) + rs/(norm(rs)^3));

%%% SRP equations of motion (acceleration)
aSRP = (-Cr*Pphi*Am*(AU^2)/((norm(Rs)^3)*c))*Rs;

%%% Defining Equations of Motion, State Variables, and A Matrix
EQM = [dx; dy; dz; a2B(1) + a3B(1) + aSRP(1); a2B(2) + a3B(2) + aSRP(2); a2B(3) + a3B(3) + aSRP(3); 0];
state = [x; y; z; dx; dy; dz; Cr];
Asym = jacobian(EQM, state);

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
rE = 6378.1363; % km
rE_SOI = 925000; % km
uE = 398600.4415; % km^3/s^2
uS = 132712440017.987; % km^3/s^2
wE = [0, 0, 7.29211585275553e-5]; % rad/s
AU =  149597870.7; % km
Am = 0.01e-6; % km^2/kg 
Pphi = 1357; % kg/s^3
c = 299792.458; % km/s
JD0 =  2456296.25; % days
Cr = 1.2; 

% ------------------------------------------------------------------------
%%% Station Positions and Velocities
% ------------------------------------------------------------------------
if recalculateStationParameters == 1
    %%% Initial lat/lon/alts of Stations
    latLonAlt = [-35.398333, 148.981944, 0.691750;
        40.427222, 4.250556, 0.834539;
        35.247164, 243.205, 1.07114904]; % deg | deg | km

    %%% Calculating Intial Positions and Velocities of Stations in ECI
    statStates = zeros(6,3,length(Times));
    for k = 1:3
        %%% Acquiring ECEF/ECI station position
        [r] = latlon2surfECEF(latLonAlt(k,1), latLonAlt(k,2), latLonAlt(k,3) + rE); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correct way of handling altitude?

        %%% Storing initial station positions and velocities
        statStates(1:3,k,1) = r'; % km
        statStates(4:6,k,1) = cross(wE',r'); % km/s
    end; clear r

    %%% Propagating Station States in ECI
    for k = 2:length(Times)
        %%% Calculating current rotation angle
        th = wE(3)*Times(k); % rad 
        for s = 1:3
            statStates(1:3,s,k) = R3(statStates(1:3,s,1),th); % km
            statStates(4:6,s,k) = R3(statStates(4:6,s,1),th); % km/s
        end
    end
    
    %%% Saving Station Parameters
    save('statStates_A.mat','statStates');
else
    %%% Loading Station Parameters
    load('statStates_A.mat')
end
warning('For 2b don''t forget to change madrid longitude')

% ------------------------------------------------------------------------
%%% Measurements
% ------------------------------------------------------------------------
%%% Measurement noise / uncertainty
pError = 0.005; % km
dpError = 0.0000005; % km/s
R = [pError^2 0; 0 dpError^2];

%%% Loading Measurements
% Time since Epoch, DSS34 Range (km), DSS65 Range, DSS13 Range, DSS34 Range-Rate (km/sec), DSS65 Range-Rate, DSS13 Range-Rate 
dataSet = csvread('Project1a_Obs.txt',2,0);
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

% ------------------------------------------------------------------------
%%% Setting Initial Satellite State                               
% ------------------------------------------------------------------------
r0_Sat = [-274096790.0; -92859240.0; -40199490.0]; % km
v0_Sat = [32.67; -8.94; -3.88]; % km/s

X0 = [r0_Sat; v0_Sat; Cr]; % (km, km/s)
stm0 = eye(7);
IC_ref = [X0; reshape(stm0,49,1)];

% ------------------------------------------------------------------------
%%% Propagating States
% ------------------------------------------------------------------------
%%% Setting time vector
time = t0:dt:tf; % (sec)

%%% Setting integrator options
tol = 1E-13;
options = odeset('RelTol',tol,'AbsTol',tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rerunRefAndSRIF == 1
for iter = 1:2
%%% Either Calculating Reference States or Loading
    %%% Propagating the Ref State (Traj & STM)
    [Times,refStatesSTMs] = ode45(@proj1Integrator,time,IC_ref,options,uS,uE,AU,Am,Pphi,c,JD0);
    
    %%% Calculating Sun Positions
    rSuns = zeros(length(Times),3);
    for k = 1:length(Times)
        [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,3,'EME2000'); % km, km/s, km^3/s^2
        rSuns(k,:) = -r;
    end

% ------------------------------------------------------------------------
%%% Preparing for Filter
% ------------------------------------------------------------------------
PBar0 = blkdiag(eye(3)*(100^2),eye(3)*(0.1^2),eye(1)*(0.1^2));

% ------------------------------------------------------------------------
%%% SRIF
% ------------------------------------------------------------------------
fprintf('SRIF Started: %3.1f sec\n',toc)
%%% Matrices to store results
yks_SRIF = zeros(2,length(Times)); % Prefit residual
Htildes_SRIF = zeros(2,7,length(Times));
bs_SRIF = zeros(7,length(Times));
Rs_SRIF = zeros(7,7,length(Times));
Ps_SRIF = zeros(7,7,length(Times));
xHats_SRIF = zeros(7,length(Times));
est_SRIF = zeros(length(Times),7); % for state estimates

%%% Setting initial values
Ps_SRIF(:,:,1) = PBar0;
Rs_SRIF(:,:,1) = chol(Ps_SRIF(:,:,1)); 
est_SRIF(1,:) = refStatesSTMs(1,1:7) + xHats_SRIF(:,1)';
Yk = [data(1,3); data(1,4)]; % First measurement set
yk = calculatePreFitResiduals(refStatesSTMs(1,1:6), statStates(:,data(1,1),1)', Yk); 
yks_SRIF(:,1) = yk;

% kEnd = 290000;
kEnd = length(Times);
for k = 2:kEnd
    %%% Grabbing current stm ... phi(tk, tk-1)
    stm = reshape(refStatesSTMs(k,8:end),7,7) * inv(reshape(refStatesSTMs(k-1,8:end),7,7));

    if any(data(:,2)==Times(k)) == 1 % If there is a measurement
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
        [ Htilde ] = calculateHtilde(refStatesSTMs(k,1:6), statStates(:,data(m,1),k)');
        Htilde = [Htilde, zeros(2,1)]; % Adding column of zeros for Cr
        %%% Calculating pre-fit residuals
        [ yk ] = calculatePreFitResiduals(refStatesSTMs(k,1:6), statStates(:,data(m,1),k)', Yk);

        %%%%%%% Whitening measurements
        Vk = chol(R);
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
    est_SRIF(k,:) = refStatesSTMs(k,1:7) + xHats_SRIF(:,k)';
    
end



%%% For iterating, recalculate a priori using estimate at time of last measurement
mEnd = find(Times==data(end,2)); % Last measurement time index

x0_iter = inv(reshape(refStatesSTMs(mEnd,8:end),7,7)) * xHats_SRIF(:,mEnd); % Mapping deviation estimate back to time 0
X0 = refStatesSTMs(1,1:7) + x0_iter';


%%% For iterating, reset a priori
IC_ref = [X0'; reshape(stm0,49,1)];
end % end iter

%%% Saving Ref State for Loading Purposes
save('RefAndSRIFResults_A.mat','kEnd','Times','refStatesSTMs','rSuns','yks_SRIF','Htildes_SRIF','bs_SRIF','Ps_SRIF','xHats_SRIF','est_SRIF','Rs_SRIF','Vk','mEnd');

else % Skip ref and SRIF and just reload them
    load('RefAndSRIFResults_A.mat')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerunRefAndSRIF == 1
    %%% Computing Covariances
    for k = 2:kEnd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not full time
        Ps_SRIF(:,:,k) = inv(Rs_SRIF(:,:,k)) * inv(Rs_SRIF(:,:,k))';
    end

    %%% Calculate Post-Fit Residuals
    [ eks_SRIF ] = calcPostFitResiduals(length(Times), yks_SRIF, Htildes_SRIF, xHats_SRIF);

    %%% Un-Whiten
    for k = 1:kEnd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not full time
        eks_SRIF(:,k) = Vk*eks_SRIF(:,k);
    end
    
    %%% Save Data
    save('PsAndEs_A.mat','Ps_SRIF','eks_SRIF')
else
    
    load('PsAndEs_A.mat')
end

%%% Plotting 3-Sig Covariances (to end of measurements)
lw = 1.5;
figure
subplot(4,2,1); hold all; ylim([-.1 .1])
plot(Times(1:kEnd),3*sqrt(squeeze(Ps_SRIF(1,1,1:kEnd))),'--r','linewidth',lw)
plot(Times(1:kEnd),-3*sqrt(squeeze(Ps_SRIF(1,1,1:kEnd))),'--r','linewidth',lw)
PlotBoi2('','X Error, km',14)
subplot(4,2,3); hold all; ylim([-10 10])
plot(Times(1:kEnd),3*sqrt(squeeze(Ps_SRIF(2,2,1:kEnd))),'--r','linewidth',lw)
plot(Times(1:kEnd),-3*sqrt(squeeze(Ps_SRIF(2,2,1:kEnd))),'--r','linewidth',lw)
PlotBoi2('','Y Error, km',14)
subplot(4,2,5); hold all; ylim([-60 60])
plot(Times(1:kEnd),3*sqrt(squeeze(Ps_SRIF(3,3,1:kEnd))),'--r','linewidth',lw)
plot(Times(1:kEnd),-3*sqrt(squeeze(Ps_SRIF(3,3,1:kEnd))),'--r','linewidth',lw)
PlotBoi2('','Z Error, km',14)
subplot(4,2,2); hold all; ylim([-2e-5 2e-5])
plot(Times(1:kEnd),3*sqrt(squeeze(Ps_SRIF(4,4,1:kEnd))),'--r','linewidth',lw)
plot(Times(1:kEnd),-3*sqrt(squeeze(Ps_SRIF(4,4,1:kEnd))),'--r','linewidth',lw)
PlotBoi2('','dX Error, km',14)
subplot(4,2,4); hold all; ylim([-2e-5 2e-5])
plot(Times(1:kEnd),3*sqrt(squeeze(Ps_SRIF(5,5,1:kEnd))),'--r','linewidth',lw)
plot(Times(1:kEnd),-3*sqrt(squeeze(Ps_SRIF(5,5,1:kEnd))),'--r','linewidth',lw)
PlotBoi2('','dY Error, km',14)
subplot(4,2,6); hold all; ylim([-2e-5 2e-5])
plot(Times(1:kEnd),3*sqrt(squeeze(Ps_SRIF(6,6,1:kEnd))),'--r','linewidth',lw)
plot(Times(1:kEnd),-3*sqrt(squeeze(Ps_SRIF(6,6,1:kEnd))),'--r','linewidth',lw)
PlotBoi2('','dZ Error, km',14)
subplot(4,2,7); hold all; ylim([-.03 .03])
plot(Times(1:kEnd),3*sqrt(squeeze(Ps_SRIF(7,7,1:kEnd))),'--r','linewidth',lw)
plot(Times(1:kEnd),-3*sqrt(squeeze(Ps_SRIF(7,7,1:kEnd))),'--r','linewidth',lw)
PlotBoi2('','Cr Error',14)

%%% Plotting Pre-Fit Residuals
PlotPostFitResiduals(Times(1:kEnd), yks_SRIF(:,1:kEnd), pError, dpError, 16) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not full time

%%% Plotting Post-Fit Residuals
PlotPostFitResiduals(Times(1:kEnd), eks_SRIF(:,1:kEnd), pError, dpError, 16) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not full time

% ------------------------------------------------------------------------
%%% Smoothing
% ------------------------------------------------------------------------
% if smoothingOn == 1
% fprintf('Smoothing Started: %3.1f sec\n',toc)
% %%% Preparing
% estSmooth_SRIF = zeros(size(est_SRIF));
% estSmooth_SRIF(end,:) = est_SRIF(end,:);
% xHats_lk = zeros(size(xHats_SRIF));
% Ps_lk = zeros(size(Ps_SRIF));
% 
% %%% CKF Smoothing
% xHats_lk(end,:) = xHats_CKF(end,:)';
% Ps_lk(:,:,end) = P_est_CKF(:,:,end);
% for l = length(Times):-1:2
%     k = l-1;
%     
%     phi_tkp1_tk = refSTMs_CKF(:,:,k+1) * inv(refSTMs_CKF(:,:,k)); % Phi(tk+1, tk)
% 
%     Sk = P_est_CKF(:,:,k) * phi_tkp1_tk' * inv(P_bar_CKF(:,:,k+1));
%     xHats_lk(k,:) = (xHats_CKF(k,:)' + Sk*(xHats_lk(k+1,:)' - phi_tkp1_tk*xHats_CKF(k,:)'))'; 
%     Ps_lk(:,:,k) = P_est_CKF(:,:,k) + Sk*(Ps_lk(:,:,k+1) - P_bar_CKF(:,:,k+1))*Sk';
%     
%     
% %     phi_tkp1_tk = refSTMs_CKF(:,:,k+1) * inv(refSTMs_CKF(:,:,k)); % Phi(tk+1, tk)
% % 
% %     Sk = P_est_CKF(:,:,k) * phi_tkp1_tk' * inv(P_bar_CKF(:,:,k+1));
% %     xHats_lk(k,:) = (xHats_CKF(k,:)' + Sk*(xHats_lk(k+1,:)' - phi_tkp1_tk*xHats_CKF(k,:)'))'; 
% %     Ps_lk(:,:,k) = P_est_CKF(:,:,k) + Sk*(Ps_lk(:,:,k+1) - P_bar_CKF(:,:,k+1))*Sk';
%     
% 
% 
% %     if processNoiseOn == 1 && noisyRangeMeas(k) ~= 0 % process noise (at tk)
% %         Sk = P_est_CKF(:,:,k) * phi_tkp1_tk' * inv(P_bar_CKF(:,:,k+1));
% %         xHats_lk(k,:) = (xHats_CKF(k,:)' + Sk*(xHats_lk(k+1,:)' - phi_tkp1_tk*xHats_CKF(k,:)'))'; 
% %         Ps_lk(:,:,k) = P_est_CKF(:,:,k) + Sk*(Ps_lk(:,:,k+1) - P_bar_CKF(:,:,k+1))*Sk';
% %     else % No process noise
% % %         phi_tk_tl = refSTMs_CKF(:,:,k) * inv(refSTMs_CKF(:,:,l)); % Phi(tk, tl)
% %         Sk = inv(phi_tkp1_tk);
% %         xHats_lk(k,:) = Sk * xHats_lk(k+1,:)';
% %         Ps_lk(:,:,k) = Sk * Ps_lk(:,:,k+1) * Sk';
% %     end
%     
%     estSmooth_SRIF(k,:) = refStatesSTM_CKF(k,1:6) + xHats_lk(k,:);
% end
% 
% 
% errorSmooth_CKF = truthStateSTM(:,1:6) - estSmooth_SRIF;
% PlotEstimateError(Times, errorSmooth_CKF, Ps_lk,fs)
% 
% 
% fprintf('Smoothing Finished: %3.1f sec\n',toc)
% end

% ------------------------------------------------------------------------
%%% B-Plane
% ------------------------------------------------------------------------
for k = 1:length(Times)
    if norm(est_SRIF(k,1:3)) < 3*rE_SOI
        % Store Velocity as vInfinity
        vInf = est_SRIF(k,4:6); % km/s
        
        % Store the index at which vInfinity is determined
        bk = k;
        break
    end
end

figure(22); hold all
figure(33); hold all
%%% Earth Equator 
th = 0:.01:2*pi;
t = rE * cos(th);
r = rE * sin(th);
plot(t, r,'b','linewidth',1.5);
for kk = 1:4
    daysOfData = 50*kk;
    ek = find(Times==daysOfData*86400); % Time index at which I desire a B-Plane estimate
    % ek = bk; 
    [ BVec, DCM, LTOF] = BPlaneCalculator(est_SRIF(ek,1:6), uE, vInf);

    %%% Finding B-Plane Cross time
    for k = ek:length(Times)
        if Times(k) > (Times(ek) + LTOF)
            crossk = k; % Time index when B-Plane is crossed
            break
        end
    end

    % rB = DCM * BVec; % km

    rCross_B = DCM*est_SRIF(crossk,1:3)';

    PCross = Ps_SRIF(:,:,crossk);
    PCross_B = DCM * PCross(1:3,1:3) * DCM';
    % figure; hold all
    figure(22)
    error_ellipse(PCross_B(2:3,2:3),'conf',0.9973)
    PlotBoi2('T-axis, km','R-axis, km',16)

%     figure; hold all
    figure(33)
    
        
    plot(rCross_B(2),-rCross_B(3),'x','markersize',10)
    axis equal
end

figure(22)
legend('50','100','150','200')

figure(33)
legend('Earth','50','100','150','200')

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


% ------------------------------------------------------------------------
%%% Numerical Integrating Reference Orbit and STM
% ------------------------------------------------------------------------
function [ dY ] = proj1Integrator(t,Y,uS,uE,AU,Am,Pphi,c,JD0)
    %%% Size dY to fit state and all of reshaped (n^2,1) STM
    dY = zeros(7+7^2,1);
    
    %%% Unpack state
    x = Y(1); % km
    y = Y(2); % km
    z = Y(3); % km
    dx = Y(4); % km/s
    dy = Y(5); % km/s
    dz = Y(6); % km/s
    Cr = Y(7);

    %%% Reshape (n^2,1) stm to (n,n)
    stm = reshape(Y(8:end),7,7);
    
    %%% Describe Sun Location
    [rSun, vSun, mu_p] = Ephem(JD0 + t/86400,3,'EME2000'); % km, km/s, km^3/s^2
    rSun = -rSun;
    xs = rSun(1); % km
    ys = rSun(2); % km
    zs = rSun(3); % km

    %%% Build A matrix and evaluate at current state
    A = zeros(7,7);
    A(1:3,4:6) = eye(3,3);
    % Asym(4,1)
    A(4,1) = (3*uE*x*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - uS*(1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(x - xs)*sign(x - xs)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2)) - uE/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (AU^2*Am*Cr*Pphi)/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (3*AU^2*Am*Cr*Pphi*abs(x - xs)*sign(x - xs)*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(4,2)
    A(4,2) = (3*uE*x*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(y - ys)*sign(y - ys)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(y - ys)*sign(y - ys)*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(4,3)
    A(4,3) = (3*uE*x*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(z - zs)*sign(z - zs)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(z - zs)*sign(z - zs)*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(4,7)
    A(4,7) = (AU^2*Am*Pphi*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
    % Asym(5,1)
    A(5,1) = (3*uE*y*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(x - xs)*sign(x - xs)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(x - xs)*sign(x - xs)*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(5,2)
    A(5,2) = (3*uE*y*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - uS*(1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(y - ys)*sign(y - ys)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2)) - uE/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (AU^2*Am*Cr*Pphi)/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (3*AU^2*Am*Cr*Pphi*abs(y - ys)*sign(y - ys)*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(5,3)
    A(5,3) = (3*uE*y*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(z - zs)*sign(z - zs)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(z - zs)*sign(z - zs)*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(5,7)
    A(5,7) = (AU^2*Am*Pphi*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
    % Asym(6,1)
    A(6,1) = (3*uE*z*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(x - xs)*sign(x - xs)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(x - xs)*sign(x - xs)*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(6,2)
    A(6,2) = (3*uE*z*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(y - ys)*sign(y - ys)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(y - ys)*sign(y - ys)*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(6,3)
    A(6,3) = (3*uE*z*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - uS*(1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(z - zs)*sign(z - zs)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2)) - uE/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (AU^2*Am*Cr*Pphi)/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (3*AU^2*Am*Cr*Pphi*abs(z - zs)*sign(z - zs)*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(6,7)
    A(6,7) = (AU^2*Am*Pphi*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));

    %%% Calculate new STM
    stm_dot = A*stm;
    if t == 0
        aa2B = [- (uE*x)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2);
            - (uE*y)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2);
            - (uE*z)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2)]
        
        aa3B = [- uS*((x - xs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + xs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
            - uS*((y - ys)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + ys/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
            - uS*((z - zs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + zs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2))]
        
        aaSRP = [(AU^2*Am*Cr*Pphi*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
            (AU^2*Am*Cr*Pphi*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2))
            (AU^2*Am*Cr*Pphi*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2))]
    end
        
        
    %%% Creating X-dot
    % Spacecraft velocities
    dY(1:3) = [dx; dy; dz];
    %%% Using u and SRP as dynamics
    % EQM(4)
    dY(4) = (AU^2*Am*Cr*Pphi*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (uE*x)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - uS*((x - xs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + xs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
    % EQM(5)
    dY(5) = (AU^2*Am*Cr*Pphi*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (uE*y)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - uS*((y - ys)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + ys/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
    % EQM(6)
    dY(6) = (AU^2*Am*Cr*Pphi*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (uE*z)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - uS*((z - zs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + zs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
    %%% Cr has no dynamics
    dY(7) = 0;

    % Filling in reshaped (7^2,1) STM to state
    dY(8:end) = reshape(stm_dot,49,1);
end