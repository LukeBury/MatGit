clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
tic
% ========================================================================
%%% Run Switches
% ========================================================================
run_sunSynchronousAnalysis_3B = 1; % 3B, normalized CR3BP
run_sunSynchronousAnalysis_2B = 1; % 2B, Inertial frame
run_sunSynchronousAnalysis_2BR = 0; % 2B, Rotating frame
run_potentialValidation = 0;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_color_palettes();

% ========================================================================
%%% Constant / Conditions
% ========================================================================
%%% Selecting bodies
earth = bodies.earth;
moon = bodies.moon;

%%% Normalizing constants
rNorm = moon.a; % km
tNorm = 1/moon.meanMot; % sec
vNorm = rNorm / tNorm; % km/sec

%%% Time vector
ti = 0; % sec
dt = 60*3; % sec
tf = 3600*24*5; % sec
989
% % tfm = [5,10,20,50]; % final time multiplier ... (days)

%%% [3B] Using validation-integrators and uE = 3.985996348071999e+05 = (1-moon.MR)*(rNorm^3)/(tNorm^2)
%     raanRate expected	|	raanRate observed
%05 - 1.9911e-07 rad/s	|	1.9978e-07 rad/s ... %Diff = 0.3374
%10 - 1.9911e-07 rad/s	|	1.9995e-07 rad/s ... %Diff = 0.4252
%20 - 1.9911e-07 rad/s	|	2.0001e-07 rad/s ... %Diff = 0.4536
%50 - 1.9911e-07 rad/s	|	1.9997e-07 rad/s ... %Diff = 0.4342

%%% [2B] Using uE = 3.985996348071999e+05 = (1-moon.MR)*(rNorm^3)/(tNorm^2)
%     raanRate expected	|	raanRate observed
%05 - 1.9911e-07 rad/s	|	1.9978e-07 rad/s ... %Diff = 0.3374
%10 - 1.9911e-07 rad/s	|	1.9995e-07 rad/s ... %Diff = 0.4252
%20 - 1.9911e-07 rad/s	|	2.0001e-07 rad/s ... %Diff = 0.4536
%50 - 1.9911e-07 rad/s	|	1.9997e-07 rad/s ... %Diff = 0.4342

% ========================================================================
%%% Preparing for Integration
% ========================================================================

%%% Determining SS state
%%% Defining Orbit
Re = earth.R;
alt = 800; % km
a =  Re + alt; % km
uE = (1-moon.MR)*(rNorm^3)/(tNorm^2);
n = sqrt(uE/(a^3)); % rad/s
e = 0;
p = a*(1-e^2);
J2 = earth.J2;
ta = 0;
w = 0;
raan = 0;

% Sun Sync Node Rate
ssnr = (360/365.2421897)*(pi/180)/86400; % rad/s

% inclination
i = acos((-2*(a^(7/2))*ssnr*(1-e^2)^2)/(3*Re*Re*J2*sqrt(uE)));

% ECI, SS initial state (+x axis crossing, y=z=0)
[r, v] = OE2ECI(a, e, i, raan, w, ta, uE);

% ========================================================================
%%% Sun Synch Analysis - BCR Frame
% ========================================================================
% ----------------------------------------------------------------
% Integration Setup
% ----------------------------------------------------------------
if run_sunSynchronousAnalysis_3B == 1
% % for ttt = 1:length(tfm)
% % tf = 3600*24*tfm(ttt);

%%% Converting ECI, SS state to normalized BCR;
r0_n = r./rNorm;
v0_n = (v-cross([0;0;moon.meanMot],r))./vNorm;

%%% Initial State
X0_n = [r0_n; v0_n];

%%% Setting time vector and normalizing 
time0_n = [ti:dt:tf] ./ tNorm;

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_J2 = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating the States with J2
% Int_CR3BP_J2_val(t,X,u,R1,J21)
[time_J2_n, X_BCR_J2_n] = ode45(@Int_CR3BP_J2_val, time0_n, X0_n, options_J2, moon.MR, earth.R/rNorm, earth.J2);
% [time_J2_n, X_BCR_J2_n] = ode45(@Int_CR3BP_J2_val, time0_n, X0_n, options_J2, moon.MR, earth.R/rNorm, moon.R_n, earth.J2, 0);

%%% Propagating the States without J2
[time_n, X_BCR_n] = ode45(@Int_CR3BP_val, time0_n, X0_n, options, moon.MR);
% [time_n, X_BCR_n] = ode45(@Int_CR3BP_val, time0_n, X0_n, options, moon.MR, moon.R_n);
    
% ----------------------------------------------------------------
%%% Rotating to Inertial States
% ----------------------------------------------------------------
X_BCI_n = zeros(size(X_BCR_n));
X_BCI_J2_n = zeros(size(X_BCR_J2_n));


for kk = 1:size(X_BCR_n,1)
    theta = time_n(kk);
    X_BCI_n(kk,1:3) = R3(X_BCR_n(kk,1:3),theta);
    X_BCI_J2_n(kk,1:3) = R3(X_BCR_J2_n(kk,1:3),theta);
    X_BCI_n(kk,4:6) = R3(X_BCR_n(kk,4:6),theta) + cross([0,0,1],X_BCI_n(kk,1:3));
    X_BCI_J2_n(kk,4:6) = R3(X_BCR_J2_n(kk,4:6),theta) + cross([0,0,1],X_BCI_J2_n(kk,1:3));
end

% ----------------------------------------------------------------
%%% Analyzing node drift rate
% ----------------------------------------------------------------
% ----------------------------------------
%%% Acquiring OEs
% ----------------------------------------
%%% Preallocating
raans = zeros(size(X_BCI_n,1),1);
raans_J2 = zeros(size(X_BCI_J2_n,1),1);
as = zeros(size(X_BCI_n,1),1);
as_J2 = zeros(size(X_BCI_J2_n,1),1);
es = zeros(size(X_BCI_n,1),1);
es_J2 = zeros(size(X_BCI_J2_n,1),1);
is = zeros(size(X_BCI_n,1),1);
is_J2 = zeros(size(X_BCI_J2_n,1),1);
ws = zeros(size(X_BCI_n,1),1);
ws_J2 = zeros(size(X_BCI_J2_n,1),1);


for kk = 1:size(X_BCR_n,1)
    [a,e,i,raan,w,ta] = ECI2OE(X_BCI_n(kk,1:3).*rNorm,X_BCI_n(kk,4:6).*vNorm,uE);
    raans(kk) = raan;
    as(kk) = a;
    es(kk) = e;
    is(kk) = i;
    ws(kk) = w;
    
    [a,e,i,raan_J2,w,ta] = ECI2OE(X_BCI_J2_n(kk,1:3).*rNorm,X_BCI_J2_n(kk,4:6).*vNorm,uE);
    raans_J2(kk) = raan_J2;
    as_J2(kk) = a;
    es_J2(kk) = e;
    is_J2(kk) = i;
    ws_J2(kk) = w;
    
    if raans(kk) > pi
        raans(kk) = raans(kk) - 2*pi;
    end
    
    if raans_J2(kk) > pi
        raans_J2(kk) = raans_J2(kk) - 2*pi;
    end
end
% ----------------------------------------
%%% Plotting OEs
% ----------------------------------------

% ---------------
%%% raan
% ---------------
figure; hold all
plot(time_n.*tNorm./86400, raans.*(180/pi))
plot(time_n.*tNorm./86400, raans_J2.*(180/pi))
raans_J2_expected = time_J2_n*tNorm*1.991e-7;
plot(time_J2_n.*tNorm./86400, raans_J2_expected.*(180/pi), 'k')
PlotBoi2('','raan, deg',14); legend('No J2','w/ J2','expected w/J2'); title('3B - raan')

%%% Printing RAAN-Rate results
raan_Rate = (raans_J2(end-12)-raans_J2(1))/((time_J2_n(end-12)-time_J2_n(1))*tNorm); % rad/s
fprintf('raanRate expected\t|\traanRate observed\n')
fprintf('%1.4e rad/s\t|\t%1.4e rad/s ... %%Diff = %2.4f\n\n',ssnr,raan_Rate,100*(raan_Rate - ssnr)/ssnr)

% % ---------------
% %%% e
% % ---------------
% figure
% subplot(2,1,1); hold all
% plot(time_n.*tNorm./86400, es)
% plot(time_n.*tNorm./86400, es_J2)
% PlotBoi2('','e',14); legend('No J2','w/ J2'); title('3B - e')
% subplot(2,1,2)
% plot(time_n.*tNorm./86400, es-es_J2)
% PlotBoi2('Time, days','\deltae',14)
% ---------------
%%% a
% ---------------
figure
subplot(2,1,1); hold all
plot(time_n.*tNorm./86400, as)
plot(time_n.*tNorm./86400, as_J2)
PlotBoi2('','a',14); legend('No J2','w/ J2'); title('3B - a')
subplot(2,1,2)
plot(time_n.*tNorm./86400, as-as_J2)
PlotBoi2('Time, days','\deltaa',14)
% ---------------
%%% i
% ---------------
figure
subplot(2,1,1); hold all
plot(time_n.*tNorm./86400, is.*(180/pi))
plot(time_n.*tNorm./86400, is_J2.*(180/pi))
PlotBoi2('','i, deg',14); legend('No J2','w/ Jw'); title('3B - i')
subplot(2,1,2)
plot(time_n.*tNorm./86400, (is-is_J2).*(180/pi))
PlotBoi2('Time, days','\deltai, deg',14)

% ----------------------------------------------------------------
%%% Plotting
% ----------------------------------------------------------------
% -----------------------------
% Rotating Plot
% -----------------------------
% figure; hold all
% title('Newtonian')
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'linewidth',2)
% plotBody3(earth.R/rNorm, [0, 0, 0], earth.color)
% PlotBoi3('X_n','Y_n','Z_n',14)
% axis equal
% view(0,90)
% 
% figure; hold all
% title('J2')
% plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'linewidth',2)
% plotBody3(earth.R/rNorm, [0, 0, 0], earth.color)
% PlotBoi3('X_n','Y_n','Z_n',14)
% axis equal
% view(0,90)

% -----------------------------
% Inertial Plot
% -----------------------------
figure;subplot(1,2,1); hold all
title('Newtonian')
plot3(X_BCI_n(:,1),X_BCI_n(:,2),X_BCI_n(:,3),'linewidth',1)
% plotBody3(earth.R/rNorm, [0, 0, 0], earth.color)
plotBodyTexture3(earth.R/rNorm, [0, 0, 0],earth.img)
PlotBoi3('X_n','Y_n','Z_n',14)
axis equal
view(0,90)

subplot(1,2,2); hold all
title('w/ J2')
plot3(X_BCI_J2_n(:,1),X_BCI_J2_n(:,2),X_BCI_J2_n(:,3),'linewidth',1)
% plotBody3(earth.R/rNorm, [0, 0, 0], earth.color)
plotBodyTexture3(earth.R/rNorm, [0, 0, 0],earth.img)
PlotBoi3('X_n','Y_n','Z_n',14)
axis equal
view(0,90)

% % tfm(ttt)
% % close all
% % end % end final time multiplier loop
end % end run_sunSynchronousAnalysis

% ========================================================================
% Sun Synch Analysis - ECI Frame
% ========================================================================
% ----------------------------------------------------------------
% Setting up integration and integrating
% ----------------------------------------------------------------
if run_sunSynchronousAnalysis_2B == 1
% % for ttt = 1:length(tfm)
% % tf = 3600*24*tfm(ttt);
%%% Storing initial condition
r0 = r;
v0 = v;

%%% Initial State
X0 = [r0; v0];

%%% Setting time vector
time0 = [ti:dt:tf];

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options2 = odeset('RelTol',tol,'AbsTol',tol);
options2_J2 = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating the States with J2
% Int_ECI_2B_J2(t,X,u,R,J2)
[time_J2, X_ECI_J2] = ode45(@Int_2BI_J2, time0, X0, options2_J2, uE, earth.R, earth.J2);

%%% Propagating the States without J2
[time, X_ECI] = ode45(@Int_2BI, time0, X0, options2, uE);

% ----------------------------------------------------------------
% Checking out node drift rate
% ----------------------------------------------------------------
%%% Preallocating
raans2 = zeros(size(X_ECI,1),1);
raans2_J2 = zeros(size(X_ECI,1),1);
as2 = zeros(size(X_ECI,1),1);
as2_J2 = zeros(size(X_ECI,1),1);
es2 = zeros(size(X_ECI,1),1);
es2_J2 = zeros(size(X_ECI,1),1);
is2 = zeros(size(X_ECI,1),1);
is2_J2 = zeros(size(X_ECI,1),1);
ws2 = zeros(size(X_ECI,1),1);
ws2_J2 = zeros(size(X_ECI,1),1);

for kk = 1:size(X_ECI,1)
    [a,e,i,raan,w,ta] = ECI2OE(X_ECI(kk,1:3),X_ECI(kk,4:6),uE);
    raans2(kk) = raan;
    as2(kk) = a;
    es2(kk) = e;
    is2(kk) = i;
    ws2(kk) = w;
    
    [a,e,i,raan_J2,w,ta] = ECI2OE(X_ECI_J2(kk,1:3),X_ECI_J2(kk,4:6),uE);
    raans2_J2(kk) = raan_J2;
    as2_J2(kk) = a;
    es2_J2(kk) = e;
    is2_J2(kk) = i;
    ws2_J2(kk) = w;
    
    if raans2(kk) > pi
        raans2(kk) = raans2(kk) - 2*pi;
    end
    
    if raans2_J2(kk) > pi
        raans2_J2(kk) = raans2_J2(kk) - 2*pi;
    end
end
% ----------------------------------------
%%% Plotting OEs
% ----------------------------------------

% ---------------
%%% raan
% ---------------
figure; hold all
plot(time/86400, raans2.*(180/pi)) 
plot(time/86400, raans2_J2.*(180/pi))
raans_J2_expected = time_J2*1.991e-7;%w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9w9
plot(time_J2./86400, raans_J2_expected.*(180/pi), 'k')
PlotBoi2('','raan, deg',14); legend('No J2','w/ J2','expected w/J2'); title('2B - raan')

%%% Printing RAAN-Rate results
raan_Rate2 = (raans2_J2(end-12)-raans2_J2(1))/((time_J2(end-12)-time_J2(1))); % rad/s
fprintf('raanRate expected\t|\traanRate observed\n')
fprintf('%1.4e rad/s\t|\t%1.4e rad/s ... %%Diff = %2.4f\n\n',ssnr,raan_Rate2,100*(raan_Rate2 - ssnr)/ssnr)
% % ---------------
% %%% e
% % ---------------
% figure
% subplot(2,1,1); hold all
% plot(time/86400, es2)
% plot(time/86400, es2_J2)
% PlotBoi2('','e',14); legend('No J2','w/ J2'); title('2B - e')
% subplot(2,1,2)
% plot(time/86400, es2-es2_J2)
% PlotBoi2('Time, days','\deltae',14)
% ---------------
%%% a
% ---------------
figure
subplot(2,1,1); hold all
plot(time/86400, as2)
plot(time/86400, as2_J2)
PlotBoi2('','a',14); legend('No J2','w/ J2'); title('2B - a')
subplot(2,1,2)
plot(time/86400, as2-as2_J2)
PlotBoi2('Time, days','\deltaa',14)
% ---------------
%%% i
% ---------------
figure
subplot(2,1,1); hold all
plot(time/86400, is2.*(180/pi))
plot(time/86400, is2_J2.*(180/pi))
PlotBoi2('','i, deg',14); legend('No J2','w/ Jw'); title('2B - i')
subplot(2,1,2)
plot(time/86400, (is2-is2_J2).*(180/pi))
PlotBoi2('Time, days','\deltai, deg',14)


% % close all
% % tfm(ttt) 
% % end
end % run_sunSynchronousAnalysis_ECI

% ========================================================================
% Sun Synch Analysis - ECR Frame
% ========================================================================
% ----------------------------------------------------------------
% Setting up integration and integrating
% ----------------------------------------------------------------
if run_sunSynchronousAnalysis_2BR == 1
    
%%% Storing initial condition
r0 = r;
v0 = v-cross([0;0;moon.meanMot],r0);

%%% Initial State
X0 = [r0; v0];

%%% Setting time vector
time0 = [ti:dt:tf];

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options2 = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating the States without J2
[time, X_ECR] = ode45(@Int_2BR, time0, X0, options2, uE, moon.meanMot);

% ----------------------------------------------------------------
%%% Rotating to Inertial States
% ----------------------------------------------------------------
X_ECRI = zeros(size(X_ECR));

for kk = 1:size(X_ECR,1)
    theta = time(kk)*moon.meanMot;
    X_ECRI(kk,1:3) = R3(X_ECR(kk,1:3),theta);
    X_ECRI(kk,4:6) = R3(X_ECR(kk,4:6),theta) + cross([0,0,moon.meanMot],X_ECRI(kk,1:3));
end


% ----------------------------------------------------------------
%%% Analyzing node drift rate
% ----------------------------------------------------------------
% ----------------------------------------
%%% Acquiring OEs
% ----------------------------------------
%%% Preallocating
raansRI = zeros(size(X_ECRI,1),1);
asRI = zeros(size(X_ECRI,1),1);
esRI = zeros(size(X_ECRI,1),1);
isRI = zeros(size(X_ECRI,1),1);
wsRI = zeros(size(X_ECRI,1),1);


for kk = 1:size(X_ECRI,1)
    [a,e,i,raan,w,ta] = ECI2OE(X_ECRI(kk,1:3),X_ECRI(kk,4:6),uE);
    raansRI(kk) = raan;
    asRI(kk) = a;
    esRI(kk) = e;
    isRI(kk) = i;
    
    if raansRI(kk) > pi
        raansRI(kk) = raansRI(kk) - 2*pi;
    end
    
end
% ----------------------------------------
%%% Plotting OEs
% ----------------------------------------
% ---------------
%%% raan
% ---------------
figure; hold all
plot(time/86400, raansRI.*(180/pi))
PlotBoi2('Time, days','raan, deg',14); title('2B, ECR - raan')

% % ---------------
% %%% e
% % ---------------
% figure; hold all
% plot(time/86400, esRI)
% PlotBoi2('Time, days','e',14); title('2B, ECR - e')

% ---------------
%%% a
% ---------------
figure; hold all
plot(time/86400, asRI)
PlotBoi2('Time, days','a',14); title('2B, ECR - a')

% ---------------
%%% i
% ---------------
figure; hold all
plot(time/86400, isRI.*(180/pi))
PlotBoi2('Time, days','i, deg',14); title('2B, ECR - i')



end % run_sunSynchronousAnalysis_2BR























% % % % % % % % % % % % ========================================================================
% % % % % % % % % % % %%% Potential Validation
% % % % % % % % % % % % ========================================================================
% % % % % % % % % % % if run_potentialValidation == 1
% % % % % % % % % % %     syms x y z x1 x2 r1 r2 u J21 J22 R1 R2 real
% % % % % % % % % % % %     r1 = ((x-x1)^2+y^2+z^2)
% % % % % % % % % % % %     r2 = ((x-x2)^2+y^2+z^2)
% % % % % % % % % % % %     r1^2 = ((x-x1)^2+y^2+z^2)
% % % % % % % % % % % %     r2^2 = ((x-x2)^2+y^2+z^2)
% % % % % % % % % % % %     r1^5 = ((x-x1)^2+y^2+z^2)
% % % % % % % % % % % %     r2^5 = ((x-x2)^2+y^2+z^2)
% % % % % % % % % % %     U = .5*(x^2+y^2) + (1-u)/r1 + u/r2 + (1-u)*J21*R1^2*(r1^2-3*z^2)/(2*r1^5) + u*J22*R2^2*(r2^2-3*z^2)/(2*r2^5);
% % % % % % % % % % %     
% % % % % % % % % % %     
% % % % % % % % % % %     
% % % % % % % % % % % 
% % % % % % % % % % % end % end run_potentialValidation

toc












