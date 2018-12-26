clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic

% ========================================================================
%%% Run Switches
% ========================================================================
run_part_a = 0;
run_part_b = 0;
run_part_c = 1;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------
% Dynamics
% -------------------------------
A0 = [0; 1; 0];
B0 = [0; 0; 0];

% -------------------------------
% Initial State
% -------------------------------
r0 = [8276; 5612; 5]; % km
v0 = [-3.142; 4.672; 0]; % km/s

scState0 = [r0; v0];

% -------------------------------
% Target Orbital Elements
% -------------------------------
a_tgt = 10000; % km
i_tgt = 0; % deg
e_tgt = 0.0001;
w_tgt = 35; % deg


% -------------------------------
% Integration options
% -------------------------------
%%% Choosing ode tolerance
tol = 2.22045e-14;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------
% Other
% -------------------------------
%%% Thrust acceleration (constant)
aThrust = 30e-3; % km/s^2 

%%% Initial guess for the burn-out time
tf_guess0 = 23; % sec

%%% System
earth = bodies.earth;














% ========================================================================
%%% Part a - Simulate burn without any guidance
% ========================================================================
% -------------------------------
% Preparing for integration
% -------------------------------
%%% Setting time vector
t_i   = 0;
n_dt = 1000;
time_burn0 = linspace(t_i,tf_guess0,n_dt);

% -------------------------------
% Integrating
% -------------------------------
%%% Integrating thrusting section
[time_burn0, scState_2BI_burn0] = ode113(@Int_2BI_BLT, time_burn0, scState0, options, earth.u, A0, B0, aThrust);

%%% Integrating post-thrust orbit
[time_postBurn0, scState_2BI_postBurn0] = ode113(@Int_2BI, linspace(t_i,20*60,1000), scState_2BI_burn0(end,:)', options, earth.u);

if run_part_a == 1
% -------------------------------
% Plotting post-burn0 trajectory
% -------------------------------
figure; hold all
plot3(scState_2BI_burn0(:,1),scState_2BI_burn0(:,2),scState_2BI_burn0(:,3),'r','linewidth',3)
plot3(scState_2BI_burn0(1,1),scState_2BI_burn0(1,2),scState_2BI_burn0(1,3),'ro','linewidth',1.5,'markersize',3)
plot3(scState_2BI_postBurn0(:,1),scState_2BI_postBurn0(:,2),scState_2BI_postBurn0(:,3),'b','linewidth',1.5)
plotBodyTexture3(earth.R, [0, 0, 0], earth.img);
PlotBoi3('x, $km$','y, $km$','z, $km$', 16,'LaTex')
axis equal
view(0,90)

% -------------------------------
% Post-burn0 orbital elements
% -------------------------------
[a,e,i,raan,w,ta] = RV2COE(scState_2BI_postBurn0(1,1:3),scState_2BI_postBurn0(1,4:6),earth.u);

fprintf('Post-burn OEs:\na = %1.3f km\ne = %1.3f\ni = %1.4f deg\nRAAN = %1.4f deg\nw = %1.3f deg\nta = %1.3f\n\n',...
    a, e, i, raan, w, ta);

fprintf('Post-burn OE Errors:\nda = %1.3f km\nde = %1.3f\ndi = %1.4f deg\ndw = %1.3f deg\n\n',...
    a-a_tgt, e-e_tgt, i-i_tgt, w-w_tgt);

end % run_part_a
















% ========================================================================
%%% Part b - Guidance
% ========================================================================
if run_part_b == 1
% -------------------------------
% Time vector to step through
% -------------------------------
dt = 0.1;
times = 0:dt:tf_guess0;

% -------------------------------
% Desired parameters
% -------------------------------
%%% Target perigee info
% Position
rp_mag = a_tgt*(1-e_tgt);
rp_vec = [rp_mag*cosd(w_tgt); rp_mag*sind(w_tgt); 0];
rpHat = rp_vec./norm(rp_vec);

% Velocity
[ vp_mag ] = visviva_v( rp_mag, a_tgt, earth.u);
vp_vec = [-sind(w_tgt); cosd(w_tgt); 0].*vp_mag;
vpHat = vp_vec./norm(vp_vec);

% Angluar momentum
h_vec_tgt = cross(rp_vec,vp_vec);
hHat_tgt = h_vec_tgt./norm(h_vec_tgt);

%%% Values to perturb states
hScale  = 1e-7; % For scaling states
hSafety = 1e-8; % To avoid 0 if state equals 0
% hScale  = 1e-11; % For scaling states
% hSafety = 1e-12; % To avoid 0 if state equals 0
% -------------------------------
% Running predictor-corrector
% -------------------------------
%%% Initial guidance state
X_i = [A0; B0; tf_guess0];

%%% Initial spacecraft state
scState0_i = scState0;

%%% Initial time
time_i = 0;

%%% Preallocating
% Norms of epsilon vectors
constraintNorms = [];
% Spacecraft states
scStates = zeros(length(times),6);
scStates(1,:) = scState0_i';
scTimes = zeros(length(times),1);
gammaWinners = [];

% Counter for iteration loops
iterCounter = 1;

%%% Beginning guidance
while time_i < X_i(7)
    %%% Changing dt if necessary
% %     if time_i > 25
% %         dt = 0.1;
% %     end
%     if time_i + dt > X_i(7)
%         dt = (X_i(7) - time_i)
%     end
%     
%     %%% Incrementing time and storing current state
%     time_i = time_i + dt;
    
    iterCounter = iterCounter + 1;
    % -------------------------------
    % Perturbing states and calculating epsilon vectors from end-states
    % -------------------------------
    %%% Perturbed integration time (current time to estimated end time)
    time_int_i = [time_i, X_i(7)];
    
    %%% Vector of perturbations
    hVec = X_i.*hScale + hSafety;
    
    %%% Integrating the un-perterbed trajectory to estimated end-state
    [epsVec_nom] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
 
    %%% Perturbing X1 and repeating
    X_pert_i = [X_i(1) + hVec(1); X_i(2:7)];
    [epsVec1] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
    
    %%% Perturbing X2 and repeating
    X_pert_i = [X_i(1); X_i(2) + hVec(2); X_i(3:7)];
    [epsVec2] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
    
    %%% Perturbing X3 and repeating
    X_pert_i = [X_i(1:2); X_i(3) + hVec(3); X_i(4:7)];
    [epsVec3] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
    
    %%% Perturbing X4 and repeating
    X_pert_i = [X_i(1:3); X_i(4) + hVec(4); X_i(5:7)];
    [epsVec4] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
    
    %%% Perturbing X5 and repeating
    X_pert_i = [X_i(1:4); X_i(5) + hVec(5); X_i(6:7)];
    [epsVec5] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
    
    %%% Perturbing X6 and repeating
    X_pert_i = [X_i(1:5); X_i(6) + hVec(6); X_i(7)];
    [epsVec6] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
    
    %%% Perturbing X7 and repeating
    X_pert_i = [X_i(1:6); X_i(7) + hVec(7)];
    [epsVec7] = returnEpsilonVector([time_i, X_i(7)+hVec(7)], scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

    % -------------------------------
    % Building Jacobian
    % -------------------------------
    J = zeros(7,7);
    J(:,1) = (epsVec1 - epsVec_nom)/hVec(1);
    J(:,2) = (epsVec2 - epsVec_nom)/hVec(2);
    J(:,3) = (epsVec3 - epsVec_nom)/hVec(3);
    J(:,4) = (epsVec4 - epsVec_nom)/hVec(4);
    J(:,5) = (epsVec5 - epsVec_nom)/hVec(5);
    J(:,6) = (epsVec6 - epsVec_nom)/hVec(6);
    J(:,7) = (epsVec7 - epsVec_nom)/hVec(7);
    
    % -------------------------------
    % Testing gamma values
    % -------------------------------
    %%% Gamma values to test
    gammas = [0.01, .1, .3, .5,.6, .7, .9, 1];
    
    %%% Preallocating vectors for norms of new epsilon vectors and for new X states
    epsNorms = zeros(1,length(gammas));
    X_news   = zeros(7,length(gammas));
    
    %%% Counter for gamma loop
    cc = 0;
    
    %%% Looping through gamma values
    for gamma = gammas
        %%% Counting loop
        cc = cc + 1;
        
        %%% Determining corrected X_i
        X_new = X_i - J\epsVec_nom.*gamma;

        %%% Calculating epsilon vector from corrected state
        [epsVec_new] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_new, aThrust,a_tgt,...
        e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);
        
        %%% Storing results from current gamma value
        epsNorms(cc) = norm(epsVec_new);
        X_news(:,cc) = X_new;
        
    end
    
    %%% Finding minimum epsilon-norm from gamma testing and grabbing the corresponding epsNorm value and X_i
    min_Index = find(epsNorms == min(epsNorms));
    epsNorm = epsNorms(min_Index);
    X_i = X_news(:,min_Index);
    
    %%% Storing winning gamma
    gammaWinners = [gammaWinners, gammas(min_Index)];
    
    %%% Tracking norm of constraint vector
    constraintNorms = [constraintNorms, epsNorm];
    
    %%% Changing dt so I end on X_i(7)
    if time_i + dt > X_i(7)
        dt = (X_i(7) - time_i);
    end
    
    %%% Update actual spacecraft state using computed guidance
    [times_dt, scState_dt]   = ode113(@Int_2BI_BLT, [time_i, time_i + dt], scState0_i, options,...
         earth.u, X_i(1:3), X_i(4:6), aThrust);
    scState0_i = scState_dt(end,:)';
    
    scStates(iterCounter,:) = scState0_i;
    scTimes(iterCounter) = times_dt(end);
    
    %%% Incrementing time and storing current state
    time_i = time_i + dt;

end

[a,e,i,raan,w,ta] = RV2COE(scStates(end,1:3),scStates(end,4:6),earth.u);
fprintf('Post-guidance OEs:\na = %1.5f km\ne = %1.5f\ni = %1.5f deg\nRAAN = %1.5f deg\nw = %1.5f deg\nta = %1.5f deg\n\n',...
    a, e, i, raan, w, ta);
fprintf('Post-guidance OE Errors:\nda = %1.5f km\nde = %1.5f\ndi = %1.5f deg\ndw = %1.5f deg\n\n',...
    a-a_tgt, e-e_tgt, i-i_tgt, w-w_tgt);


figure
plot(scTimes(2:end),constraintNorms)
PlotBoi2('Time, $sec$','Error Vector Norm',16,'LaTex')
setLogPlot
end % run_part_b














% ========================================================================
%%% Part c - Monte Carlo
% ========================================================================
if run_part_c == 1
%%% Number of monte carlo runs
n_MonteCarlos = 100;

%%% Determining thrust errors
sig = aThrust * 0.05 / 3;
thrustDispersions = randn(1,n_MonteCarlos)*sig;
thrustDispersions = [0, thrustDispersions]; % First case will be nominal

%%% Preallocating vectors for results
MCresults_da = zeros(1,n_MonteCarlos); % semi major axis
MCresults_de = zeros(1,n_MonteCarlos); % eccentricity
MCresults_di = zeros(1,n_MonteCarlos); % inclination
MCresults_dw = zeros(1,n_MonteCarlos); % argument of periapsis
MCresults_ta = zeros(1,n_MonteCarlos); % true anomaly
MCresults_bt = zeros(1,n_MonteCarlos); % burn time

DispersionLoopCounter = 1;
for thrustDispersion = thrustDispersions
    % -------------------------------
    % Time vector to step through
    % -------------------------------
    dt = 0.1;
    times = 0:dt:tf_guess0;

    % -------------------------------
    % Desired parameters
    % -------------------------------
    %%% Target perigee info
    % Position
    rp_mag = a_tgt*(1-e_tgt);
    rp_vec = [rp_mag*cosd(w_tgt); rp_mag*sind(w_tgt); 0];
    rpHat = rp_vec./norm(rp_vec);

    % Velocity
    [ vp_mag ] = visviva_v( rp_mag, a_tgt, earth.u);
    vp_vec = [-sind(w_tgt); cosd(w_tgt); 0].*vp_mag;
    vpHat = vp_vec./norm(vp_vec);

    % Angluar momentum
    h_vec_tgt = cross(rp_vec,vp_vec);
    hHat_tgt = h_vec_tgt./norm(h_vec_tgt);

    %%% Values to perturb states
    hScale  = 1e-7; % For scaling states
    hSafety = 1e-8; % To avoid 0 if state equals 0
    % hScale  = 1e-11; % For scaling states
    % hSafety = 1e-12; % To avoid 0 if state equals 0
    % -------------------------------
    % Running predictor-corrector
    % -------------------------------
    %%% Initial guidance state
    X_i = [A0; B0; tf_guess0];

    %%% Initial spacecraft state
    scState0_i = scState0;

    %%% Initial time
    time_i = 0;

    %%% Preallocating
    % Norms of epsilon vectors
    constraintNorms = [];
    % Spacecraft states
    scStates = zeros(length(times),6);
    scStates(1,:) = scState0_i';
    scTimes = zeros(length(times),1);
    gammaWinners = [];

    % Counter for iteration loops
    iterCounter = 1;

    %%% Beginning guidance
    while time_i < X_i(7)
        %%% Changing dt if necessary
    % %     if time_i > 25
    % %         dt = 0.1;
    % %     end
    %     if time_i + dt > X_i(7)
    %         dt = (X_i(7) - time_i)
    %     end
    %     
    %     %%% Incrementing time and storing current state
    %     time_i = time_i + dt;

        iterCounter = iterCounter + 1;
        % -------------------------------
        % Perturbing states and calculating epsilon vectors from end-states
        % -------------------------------
        %%% Perturbed integration time (current time to estimated end time)
        time_int_i = [time_i, X_i(7)];

        %%% Vector of perturbations
        hVec = X_i.*hScale + hSafety;

        %%% Integrating the un-perterbed trajectory to estimated end-state
        [epsVec_nom] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        %%% Perturbing X1 and repeating
        X_pert_i = [X_i(1) + hVec(1); X_i(2:7)];
        [epsVec1] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        %%% Perturbing X2 and repeating
        X_pert_i = [X_i(1); X_i(2) + hVec(2); X_i(3:7)];
        [epsVec2] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        %%% Perturbing X3 and repeating
        X_pert_i = [X_i(1:2); X_i(3) + hVec(3); X_i(4:7)];
        [epsVec3] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        %%% Perturbing X4 and repeating
        X_pert_i = [X_i(1:3); X_i(4) + hVec(4); X_i(5:7)];
        [epsVec4] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        %%% Perturbing X5 and repeating
        X_pert_i = [X_i(1:4); X_i(5) + hVec(5); X_i(6:7)];
        [epsVec5] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        %%% Perturbing X6 and repeating
        X_pert_i = [X_i(1:5); X_i(6) + hVec(6); X_i(7)];
        [epsVec6] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        %%% Perturbing X7 and repeating
        X_pert_i = [X_i(1:6); X_i(7) + hVec(7)];
        [epsVec7] = returnEpsilonVector([time_i, X_i(7)+hVec(7)], scState0_i, options, earth.u, X_pert_i, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

        % -------------------------------
        % Building Jacobian
        % -------------------------------
        J = zeros(7,7);
        J(:,1) = (epsVec1 - epsVec_nom)/hVec(1);
        J(:,2) = (epsVec2 - epsVec_nom)/hVec(2);
        J(:,3) = (epsVec3 - epsVec_nom)/hVec(3);
        J(:,4) = (epsVec4 - epsVec_nom)/hVec(4);
        J(:,5) = (epsVec5 - epsVec_nom)/hVec(5);
        J(:,6) = (epsVec6 - epsVec_nom)/hVec(6);
        J(:,7) = (epsVec7 - epsVec_nom)/hVec(7);

        % -------------------------------
        % Testing gamma values
        % -------------------------------
        %%% Gamma values to test
        gammas = [0.01, .1, .3, .5,.6, .7, .9, 1];

        %%% Preallocating vectors for norms of new epsilon vectors and for new X states
        epsNorms = zeros(1,length(gammas));
        X_news   = zeros(7,length(gammas));

        %%% Counter for gamma loop
        cc = 0;

        %%% Looping through gamma values
        for gamma = gammas
            %%% Counting loop
            cc = cc + 1;

            %%% Determining corrected X_i
            X_new = X_i - J\epsVec_nom.*gamma;

            %%% Calculating epsilon vector from corrected state
            [epsVec_new] = returnEpsilonVector(time_int_i, scState0_i, options, earth.u, X_new, aThrust,a_tgt,...
            e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag);

            %%% Storing results from current gamma value
            epsNorms(cc) = norm(epsVec_new);
            X_news(:,cc) = X_new;

        end

        %%% Finding minimum epsilon-norm from gamma testing and grabbing the corresponding epsNorm value and X_i
        min_Index = find(epsNorms == min(epsNorms));
        epsNorm = epsNorms(min_Index);
        X_i = X_news(:,min_Index);

        %%% Storing winning gamma
        gammaWinners = [gammaWinners, gammas(min_Index)];

        %%% Tracking norm of constraint vector
        constraintNorms = [constraintNorms, epsNorm];

        %%% Changing dt so I end on X_i(7)
        if time_i + dt > X_i(7)
            dt = (X_i(7) - time_i);
        end

        %%% Update actual spacecraft state using computed guidance
        [times_dt, scState_dt]   = ode113(@Int_2BI_BLT_partC, [time_i, time_i + dt], scState0_i, options,...
             earth.u, X_i(1:3), X_i(4:6), aThrust, thrustDispersion);
        scState0_i = scState_dt(end,:)';

        scStates(iterCounter,:) = scState0_i;
        scTimes(iterCounter) = times_dt(end);

        %%% Incrementing time and storing current state
        time_i = time_i + dt;

    end

    [a,e,i,raan,w,ta] = RV2COE(scStates(end,1:3),scStates(end,4:6),earth.u);
%     fprintf('Post-guidance OEs:\na = %1.5f km\ne = %1.5f\ni = %1.5f deg\nRAAN = %1.5f deg\nw = %1.5f deg\nta = %1.5f deg\n\n',...
%         a, e, i, raan, w, ta);
%     fprintf('Post-guidance OE Errors:\nda = %1.5f km\nde = %1.5f\ndi = %1.5f deg\ndw = %1.5f deg\n\n',...
%         a-a_tgt, e-e_tgt, i-i_tgt, w-w_tgt);
% 
% 
%     figure
%     plot(scTimes(2:end),constraintNorms)
%     PlotBoi2('Time, $sec$','Error Vector Norm',16,'LaTex')
%     setLogPlot

%%% if DispersionLoopCounter == 1, then there was no dispersion... nominal values were used
if DispersionLoopCounter == 1
    nomResults_da = a-a_tgt; % semi major axis;
    nomResults_de = e-e_tgt; % eccentricity;
    nomResults_di = i-i_tgt; % inclination;
    nomResults_dw = w-w_tgt; % argument of periapsis;
    nomResults_ta = ta; % true anomaly;
    nomResults_bt = X_i(7); % burn time;
elseif DispersionLoopCounter > 1
    MCresults_da(DispersionLoopCounter-1) = a-a_tgt; % semi major axis
    MCresults_de(DispersionLoopCounter-1) = e-e_tgt; % eccentricity
    MCresults_di(DispersionLoopCounter-1) = i-i_tgt; % inclination
    MCresults_dw(DispersionLoopCounter-1) = w-w_tgt; % argument of periapsis
    MCresults_ta(DispersionLoopCounter-1) = ta; % true anomaly
    MCresults_bt(DispersionLoopCounter-1) = X_i(7); % burn time
end

%%% printing status
clc
fprintf('%1.1f%% complete ... %1.2f minutes\n',DispersionLoopCounter*100/length(thrustDispersions), toc/60)

%%% Incrementing counter
DispersionLoopCounter = DispersionLoopCounter + 1;
end % thrustDispersion = thrustDispersions

%%% Calculating standard deviations and means of results
da_std  = std(MCresults_da);
da_mean = mean(MCresults_da);
de_std  = std(MCresults_de);
de_mean = mean(MCresults_de);
di_std  = std(MCresults_di);
di_mean = mean(MCresults_di);
dw_std  = std(MCresults_dw);
dw_mean = mean(MCresults_dw);
ta_std  = std(MCresults_ta);
ta_mean = mean(MCresults_ta);
bt_std  = std(MCresults_bt);
bt_mean = mean(MCresults_bt);

%%% Plotting histograms of results
figure
histfit(MCresults_da)
PlotBoi2('Semi-Major Axis, $km$','Frequency',16,'LaTex')
figure
histfit(MCresults_de)
PlotBoi2('Eccentricity','Frequency',16,'LaTex')
figure
histfit(MCresults_di)
PlotBoi2('Inclination, $^\circ$','Frequency',16,'LaTex')
figure
histfit(MCresults_dw)
PlotBoi2('Argument of Periapsis, $^\circ$','Frequency',16,'LaTex')
figure
histfit(MCresults_ta)
PlotBoi2('Injection True Anomaly, $^\circ$','Frequency',16,'LaTex')
figure
histfit(MCresults_bt)
PlotBoi2('Burn Time, $sec$','Frequency',16,'LaTex')
end % run_part_c


















toc
% ========================================================================
%%% Functions
% ========================================================================
function [epsVec] = returnEpsilonVector(time_int, scState0_i, options, u, X_pert_i, aThrust,a_tgt,e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag)
%%% Integrate spacecraft state forward from current time to predicted end time ... perturb X(1)
[~, scState_pert_i] = ode113(@Int_2BI_BLT, time_int, scState0_i, options, u, X_pert_i(1:3), X_pert_i(4:6), aThrust);

%%% Grab final state (tf)
scState_tf_pert_i = scState_pert_i(end,:)';

%%% Calculate epsilon vector with values from this final state
epsVec = calculateEpsilonVector(u,a_tgt,e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag,scState_tf_pert_i,X_pert_i(1:3),X_pert_i(4:6),X_pert_i(7));  
end

function [epsVec] = calculateEpsilonVector(u,a_tgt,e_tgt,hHat_tgt,rp_vec,rpHat,rp_mag,vp_vec,vpHat,vp_mag,scState_tf,A,B,tf_guess)
%%% Necessary parameters
r_tf = [scState_tf(1); scState_tf(2); scState_tf(3)];
v_tf = [scState_tf(4); scState_tf(5); scState_tf(6)];
ta = atan(dot(r_tf,vpHat)/dot(r_tf,rpHat));
p = a_tgt*(1-e_tgt^2);
rR = p/(1 + e_tgt*cos(ta));
vR = (-1/rp_mag)*sqrt(u/p)*sin(ta).*rp_vec + (1-rp_mag*(1-cos(ta))/p).*vp_vec;

%%% Building Lambda
lam_tf = A + B.*tf_guess;

%%% Gravity vector
g_tf = -u*r_tf./(norm(r_tf)^3);

%%% Rotatation matrix for lambda and B to go from R-Theta-H frame to Inertial
rHat_tf = r_tf./norm(r_tf);
h_tf = cross(r_tf,v_tf);
hHat_tf = h_tf./norm(h_tf);
thHat_tf = cross(hHat_tf,rHat_tf)./norm(cross(hHat_tf,rHat_tf));
RThetaH_to_Inertial = [rHat_tf, thHat_tf, hHat_tf];
if size(RThetaH_to_Inertial) ~= [3,3]
    warning('Dimension Issues')
end

%%% Building epsilon vector
eps1 = vR(1) - v_tf(1);
eps2 = vR(2) - v_tf(2);
eps3 = vR(3) - v_tf(3);
eps4 = rR - norm(r_tf);
eps5 = dot(hHat_tgt,rHat_tf);
eps6 = dot(RThetaH_to_Inertial*lam_tf, g_tf) - dot(RThetaH_to_Inertial*B,scState_tf(4:6));
eps7 = sqrt(dot(lam_tf,lam_tf)) - 1;

epsVec = [eps1; eps2; eps3; eps4; eps5; eps6; eps7];
end


function [dX] = Int_2BI_BLT(t,X,u,A,B,tMag)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity and constant in-track acceleration
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body
%          tMag - in-track thrust value

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances from body to spacecraft
r = sqrt(x^2 + y^2 + z^2);

%%% Finding uHat
uHat = (A + B.*t)./norm(A + B.*t);

%%% Finding rotation matrix to go from R-Theta-H frame to Inertial
rHat = [x;y;z]./norm([x;y;z]);
hHat = cross([x;y;z],[dx;dy;dz])./norm(cross([x;y;z],[dx;dy;dz]));
thHat = cross(hHat,rHat)./norm(cross(hHat,rHat));
RThetaH_to_Inertial = [rHat, thHat, hHat];

%%% Rotating uHat into inertial
uHat = RThetaH_to_Inertial * uHat;

%%% Finding acceleration from thrust
aThrust = uHat .* tMag;

%%% Equations of Motion
ddx = x*(-u/r^3) + aThrust(1);
ddy = y*(-u/r^3) + aThrust(2);
ddz = z*(-u/r^3) + aThrust(3);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end


function [dX] = Int_2BI_BLT_partC(t,X,u,A,B,tMag,thrustDispersion)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity and constant in-track acceleration
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body
%          tMag - in-track thrust value

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances from body to spacecraft
r = sqrt(x^2 + y^2 + z^2);

%%% Finding uHat
uHat = (A + B.*t)./norm(A + B.*t);

%%% Finding rotation matrix to go from R-Theta-H frame to Inertial
rHat = [x;y;z]./norm([x;y;z]);
hHat = cross([x;y;z],[dx;dy;dz])./norm(cross([x;y;z],[dx;dy;dz]));
thHat = cross(hHat,rHat)./norm(cross(hHat,rHat));
RThetaH_to_Inertial = [rHat, thHat, hHat];

%%% Rotating uHat into inertial
uHat = RThetaH_to_Inertial * uHat;

%%% Finding acceleration from thrust
aThrust = uHat .* (tMag + thrustDispersion);

%%% Equations of Motion
ddx = x*(-u/r^3) + aThrust(1);
ddy = y*(-u/r^3) + aThrust(2);
ddz = z*(-u/r^3) + aThrust(3);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end









