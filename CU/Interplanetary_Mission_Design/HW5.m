clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================
run_p1 = 0;
run_p2 = 1;

% ========================================================================
%%% HW5 - Setup and P1
% ========================================================================
% ------------------------------------
%%% Setup
% ------------------------------------
%%% Initial State
r0 = [546507.344255845; -527978.380486028; 531109.066836708];  % km
v0 = [-4.9220589268733; 5.36316523097915; -5.22166308425181]; % km/s
X0 = [r0; v0];

%%% Body parameters
Sun = bodies.sun;
Earth = bodies.earth;

Earth.u = 398600.4415; % km^3/s^2

% ------------------------------------
%%% Setup
% ------------------------------------
%%% B-Plane related quantities
vInfM = v0; % km/s
VInf = norm(vInfM); % km/s

BT_des = 13135.7982982557; % km
BR_des = 5022.26511510685; % km

% ------------------------------------
%%% Integration 
% ------------------------------------
t_i = 0; % sec
t_f = 28.445*3600*2;
n_dt = 10000;
time0 = linspace(t_i,t_f,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Integrating state
[time, X_2BI] = ode113(@Int_2BI, time0, X0, options, Earth.u);

if run_p1 == 1
    
figure; hold all
plot3(X_2BI(:,1),X_2BI(:,2),X_2BI(:,3),'color',colors.std.blue,'linewidth',2)
plotBodyTexture3(Earth.R, [0,0,0],Earth.img)
PlotBoi3('X, $km$','Y, $km$','Z, $km$',18,'LaTex')
axis equal

% ------------------------------------
%%% Initial B-Plane calculations
% ------------------------------------
%%% Gathering a v-infinity-plus
vInfP = X_2BI(end,4:6)'; % km/s

%%% Calculating B-Plane parameters
% [BT, BR, B, theta_rad] = BPlaneTargeter_vInfs(vInfM, vInfP, Earth.u, 2);
[BT, BR, B, theta_rad] = BPlaneTargeter_state(X0, Earth.u);


% ------------------------------------
%%% Targeting desired B-plane parameters
% ------------------------------------
%%% Testing perturbation sizes
n_perturbations = 1000;
perturbations = logspace(-16, 1, n_perturbations)';
perturbations = perturbations(43:end);
n_perturbations = length(perturbations);

%%% Preallocating
dBT_dvInfX = zeros(n_perturbations,1);
dBT_dvInfY = zeros(n_perturbations,1);
dBR_dvInfX = zeros(n_perturbations,1);
dBR_dvInfY = zeros(n_perturbations,1);

%%% Looping through perturbations and targeting desired conditions
for pp = 1:n_perturbations
    %%% Setting current perturbation
    pert = perturbations(pp);
    
    %%% Perturb vInf_x
    vInfM_pertX = vInfM + [pert; 0; 0];
    X0_pertX = [r0; vInfM_pertX];
    [BT_pertX, BR_pertX, B_pertX, theta_rad_pertX] = BPlaneTargeter_state(X0_pertX, Earth.u);

    dBT_dvInfX(pp) = (BT_pertX - BT)/pert;
    dBR_dvInfX(pp) = (BR_pertX - BR)/pert;
    
    
    %%% Perturb vInf_y
    vInfM_pertY = vInfM + [0; pert; 0];
    X0_pertY = [r0; vInfM_pertY];
    [BT_pertY, BR_pertY, B_pertY, theta_rad_pertY] = BPlaneTargeter_state(X0_pertY, Earth.u);

    dBT_dvInfY(pp) = (BT_pertY - BT)/pert;
    dBR_dvInfY(pp) = (BR_pertY - BR)/pert;
    
end

lineColor1 = colors.sch.d4_1(3,:);
lineColor2 = colors.sch.d4_1(2,:);

figure; hold all
p1 = plot(perturbations,dBT_dvInfX,'color',lineColor1,'linewidth',2);
p2 = plot(perturbations,dBT_dvInfY,'color',lineColor2,'linewidth',2);
setLogPlot('x') 
PlotBoi2('Perturbation Magnitude, $km/s$','Numerical Partial, $s$',18,'LaTex')
legend([p1 p2], {'$\frac{dB_T}{d\Delta V_x}$','$\frac{dB_T}{d\Delta V_y}$'},'Interpreter','LaTex','fontsize',20)


figure; hold all
p1 = plot(perturbations,dBR_dvInfX,'color',lineColor1,'linewidth',2);
p2 = plot(perturbations,dBR_dvInfY,'color',lineColor2,'linewidth',2);
setLogPlot('x') 
PlotBoi2('Perturbation Magnitude, $km/s$','Numerical Partial, $s$',18,'LaTex')
legend([p1 p2], {'$\frac{dB_R}{d\Delta V_x}$','$\frac{dB_R}{d\Delta V_y}$'},'Interpreter','LaTex','fontsize',20)

end % run_p1 == 1

% ========================================================================
%%% HW5 - P2
% ========================================================================
if run_p2 == 1
% ------------------------------------
%%% Setup
% ------------------------------------
%%% Setting a tolerance for the differential corrector 
tol_diffCorr = 1e-5;

% ------------------------------------
%%% Differential correction
% ------------------------------------
%%% Prealloctating
dvCorrect = [0; 0];
vInfM_new = v0;
currentError = 1000;
twoErrorsAgo = 0;

%%% Initial perturbation
perturbation = 2e-02;

%%% Initializing figure for results
plots = [];
names = {};
figure(10); hold all
ax = gca;
ax.YDir = 'reverse';
p = plot(BT_des, BR_des,'.','markersize',10);
plots = [plots, p];
names{end+1} = 'Desired';

%%% Entering 'while' loop that tracks error
iter = -1;
while currentError > tol_diffCorr
    iter = iter + 1;
    
    %%% Setting new initial condition
    vInfM_old = vInfM_new;
    vInfM_new = vInfM_old + [dvCorrect(1); dvCorrect(2); 0];
    X0_new = [r0; vInfM_new];
    
    % ------------------------------------
    %%% Current B-Plane calculations
    % ------------------------------------
    %%% Calculating B-Plane parameters
    [BT_current, BR_current, B_nominal, theta_rad_nominal] = BPlaneTargeter_state(X0_new, Earth.u);
    
    %%% Plotting B-plane
    figure(10)
    p = plot(BT_current, BR_current,'^','markersize',10);
    plots = [plots, p];
    names{end+1} = sprintf('Iteration %d',iter);
    
    % ------------------------------------
    %%% Perturbing and calculating partials
    % ------------------------------------
    %%% Perturb vInf_x
    vInfM_pertX_new = vInfM_new + [perturbation; 0; 0];
    X0_pertX_new = [r0; vInfM_pertX_new];
    [BT_pertX_new, BR_pertX_new, B_pertX, theta_rad_pertX] = BPlaneTargeter_state(X0_pertX_new, Earth.u);
    
    dBT_dvInfX = (BT_pertX_new - BT_current)/perturbation;
    dBR_dvInfX = (BR_pertX_new - BR_current)/perturbation;

    
    %%% Perturb vInf_y
    vInfM_pertY_new = vInfM_new + [0; perturbation; 0];
    X0_pertY_new = [r0; vInfM_pertY_new];
    [BT_pertY_new, BR_pertY_new, B_pertY, theta_rad_pertY] = BPlaneTargeter_state(X0_pertY_new, Earth.u);

    dBT_dvInfY = (BT_pertY_new - BT_current)/perturbation;
    dBR_dvInfY = (BR_pertY_new - BR_current)/perturbation;
    
    % ------------------------------------
    %%% Make corrections
    % ------------------------------------
    dB_dv = [dBR_dvInfX, dBR_dvInfY;...
             dBT_dvInfX, dBT_dvInfY];

    deltaBT = BT_des - BT_current;
    deltaBR = BR_des - BR_current;

    deltaB = [deltaBR; deltaBT];
         
    dvCorrect = inv(dB_dv) * deltaB;
    
    % ------------------------------------
    %%% Check errors
    % ------------------------------------
    BT_Error = abs(BT_current - BT_des);
    BR_Error = abs(BR_current - BR_des);
    
    oldError = currentError;
    currentError = sqrt(BT_Error + BR_Error);
    
    if currentError > oldError
        dvCorrect = -dvCorrect;
    end
        
end % currentError > tol_diffCorr

figure(10)
legend(plots,names)
PlotBoi2('B$_T$, $km$','B$_R$, $km$',18,'LaTex')

end % run_p2 == 1




































