clear
clc
close all

for bigLoop_kk = 3
tic
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
% ========================================================================
%%% Run Switches
% ========================================================================

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
% System info
% -------------------------------
%%% bodies
primary = bodies.jupiter;
secondary = bodies.europa;

% primary = bodies.saturn;
% secondary = bodies.enceladus;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates


% -------------------------------
% Creating Rerence Trajectory
% -------------------------------
%%% Setting time vector
t_0 = 0;
% t_f = 4*pi;
% t_f = 2.677618658082257; % Europa impact
% t_f = 2.287449542144822; % Enceladus impact
t_f = 5;
% n_dt = 5000;
% time0_n = linspace(t_0,t_f,n_dt);
% dt = time0_n(2)-time0_n(1);

dt = 30/tNorm; % 30 seconds
time0_n = [t_0:dt:t_f]';

% dt = 0.001;
% time0_n = [0:dt:5]';

%%% To impact Europa plume site at low angle
r0ref_n = [1.0204617015266166, -0.0054921659260262, -0.0099525155929432]';
v0ref_n = [-0.0111591144905944, -0.0104200955172558, 0.0145837924677663]';

% %%% Enceladus south pole
% r0ref_n = [1.0039914454138370, 0.0006268291981895, 0.0032971880321676]';
% v0ref_n = [-0.0012896112019572, -0.0043048383232825, -0.0030710153996081]';

%%% Setting initial state of reference
X0ref_n = [r0ref_n; v0ref_n];

%%% extra parameters
prms = struct();
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);

%%% Choosing ode tolerance
% tol = 2.22045e-14;
tol = 1e-12;

%%% Setting Options
options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Integrating reference
[time_n, Xref_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, X0ref_n, options_ImpactEscape, prms);
 
%%% Creating a backwards time vector for integration of K
time_bkwds_n = flipud(time_n);

% -------------------------------
% Creating B matrix and selecting weights for R, Q, and S
% -------------------------------
%%% B matrix (control inputs)
B = [zeros(3,3); eye(3)];

%%% Weights
% a_rf = 1000; % Final position error weights (S,1-3)
% a_vf = 1000; % Final velocity error weights (S,4-6)
% a_r  = 1; % Position error weights (Q,1-3)
% a_v  = 1; % Velocity error weights (Q,4-6)
% a_u  = 1e-6; % Control weights (R)

a_rf = 1; % Final position error weights (S,1-3)
a_vf = 1; % Final velocity error weights (S,4-6)
a_r  = 1000; % Position error weights (Q,1-3)
a_v  = 10; % Velocity error weights (Q,4-6)
a_u  = 1; % Control weights (R)


%%% Weights on controls (u' R u)
R = eye(3).*a_u;
%%% Weights on states (x' Q x)
Q = blkdiag(eye(3).*a_r, eye(3).*a_v);
%%% Weights on final states (xf' S xf
S = blkdiag(eye(3).*a_rf, eye(3).*a_vf);

% -------------------------------
% Creating time-history of S
% -------------------------------
%%% Choosing initial (final) conditions
K_tf = S;
K_tf = reshape(K_tf,6,6);

%%% Integrating
[t_bkwds, K_out_bkwds] = ode45(@Int_K, time_bkwds_n, K_tf, options,B,R,Q,Xref_n,time_n,secondary.MR);
K_out_fwds = flipud(K_out_bkwds);

% -------------------------------
% Setting up Monte Carlo and Creating Perturbations
% -------------------------------
%%% Number of MC conditons
n_MonteCarlos = 75;

%%% Position and Velocity perturbation sigmas and means
if bigLoop_kk == 1
    r_sig  = 5/rNorm; % Normalizes km
    v_sig  = 0.01/vNorm; % Normalizes km/s
elseif bigLoop_kk == 2
    r_sig  = 1000/rNorm; % Normalizes km
    v_sig  = 0.250/vNorm; % Normalizes km/s
elseif bigLoop_kk == 3
    r_sig  = 0.05/rNorm; % Normalizes km
    v_sig  = 0.001/vNorm; % Normalizes km/s
end
r_mean = zeros(2,1)./rNorm;
v_mean = zeros(3,1)./vNorm;

%%% Creating vectors of perturbations
rPerts = randn(2,n_MonteCarlos)*r_sig + r_mean; % Not perturbing x
vPerts = randn(3,n_MonteCarlos)*v_sig + v_mean;

%%% Preallocating data structures
positionErrors_component_m   = zeros(n_MonteCarlos,3);
velocityErrors_component_mps = zeros(n_MonteCarlos,3);
positionErrors_m   = zeros(n_MonteCarlos,1);
velocityErrors_mps = zeros(n_MonteCarlos,1);
totalDeltaV_mps    = zeros(n_MonteCarlos,1);

%%% Beginning Monte Carlo!
for mm = 1:n_MonteCarlos
rPert = rPerts(:,mm);
vPert = vPerts(:,mm);
% -------------------------------
% Integrating state with guidance algorithm
% -------------------------------
%%% Perturbation
dX0 = [0; rPert; vPert];
% dX0 = [10/rNorm; 0;0;0;0;0];

%%% Creating perturbed X0
X0_n = X0ref_n + dX0;

%%% Integrating
[t_out, X_guidance] = ode45(@Int_State, time_n,X0_n,options,B,R,K_out_fwds,Xref_n,time_n,secondary.MR,dt);

% -------------------------------
% Recreating controls
% -------------------------------
%%% Preallocating
us_mpss = zeros(length(t_out)-1,3);
for kk = 1:(length(t_out)-1)
    %%% Current K matrix
    K = reshape(K_out_fwds(kk,:),6,6);
    
    %%% Computing Controls
    u_t = -(inv(R)*B'*K)*(X_guidance(kk,:)' - Xref_n(kk,:)');
    
    %%% In meters/second^2
    us_mpss(kk,:) = u_t' .*(rNorm/(tNorm^2))*1000;
end

%%% Storing results of current run
positionErrors_component_m(mm,:)   = (X_guidance(end,1:3) - Xref_n(end,1:3)).*rNorm*1000;
velocityErrors_component_mps(mm,:) = (X_guidance(end,4:6) - Xref_n(end,4:6)).*vNorm*1000;

positionErrors_m(mm)   = norm(positionErrors_component_m(mm,:));
velocityErrors_mps(mm) = norm(velocityErrors_component_mps(mm,:));
totalDeltaV_mps(mm)    = trapz(rowNorm(us_mpss));

fprintf('Monte Carlo Completion Status: %1.0f%%\n',mm*100/n_MonteCarlos)
end % Monte Carlo
% 
% figure
% histfit(positionErrors_m)
% PlotBoi2('Final Position Error, $m$','Frequency',16,'LaTex')
% 
% figure
% histfit(velocityErrors_mps)
% PlotBoi2('Final Velocity Error, $m/s$','Frequency',16,'LaTex')
% 
% figure
% histfit(totalDeltaV_mps)
% PlotBoi2('Total $\Delta$V, $m/s$','Frequency',16,'LaTex')

% figure
% subplot(3,1,1)
% histfit(positionErrors_component_m(:,1))
% PlotBoi2('$\Delta x_f$, $m$','Frequency',16,'LaTex')
% subplot(3,1,2)
% histfit(positionErrors_component_m(:,2))
% PlotBoi2('$\Delta y_f$, $m$','Frequency',16,'LaTex')
% subplot(3,1,3)
% histfit(positionErrors_component_m(:,3))
% PlotBoi2('$\Delta z_f$, $m$','Frequency',16,'LaTex')

std(positionErrors_m)
mean(positionErrors_m)
std(velocityErrors_mps)
mean(velocityErrors_mps)
std(totalDeltaV_mps)
mean(totalDeltaV_mps)

if bigLoop_kk == 1
    save('Guidance_MC_smallDist.mat','positionErrors_component_m','velocityErrors_component_mps','positionErrors_m','velocityErrors_mps','totalDeltaV_mps','rPerts','vPerts')
elseif bigLoop_kk == 2
    save('Guidance_MC_hugeDist.mat','positionErrors_component_m','velocityErrors_component_mps','positionErrors_m','velocityErrors_mps','totalDeltaV_mps','rPerts','vPerts')
elseif bigLoop_kk == 3
    save('Guidance_MC_RealSmallDist.mat','positionErrors_component_m','velocityErrors_component_mps','positionErrors_m','velocityErrors_mps','totalDeltaV_mps','rPerts','vPerts')
end

toc
clear

end % bigLoop_kk = 1:2

% ========================================================================
% ========================================================================
%%% Functions
% ========================================================================
% ========================================================================
function [dK] = Int_K(t,K,B,R,Q,Xref,tref,u)

%%% Interpolating reference state for current values
Xref_now = interp1(tref,Xref,t);

%%% Getting A-Matrix
[A] = get_Amat_CR3BP(u, Xref_now(1), Xref_now(2), Xref_now(3));

K = reshape(K,6,6);

%%% Preallocate state output
dK = -K*A + K*B*inv(R)*B'*K - Q - A'*K;
dK = reshape(dK,36,1);
end


function [dX] = Int_State(t,X,B,R,K,Xref,tref,u,dt)

if t < (tref(end)-dt)
    % -------------------------------
    %%% Interpolating reference and calculating A matrix
    % -------------------------------
    %%% Interpolating reference state for current values
    Xref_now = interp1(tref,Xref,t);
    Xref_now = Xref_now';
    % 
    % %%% Getting A-Matrix
    % [A] = get_Amat_CR3BP(u, Xref_now(1), Xref_now(2), Xref_now(3));

    % -------------------------------
    %%% Creating controls
    % -------------------------------
    K_now = interp1(tref,K,t);
    K = reshape(K_now,6,6);

    deltaX = X - Xref_now;
    u_t = -(inv(R)*B'*K)*deltaX;
    
    % -------------------------------
    %%% Integrating state
    % -------------------------------
    %%% Preallocate state output
    dX = zeros(6,1);

    %%% Distances to primary (1) and secondary (2) bodies
    r1 = sqrt((X(1)+u)^2 + X(2)^2 + X(3)^2);
    r2 = sqrt((X(1)+u-1)^2 + X(2)^2 + X(3)^2);

    %%% Equations of Motion
    ddx = 2*X(5) + X(1) - (1-u)*(X(1)+u)/(r1^3) - u*(X(1)+u-1)/(r2^3) + u_t(1);
    ddy = -2*X(4) + X(2) -((1-u)/(r1^3) + u/(r2^3))*X(2) + u_t(2);
    ddz = -((1-u)/(r1^3) + u/(r2^3))*X(3) + u_t(3);

    %%% Output the derivative of the state
    dX(1:3) = [X(4); X(5); X(6)]; % km/s
    dX(4:6) = [ddx; ddy; ddz]; % km/s^2
else
    % -------------------------------
    %%% Integrating state
    % -------------------------------
    %%% Preallocate state output
    dX = zeros(6,1);

    %%% Distances to primary (1) and secondary (2) bodies
    r1 = sqrt((X(1)+u)^2 + X(2)^2 + X(3)^2);
    r2 = sqrt((X(1)+u-1)^2 + X(2)^2 + X(3)^2);

    %%% Equations of Motion
    ddx = 2*X(5) + X(1) - (1-u)*(X(1)+u)/(r1^3) - u*(X(1)+u-1)/(r2^3);
    ddy = -2*X(4) + X(2) -((1-u)/(r1^3) + u/(r2^3))*X(2);
    ddz = -((1-u)/(r1^3) + u/(r2^3))*X(3);

    %%% Output the derivative of the state
    dX(1:3) = [X(4); X(5); X(6)]; % km/s
    dX(4:6) = [ddx; ddy; ddz]; % km/s^2
end
end













