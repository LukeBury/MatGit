clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
tic
% ========================================================================
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
% Choose 3B system and Lagrange point
% -------------------------------------------------
%%% 3B system
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Collinear lagrange point of interest
Lpoint = 2;  % 1 or 2

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec
% -------------------------------------------------
% Equilibrium Points and Jacobi Constants
% -------------------------------------------------
%%% Acquire Collinear Lagrange points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Jacobi constant of Lagrange point
[JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(Lpoint,:),[0,0,0]);

% -------------------------------------------------
% Specifying set of initial conditions
% -------------------------------------------------
%%% Choosing Jacobi constant values
% dJC0_range = [0.000012:0.000001:0.000023];
% for dJC0 = dJC0_range
% dJC0 = 0.000035; % 0.000017
JCs_cc = [0.000055,0.00005,0.000045,0.00004,0.000035,0.00003,0.000025,0.00002,0.000017];
for cc = 1:9
dJC0 = JCs_cc(cc);

JC0 = JC_Lp-dJC0;

%%% Initial positions and velocities
rn = 100;
vn = 100;
r0s = zeros(rn, 3);
v0s = zeros(vn, 3);
X0s = zeros(rn*vn, 6);
X0i = 0;

for ri = 1:rn
%     r0s(kk,:) = [L123(Lpoint,1)-.001, (1/n)*secondary.R_n*kk, 0];
%     r0s(ri,:) = [L123(Lpoint,1), 0.05*secondary.R_n, 0];
    r0s(ri,:) = [L123(Lpoint,1), (1/rn)*secondary.R_n*ri, 0];
%     r0s(kk,:) = [L123(Lpoint,1), 0, 0];

    %%% Determining initial velocity
    % JC of new starting position
    JC_r0i = JacobiConstantCalculator(secondary.MR,r0s(ri,:),[0,0,0]);

    % Difference between JC of new spot and reference spot (L-point)
    dJC = JC_r0i - JC0;
        
    for vi = 1:vn
        %%% Find necessary velocity magnitude
        if dJC < 0
            warning('Bad initial conditions')
            return
        elseif dJC == 0
            v0_n = [0, 0, 0];
        elseif dJC > 0
            v0i = sqrt(abs(dJC));

            % Find direction for velocity
%             vHat = r0s(ri,:)./norm(r0s(ri,:));
    %         vHat = R3(-vHat,(2*pi/n*kk));
            vHat = R3([0, -1, 0], -pi*vi*(1/vn));

            % Create initial velocity vector
            v0_n = vHat .* v0i;
            v0s(vi,:) = v0_n;
        end
        
        %%% Store initial conditions
        X0i = X0i + 1;
        X0s(X0i, :) = [r0s(ri,:), v0s(vi,:)];
        
    end % vi
end % ri

% -------------------------------------------------
% Integration options
% -------------------------------------------------
%%% Setting time vector
ti = 0; % sec
tf = 0.65*pi;
dt = tf/10000;
time0_n = ti:dt:tf;

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options_Apsis = odeset('Events',@event_Apsis,'RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Looping through integrations and storing periapses
% ========================================================================
% -------------------------------------------------
% Preallocation
% -------------------------------------------------
X_apses_SCR_n = [];
time_apses = [];

% for kk = 1:size(r0s,1)
for kk = 1:(rn*vn)
    %%% Creating initial state
%     X0_n = [r0s(kk,:), v0s(kk,:)];
    X0_n = X0s(kk,:);
    
    %%% Propagating the States
    [time_n, X_BCR_n, time_event, X_event, index_event] = ode45(@Int_CR3Bn, time0_n, X0_n,...
                                                       options_Apsis, secondary.MR, secondary.R_n);
%     %%% Creating SCR positions
%     r_SCR_n = X_BCR_n(:,1:3) - [1-secondary.MR, 0, 0];
    
    %%% Storing periapses and apoapses
    X_apses_SCR_n = [X_apses_SCR_n; X_event];
    time_apses = [time_apses; time_event];
        
%     %%% Plot current figure
%     figure(1); hold all
%     plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3))

end



% ========================================================================
%%% Plots
% ========================================================================
% -------------------------------------------------
% Creating countour for plot
% -------------------------------------------------
xs = linspace(L123(1,1)-2*secondary.R_n,L123(2,1)+2*secondary.R_n,600);
ys = linspace(-6*secondary.R_n,6*secondary.R_n,100);
[X, Y] = meshgrid(xs,ys);
clear xs ys

Z = zeros(size(X));
for xk = 1:size(X,1)
    for yk = 1:size(X,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[X(xk,yk), Y(xk,yk), 0] ,[0, 0, 0]);
        Z(xk,yk) = zv;
    end
end

figure(1); hold all
% % Plotting Lagrange Points
plot3(L123(1,1), L123(1,2), L123(1,3), '^', 'markeredgecolor', colors.std.maglt,...
                                       'markerfacecolor',colors.std.mag,'markersize',10);
plot3(L123(2,1), L123(2,2), L123(2,3), '^', 'markeredgecolor', colors.std.maglt,...
                                       'markerfacecolor',colors.std.mag,'markersize',10);

% Plotting contours
[C,h] = contour(X,Y,Z,[JC0, JC0],'color',colors.sch.r9(cc,:),'linewidth',1.5);

% Plotting trajectory
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3))

% Plotting apses
plot3(X_apses_SCR_n(:,1),X_apses_SCR_n(:,2),X_apses_SCR_n(:,3),'.','markersize',10,...
                                        'color',colors.sch.r9(cc,:))

% Plotting secondary
plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1,1,1],colors.std.black,1.5,0)


% figure; hold all
% plot3(X_periapses_SCR_n(:,1),X_periapses_SCR_n(:,2),X_periapses_SCR_n(:,3),'o')
% plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1,1,1],colors.std.black,1,0)



% % % find a periapsis with the X_BCR_n data and see how close it is with the one from the event function.
% % % If close enough, we can just use that method


% end


end % cc
toc
















