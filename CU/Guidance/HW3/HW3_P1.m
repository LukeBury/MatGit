clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Run Switches
% ========================================================================
%%% 1-a
run_ProNav                = 0;
run_AugmentedProNav       = 0;
run_SoftConstraintOptimal = 0;
%%% 1-b
run_AugmentedProNav_1b    = 1;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% P1-a
% ========================================================================
% -------------------------------------------------
% Setup
% -------------------------------------------------
%%% Dynamics
A = [0 1 0; 0 0 1; 0 0 0];
B = [0 -1 0]';

%%% Constants
VM = 3000; % ft/s
VT = 1000; % ft/s

VC = VM + VT; % ft/s
R0 = 40000; % ft

y0 = 0;

%%% Scenario top level parameters
Nprimes  = [3, 4, 5];
Ks = [1, 10, 100];


%%% Scenario bottom level parameters
Theta_HEs = [30, 0, 30]; % deg
nTs = [0, 3, 3].*32; % Gs

%%% Propagation time
t0 = 0;
dt = 0.1;
tf = R0/VC;
time = t0:dt:tf;

%%% Propagation time
cols = [colors.std.red; colors.std.blue; colors.std.grn];

% -------------------------------------------------
% ProNav
% -------------------------------------------------
if run_ProNav == 1
%%% Preallocating
X_PN{1} = zeros(length(time),3);
X_PN{2} = zeros(length(time),3);
X_PN{3} = zeros(length(time),3);

nC_PN{1} = zeros(length(time),1);
nC_PN{2} = zeros(length(time),1);
nC_PN{3} = zeros(length(time),1);

%%% Looping through N-prime valus
for Ni = 1:3
    %%% Preallocating and setting ICs
%     warning('Not sure which of these')
    
    %%% Looping through parameters
    for pii = 1:3
%         ydot0 = sind(Theta_HEs(pii))*VM;
        ydot0 = tand(Theta_HEs(pii))*VC;
        X_PN{Ni}(1,:) = [y0, ydot0, nTs(pii)];

        %%% Propagating
        for ti = 1:length(time)-1
            %%% Linearized ZEM
            tgo = tf - time(ti);
    %         R = VC * tgo;
            ZEM = X_PN{Ni}(ti,1) + X_PN{Ni}(ti,2)*tgo;


            %%% Calculating nC
            nC_PN{Ni}(ti) = Nprimes(Ni)*ZEM/(tgo^2);

            [X_p] = rk4_proNavGuidance(X_PN{Ni}(ti,:)',dt,nC_PN{Ni}(ti));

            %%% Storing
            X_PN{Ni}(ti+1,:) = X_p;
        end
        
        figure(pii)
        subplot(1,2,1); hold all
        plot(time,X_PN{Ni}(:,1),'linewidth',2,'color',cols(Ni,:));
        PlotBoi2('Time, sec','y',16,'LaTex')
%         title(sprintf('ProNav, N'' = %1.0f',Nprimes(Ni)));
        title(sprintf('\\theta_H_E = %1.0f°, n_T = %1.0f Gs',Theta_HEs(pii),nTs(pii)/32));
        
        subplot(1,2,2); hold all
        plot(time,nC_PN{Ni}(:),'linewidth',2,'color',cols(Ni,:));
        PlotBoi2('Time, sec','$n_c$',16,'LaTex')
    end
end

for pi = 1:3
    figure(pi)
%     legend(,sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(2),nTs(2)),sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(3),nTs(3)))
    legend(sprintf('ProNav, N'' = %1.0f',Nprimes(1)),sprintf('ProNav, N'' = %1.0f',Nprimes(2)),sprintf('ProNav, N'' = %1.0f',Nprimes(3)))
end

end % run_ProNav
% -------------------------------------------------
% Augmented ProNav
% -------------------------------------------------
if run_AugmentedProNav == 1
%%% Preallocating
X_APN{1} = zeros(length(time),3);
X_APN{2} = zeros(length(time),3);
X_APN{3} = zeros(length(time),3);

nC_APN{1} = zeros(length(time),1);
nC_APN{2} = zeros(length(time),1);
nC_APN{3} = zeros(length(time),1);

%%% Looping through N-prime valus
for Ni = 1:3
    %%% Preallocating and setting ICs
%     warning('Not sure which of these')
    
    %%% Looping through parameters
    for pii = 1:3
%         ydot0 = sind(Theta_HEs(pii))*VM;
        ydot0 = tand(Theta_HEs(pii))*VC;
        X_APN{Ni}(1,:) = [y0, ydot0, nTs(pii)];

        %%% Propagating
        for ti = 1:length(time)-1
            %%% Linearized ZEM
            tgo = tf - time(ti);
    %         R = VC * tgo;
            ZEM = X_APN{Ni}(ti,1) + X_APN{Ni}(ti,2)*tgo + 0.5*X_APN{Ni}(ti,3)*tgo*tgo;


            %%% Calculating nC
            nC_APN{Ni}(ti) = Nprimes(Ni)*ZEM/(tgo^2);

            [X_p] = rk4_proNavGuidance(X_APN{Ni}(ti,:)',dt,nC_APN{Ni}(ti));

            %%% Storing
            X_APN{Ni}(ti+1,:) = X_p;
        end
        
        figure(pii+3)
        subplot(1,2,1); hold all
        plot(time,X_APN{Ni}(:,1),'linewidth',2,'color',cols(Ni,:));
        PlotBoi2('Time, sec','y',16,'LaTex')
%         title(sprintf('Augmented ProNav, N'' = %1.0f',Nprimes(Ni)));
        title(sprintf('\\theta_H_E = %1.0f°, n_T = %1.0f Gs',Theta_HEs(pii),nTs(pii)/32));
        
        subplot(1,2,2); hold all
        plot(time,nC_APN{Ni}(:),'linewidth',2,'color',cols(Ni,:));
        PlotBoi2('Time, sec','$n_c$',16,'LaTex')
    end
end

for Ni = 1:3
    figure(Ni+3)
%     legend(sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(1),nTs(1)),sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(2),nTs(2)),sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(3),nTs(3)))
    legend(sprintf('Augmented ProNav, N'' = %1.0f',Nprimes(1)),sprintf('Augmented ProNav, N'' = %1.0f',Nprimes(2)),sprintf('Augmented ProNav, N'' = %1.0f',Nprimes(3)))
end



end % run_AugmentedProNav

% -------------------------------------------------
% Soft-Constrained Optimal Control
% -------------------------------------------------
if run_SoftConstraintOptimal == 1
%%% Preallocating
X_SCO{1} = zeros(length(time),3);
X_SCO{2} = zeros(length(time),3);
X_SCO{3} = zeros(length(time),3);

nC_SCO{1} = zeros(length(time),1);
nC_SCO{2} = zeros(length(time),1);
nC_SCO{3} = zeros(length(time),1);

%%% Looping through N-prime valus
for Ki = 1:3
    %%% Preallocating and setting ICs
%     warning('Not sure which of these')
    
    %%% Looping through parameters
    for pii = 1:3
%         ydot0 = sind(Theta_HEs(pii))*VM;
        ydot0 = tand(Theta_HEs(pii))*VC;
        X_SCO{Ki}(1,:) = [y0, ydot0, nTs(pii)];

        %%% Propagating
        for ti = 1:length(time)-1
            %%% Linearized ZEM
            tgo = tf - time(ti);
    %         R = VC * tgo;
%             ZEM = X_SCO{Ki}(ti,1) + X_SCO{Ki}(ti,2)*tgo + 0.5*X_SCO{Ki}(ti,3)*tgo*tgo;


            %%% Calculating nC
            nC_SCO{Ki}(ti) = 3*(X_SCO{Ki}(ti,1) + X_SCO{Ki}(ti,2)*tgo + 0.5*X_SCO{Ki}(ti,3)*tgo*tgo)*tgo/(3/Ks(Ki) + tgo^3);

            [X_p] = rk4_proNavGuidance(X_SCO{Ki}(ti,:)',dt,nC_SCO{Ki}(ti));

            %%% Storing
            X_SCO{Ki}(ti+1,:) = X_p;
        end
        
        figure(pii+6)
        subplot(1,2,1); hold all
        plot(time,X_SCO{Ki}(:,1),'linewidth',2,'color',cols(Ki,:));
        PlotBoi2('Time, sec','y',16,'LaTex')
%         title(sprintf('Soft-Constraint Optimal Control, K = %1.0f',Ks(Ki)));
        title(sprintf('\\theta_H_E = %1.0f°, n_T = %1.0f Gs',Theta_HEs(pii),nTs(pii)/32));
        
        subplot(1,2,2); hold all
        plot(time,nC_SCO{Ki}(:),'linewidth',2,'color',cols(Ki,:));
        PlotBoi2('Time, sec','$n_c$',16,'LaTex')
    end
end

for Ki = 1:3
    figure(Ki+6)
%     legend(sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(1),nTs(1)),sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(2),nTs(2)),sprintf('\\theta_HE = %1.0f°, n_T = %1.0f Gs',Theta_HEs(3),nTs(3)))
    legend(sprintf('Soft-Constraint Optimal Control, K = %1.0f',Ks(1)),sprintf('Soft-Constraint Optimal Control, K = %1.0f',Ks(2)),sprintf('Soft-Constraint Optimal Control, K = %1.0f',Ks(3)))
end


end % run_SoftConstraintOptimal

% ========================================================================
%%% P1-b
% ========================================================================
if run_AugmentedProNav_1b == 1
%%% Preallocating
X_APNb = zeros(length(time),3);
nC_APNb = zeros(length(time),1);

%%% Constants
Theta_HE = 30;
nT = 3*32;
Nprime = 4;
sigma = sqrt(0.3);

%%% Preallocating and setting ICs
% warning('Not sure which of these')

%%% Looping through parameters
% ydot0 = sind(Theta_HE)*VM;
ydot0 = tand(Theta_HE)*VC;
    X_APNb(1,:) = [y0, ydot0, nT];

%%% Propagating
engineCutoff = 0;
for ti = 1:length(time)-1
    %%% Linearized ZEM
    tgo = tf - time(ti) + randn(1,1)*sigma;
%     tgo = tf - time(ti) + randn(1,1)*sigma;
%     tgo = tf - time(ti);
%     R = VC * tgo;
    ZEM = X_APNb(ti,1) + X_APNb(ti,2)*tgo + 0.5*X_APNb(ti,3)*tgo*tgo;


    %%% Calculating nC
    nC_APNb(ti) = Nprime*ZEM/(tgo^2);
%     if time(ti) > 8
%         nC_APNb(ti) = 0;
%     end
%     if (X_APNb(ti,1) < 0.01 && time(ti) > 4) || engineCutoff == 1
%         engineCutoff = 1;
%         nC_APNb(ti) = 0;
%     end
    
    [X_p] = rk4_proNavGuidance(X_APNb(ti,:)',dt,nC_APNb(ti));

    %%% Storing
    X_APNb(ti+1,:) = X_p;
end

figure(10)
subplot(1,2,1); hold all
plot(time,X_APNb(:,1),'linewidth',2,'color',cols(1,:));
PlotBoi2('Time, sec','y',16,'LaTex')
title(sprintf('Augmented ProNav, N'' = %1.0f',Nprime));

subplot(1,2,2); hold all
plot(time,nC_APNb(:),'linewidth',2,'color',cols(1,:));
PlotBoi2('Time, sec','$n_c$',16,'LaTex')
legend(sprintf('\\theta_H_E = %1.0f°, n_T = %1.0f Gs',Theta_HE,nT/32))



end % run_AugmentedProNav_1b









% ========================================================================
% ========================================================================
%%% Functions
% ========================================================================
% ========================================================================

function [X_p] = rk4_proNavGuidance(X,dt,nC)
%%%
%%% Inputs:
%           1) X  - State, [3x1]
%           2) dt - step size, [1x1]
%           3) nC
%           4) 
%           5) 
%           6)
%%% Outputs:
%           1) X_p - State at time t+dt, [3x1]
% ========================================================================
%%% EOM
dXdt = @(y, ydot, nT, nC) [ydot; nT - nC; 0];     
%                        dXdt(X(1:3),nC)
%%% Runge Kutta Time!
k1 = dXdt(X(1), X(2), X(3),nC);
X1 = X + k1.*(dt/2);
k2 = dXdt(X1(1),X1(2),X1(3),nC);
X2 = X + k2.*(dt/2);
k3 = dXdt(X2(1),X2(2),X2(3),nC);
X3 = X + k3.*dt;
k4 = dXdt(X3(1),X3(2),X3(3),nC);

X_p = X + (k1 + 2*k2 + 2*k3 + k4)*(dt/6);

end










