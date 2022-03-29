% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 06/15/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath    = '~/CU_Google_Drive/Documents/MatGit/mbin';
% dataBinpath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs/ShallowImpactDataSets';
dataBinpath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs/ManifoldShallowImpactDataSets';
addpath(genpath(mbinPath))
ticWhole = tic;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Periodic orbit ICs
PO_ICs = get_PO_ICs();

% ========================================================================
%%% Run Switches
% ========================================================================



% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Choose data file
% -------------------------------------------------
dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_Lyapunov.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_Vertical.nodes2000.txt'; % n_shallowImpacts = 197;
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_SHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_NHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_All.nodes2000.txt';

%%% Create path to dataFile
dataFilePath = [dataBinpath, '/',dataFile];

% -------------------------------------------------
%%% Set up system
% -------------------------------------------------
%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(dataFile, bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% Setting parameter structure
prms.u     = secondary.MR;
prms.R1    = primary.R / rNorm;
prms.R2    = secondary.R_n;
prms.rNorm = rNorm;
prms.n = 1;

% -------------------------------------------------
%%% Load PO family data
% -------------------------------------------------
%%% Load shallow angle impact data
shallowAngleImpactData = dlmread(dataFilePath,',',1,0);

%%% Column specifiers
c_x0              = 1;
c_y0              = 2;
c_z0              = 3;
c_xd0             = 4;
c_yd0             = 5;
c_zd0             = 6;
c_Tp              = 7;
c_lat               = 8;
c_lon               = 9;
c_impactAngle_deg   = 10;
c_impactHeading_lat = 11;
c_impactHeading_lon = 12;
c_JC                = 13;
c_landingVel_mps    = 14;

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);
options_impact           = odeset('Events', @event_Impact_CR3Bn, 'RelTol',tol,'AbsTol',tol);

stm0_colVec = reshape(eye(6),36,1);

% ========================================================================
%%% Run
% ========================================================================

% index = 200;
% [T_PO, X_PO] = ode113(@Int_CR3Bn, [0, shallowAngleImpactData(index,c_Tp)], [shallowAngleImpactData(index,c_x0:c_zd0)'], options, prms);
% [T_PO_impact, X_PO_impact] = ode113(@Int_CR3Bn, [0, shallowAngleImpactData(index,c_Tp)].*2, [shallowAngleImpactData(index,c_x0:c_zd0)'], options_impact, prms);
% 
% 
% figure; hold all
% plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3), 'linewidth', 1, 'color', colors.black)
% plot3(X_PO_impact(:,1),X_PO_impact(:,2),X_PO_impact(:,3), 'linewidth', 1, 'color', colors.mag)
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(26)
% view(0,90)

% ========================================================================
%%% Loop through states, add burns, re-integrate, and plot
% ========================================================================
burnMagnitude_mps = 10;

minutesBeforeLandingForFinalBurn = 60; % 60


% impactLatLons_plus  = [];
% impactLatLons_minus = [];
989
shallowAngleImpactData = shallowAngleImpactData(1:6:6273, :);
impactLatLons_plus  = cell(size(shallowAngleImpactData,1),1);
impactLatLons_minus = cell(size(shallowAngleImpactData,1),1);

parfor solution_i = 1:size(shallowAngleImpactData,1)
% for solution_i = 1:20
    [T_PO, X_PO] = ode113(@Int_CR3Bn, [0, shallowAngleImpactData(solution_i,c_Tp)], [shallowAngleImpactData(solution_i,c_x0:c_zd0)'], options, prms);

    stateIndex_firstBurn = find(abs((T_PO - T_PO(end)/3)) == min(abs(T_PO - T_PO(end)/3)));


    tn_finalBurn = shallowAngleImpactData(solution_i,c_Tp) - (minutesBeforeLandingForFinalBurn*60)/tNorm;

    stateIndex_finalBurn = find(abs((T_PO - tn_finalBurn)) == min(abs(T_PO - tn_finalBurn)));

    T_end = shallowAngleImpactData(solution_i,c_Tp).*1.2;

%     impactLatLons_minus = [];
%     impactLatLons_plus  = [];
% %     % for kk = 1:size(X_PO,1)
% %     % for kk = 100:302 
    for kk = stateIndex_firstBurn:stateIndex_finalBurn 
        newState = X_PO(kk,1:6)';

        unitVec_vel = newState(4:6)./norm(newState(4:6));

        newState_plusVel = newState + [0; 0; 0; unitVec_vel].*(burnMagnitude_mps/(vNorm *1000));
        newState_minusVel = newState - [0; 0; 0; unitVec_vel].*(burnMagnitude_mps/(vNorm *1000));

        T_int_kk = T_end - T_PO(kk);

        [T_kk_plus_impact, X_kk_plus_impact]   = ode113(@Int_CR3Bn, [0, T_int_kk], newState_plusVel, options_impact, prms);
        [T_kk_minus_impact, X_kk_minus_impact] = ode113(@Int_CR3Bn, [0, T_int_kk], newState_minusVel, options_impact, prms);

        if T_kk_plus_impact(end) ~= T_int_kk % There was an impact
            [latLon_kk_plus_deg, impactAngle_kk_plus_deg, latLonHeadingHat_kk_plus, impactColor_kk_plus] = getImpactConditions(X_kk_plus_impact(end,1:6), prms, [0, 3, 90], [colors.blue2; colors.black]);
            if impactAngle_kk_plus_deg <= 3
%                 plot3(X_kk_plus_impact(:,1),X_kk_plus_impact(:,2),X_kk_plus_impact(:,3), 'linewidth', 1, 'color', colors.blue2)
%                 impactLatLons_plus = [impactLatLons_plus; latLon_kk_plus_deg];
                impactLatLons_plus{solution_i} = latLon_kk_plus_deg;
            end
        end

        if T_kk_minus_impact(end) ~= T_int_kk % There was an impact
            [latLon_kk_minus_deg, impactAngle_kk_minus_deg, latLonHeadingHat_kk_minus, impactColor_kk_minus] = getImpactConditions(X_kk_minus_impact(end,1:6), prms, [0, 3, 90], [colors.blue2; colors.black]);
            if impactAngle_kk_minus_deg <= 3
%                 plot3(X_kk_minus_impact(:,1),X_kk_minus_impact(:,2),X_kk_minus_impact(:,3), 'linewidth', 1, 'color', colors.red2)
%                 impactLatLons_minus = [impactLatLons_minus; latLon_kk_minus_deg];
                impactLatLons_minus{solution_i} = latLon_kk_minus_deg;
            end
        end


    end

%     plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3), 'linewidth', 1.5, 'color', colors.black)

    
end

% figure; hold all
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% xlim([-180 180])
% ylim([-90 90])
% 
% 
% if ~isempty(impactLatLons_plus)
%     plot(impactLatLons_plus(:,2), impactLatLons_plus(:,1), 'x', 'markersize', 4, 'linewidth', 2, 'color', colors.blue2)
% end
% if ~isempty(impactLatLons_minus)
%     plot(impactLatLons_minus(:,2), impactLatLons_minus(:,1), 'x', 'markersize', 4, 'linewidth', 2, 'color', colors.red2)
% end
% 
% plot(shallowAngleImpactData(:,c_lon), shallowAngleImpactData(:,c_lat), '.', 'markersize', 3, 'color', colors.black)

    

figure; hold all
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
xlim([-180 180])
ylim([-90 90])
for kk = 1:size(shallowAngleImpactData,1)
    if ~isempty(impactLatLons_plus{kk})
        plot(impactLatLons_plus{kk}(:,2), impactLatLons_plus{kk}(:,1), 'x', 'markersize', 4, 'linewidth', 2, 'color', colors.blue2)
        plot(impactLatLons_plus{kk}(:,2), -impactLatLons_plus{kk}(:,1), 'x', 'markersize', 4, 'linewidth', 2, 'color', colors.blue2)
    end
    if ~isempty(impactLatLons_minus{kk})
        plot(impactLatLons_minus{kk}(:,2), impactLatLons_minus{kk}(:,1), 'x', 'markersize', 4, 'linewidth', 2, 'color', colors.red2)
        plot(impactLatLons_minus{kk}(:,2), -impactLatLons_minus{kk}(:,1), 'x', 'markersize', 4, 'linewidth', 2, 'color', colors.red2)
    end
end
plot(shallowAngleImpactData(:,c_lon), shallowAngleImpactData(:,c_lat), '.', 'markersize', 3, 'color', colors.black)
plot(shallowAngleImpactData(:,c_lon), -shallowAngleImpactData(:,c_lat), '.', 'markersize', 3, 'color', colors.black)


% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)






