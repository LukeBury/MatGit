clear
clc
% close all
tic


% -------------------------------------------------
% Choosing data files
% -------------------------------------------------
for kk = [1:7]
    if kk == 1
        %%% 50 mps
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_50mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_50mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_50mps_50km_149v0s_land.txt';
    end
    if kk == 2
        %%% 100 mps
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_100mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_100mps_50km_149v0s_land.txt';
    end
    if kk == 3
        %%% 150 mps
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_150mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_150mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_150mps_50km_149v0s_land.txt';
    end
    if kk == 4
        %%% 200 mps
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_200mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_200mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_200mps_50km_149v0s_land.txt';
    end
    if kk == 5
        %%% 250 mps
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_250mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_250mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_250mps_50km_149v0s_land.txt';
    end
    if kk == 6
        %%% 300 mps
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_300mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_300mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_300mps_50km_149v0s_land.txt';
    end
    if kk == 7
        %%% 350 mps
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_350mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_350mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_350mps_50km_149v0s_land.txt';
    end
    
    if kk == 8
        %%% 50 mps - ZH                                                                      
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_50mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_50mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_50mps_50km_149v0s_land.txt';
    end
    
    if kk == 9
        %%% 100 mps - ZH
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_100mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_100mps_50km_149v0s_land.txt';
    end
    
    if kk == 10
        %%% 150 mps - ZH
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_150mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_150mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_150mps_50km_149v0s_land.txt';
    end
    
    if kk == 11
        %%% 200 mps - ZH
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_200mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_200mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_200mps_50km_149v0s_land.txt';
    end
    
    if kk == 12
        %%% 250 mps - ZH
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_250mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_250mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_250mps_50km_149v0s_land.txt';
    end
    
    if kk == 13
        %%% 300 mps - ZH
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_300mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_300mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_300mps_50km_149v0s_land.txt';
    end
    
    if kk == 14
        %%% 350 mps - ZH
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_350mps_50km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_350mps_50km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_J21_350mps_50km_149v0s_land.txt';
    end
    
    
    if kk == 15
        %%% 13 mps (Enc)
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_13mps_4km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_13mps_4km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_13mps_4km_149v0s_land.txt';
    end
    
    if kk == 16
        %%% 26 mps (Enc)
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_26mps_4km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_26mps_4km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_26mps_4km_149v0s_land.txt';
    end
    
    if kk == 17
        %%% 39 mps (Enc)
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_39mps_4km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_39mps_4km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_39mps_4km_149v0s_land.txt';
    end
    
    if kk == 18
        %%% 52 mps (Enc)
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_52mps_4km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_52mps_4km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_52mps_4km_149v0s_land.txt';
    end
    
    if kk == 19
        %%% 65 mps (Enc)
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_65mps_4km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_65mps_4km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_65mps_4km_149v0s_land.txt';
    end
    
    if kk == 20
        %%% 78 mps (Enc)
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_78mps_4km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_78mps_4km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_78mps_4km_149v0s_land.txt';
    end
    
    if kk == 21
        %%% 91 mps (Enc)
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_91mps_4km_149v0s_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_91mps_4km_149v0s_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_encL2_91mps_4km_149v0s_land.txt';
    end
    
    
    %%% Test case
    if kk == 989
        logFile            = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/testFile_log.txt';
        impactFile         = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/testFile_data.txt';
        lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/testFile_land.txt';
    end

    




% -------------------------------------------------
% Running
% -------------------------------------------------
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
bodies = getBodyData(mbinPath);

%%% Was this file run with zonal harmonics?
if contains(logFile,'J21') == 0
    on_J21 = 0;
elseif contains(logFile,'J21') == 1
    on_J21 = 1;
end

iGS_LoadPlot2(logFile,impactFile,lowImpactAngleFile,mbinPath,on_J21)

end
toc





% ========================================================================
%%% Load/Plot Function
% ========================================================================
function iGS_LoadPlot2(logFile,impactFile,lowImpactAngleFile,mbinPath,on_J21)
% function iGS_LoadPlot2(logFile,impactFile,lowImpactAngleFile,primary,secondary,bins_impactAngles,bins_neckSections,JC_scInitial,Lpoint,y_neck_upper,z_neck_upper)
%%% For numerical integration in the normalized CR3BP with J2 of each body
%%% Inputs:
%          filePath -
%          logFile - 
%          impactFile -
%          lowImpactAngleFile
%          primary - 
%          secondary - 
%          L123 - 
%          bins_neckSections
%          bins_impactAngles
%          yzContourPoints
%          
%          
%          
%          
% ========================================================================
%%% Settings and Run Switches
% ========================================================================
% -------------------------------------------------
% Run Switches
% -------------------------------------------------
on_symmetricImpactResults = 0;

plot_fullLowTrajectories  = 0;
plot_lowTrajectories_r0   = 0;
plot_allTrajectories_r0   = 0;
plot_fullSystemContour    = 0;
plot_sectionColors        = 0;
plot_allTrajectories      = 0; % Careful now!!

plot_binImpactAngles      = 0;
plot_binNeckSections_rad  = 0;
plot_binNeckSections_vert = 0;
plot_binNeckSections_horz = 0;
plot_metaData             = 0;
plot_latitudeV0AngleCorr  = 0;
plot_binTOF               = 0;
calc_maxLatitudeAnalysis  = 1;

% ========================================================================
%%% Setting data
% ========================================================================
% -------------------------------------------------
% Importing data structures
% -------------------------------------------------
%%% Color options/schemes
colors = get_colors();

%%% General data on solar system bodies
bodies = getBodyData(mbinPath);
% -------------------------------------------------
% Loading data from log file
% -------------------------------------------------
flogFile = fopen(logFile);
logData = textscan(flogFile, '%s %s','delimiter',':');
fclose(flogFile);

for kk = 1:size(logData{1})
if isequal(logData{1}{kk},'total trajs') == 1
nTraj = str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'primary') == 1
primary = getfield(bodies,logData{2}{kk});
end

if isequal(logData{1}{kk},'secondary') == 1
secondary = getfield(bodies,logData{2}{kk});
end

if isequal(logData{1}{kk},'bins_impactAngles') == 1
    bins_impactAngles = str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'bins_neckSectionScalars') == 1
    bins_neckSectionScalars =str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'dvLp_mps') == 1
    dvLp_mps =str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'JC_scInitial') == 1
    JC_scInitial = str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'Lpoint') == 1
    Lpoint = str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'y_neck_upper') == 1
    y_neck_upper = str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'z_neck_upper') == 1
    z_neck_upper = str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'maxLat') == 1
    maxLat = str2num(logData{2}{kk});
end

if isequal(logData{1}{kk},'maxLowLat') == 1
    maxLowLat = str2num(logData{2}{kk});
end
end

% -------------------------------------------------
% Getting constants
% -------------------------------------------------
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Find equilibrium points
if on_J21 == 0
    L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
elseif on_J21 == 1
    J21 = primary.J2;
    J22 = 0;
    R1_n = primary.R/rNorm;
    R2_n = secondary.R_n;
    L123 = EquilibriumPoints_J2(secondary.MR,J21,J22,R1_n,R2_n,1:3);
end



% ========================================================================
%%% Loading General Impact Data
% ========================================================================
%%% Load impact data
impactData = dlmread(impactFile,',',1,0);
nImpactTraj = size(impactData,1);
%%% set column specifiers
c_full_impactBin       = 2;
c_full_neckBin         = 3;
c_full_impactLatitude  = 4;
c_full_impactLongitude = 5;
c_full_impactAng       = 6;
c_full_endTime         = 7;
c_full_y0_n            = 8;
c_full_z0_n            = 9;
c_full_v0_az           = 10;
c_full_v0_el           = 11;
c_full_maxLat          = 12;


%%% Bin counts
binCount_ImpactAngles = length(bins_impactAngles) - 1;
binCount_NeckSections = length(bins_neckSectionScalars);

%%% Preallocate
binData_impactAngles(binCount_ImpactAngles).latLons = [];
binData_neckSections(length(bins_neckSectionScalars)).latLons = [];

%%% Store impact angle bins
for kk = 1:binCount_ImpactAngles
    %%% Indices of all impact trajectories that belong to the current bin
    currentBinIndices = find(impactData(:,c_full_impactBin)==kk);
    
    %%% Grabbing those trajectories and storing their impact lat/lon
    binData_impactAngles(kk).latLons = impactData(currentBinIndices,[c_full_impactLatitude,c_full_impactLongitude]);
end

%%% Store neck section bins
for kk = 1:binCount_NeckSections
    %%% Indices of all impact trajectories that belong to the current bin
    currentBinIndices = find(impactData(:,c_full_neckBin)==kk);
    
    %%% Grabbing those trajectories and storing their impact lat/lon
    binData_neckSections(kk).latLons = impactData(currentBinIndices,[c_full_impactLatitude,c_full_impactLongitude]);
end

if plot_allTrajectories == 1
    allTrajs_r = [];
    for kk = 1:size(impactData,1)
        v0_n = azEl_2_Vec( [1; 0; 0], impactData(kk,c_full_v0_az)*180/pi, impactData(c_full_v0_el)*180/pi );
        r0_kk = [L123(2,1); impactData(kk,c_full_y0_n); impactData(kk,c_full_z0_n)];
        if on_J21 == 0
            JC_r0 = JacobiConstantCalculator(secondary.MR,r0_kk(1:3)' ,[0, 0, 0]);
        elseif on_J21 == 1
            JC_r0 = JacobiConstantCalculator_J2(secondary.MR,r0_kk(1:3),[0,0,0], R1_n, R2_n, J21, J22);
        end
        v0_n = v0_n .* sqrt(JC_r0 - JC_scInitial);
        X0_kk = [r0_kk; v0_n];
        endTime_kk = impactData(kk,c_full_endTime);
        %%% Setting time vector
        t_i = 0; % sec
        t_f = endTime_kk; % Long bc events are watching for impact or escape
        n_dt = 10000;
        time0_n = linspace(t_i,t_f,n_dt);

        %%% Choosing ode tolerance
        tol = 1e-13;

        %%% Setting integrator options
        options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

        %%% Setting necessary parameters for integration
        prms.u = secondary.MR;
        prms.R2_n = secondary.R_n;
        prms.L1x = L123(1,1);
        prms.L2x = L123(2,1);
        if on_J21 == 1
            prms.R1_n = primary.R/rNorm;
            prms.J21  = primary.J2;
        end
        %%% Propagating trajectory
        if on_J21 == 0
            [time_n, X_BCR_n, ~, ~, ~] = ode113(@Int_CR3Bn,...
                time0_n, X0_kk, options_ImpactEscape, prms);
        elseif on_J21 == 1
            [time_n, X_BCR_n, ~, ~, ~] = ode113(@Int_CR3Bn_ZH,...
                time0_n, X0_kk, options_ImpactEscape, prms);
        end

        allTrajs_r = [allTrajs_r; NaN(1,3); X_BCR_n(:,1:3)];
    end
    figure; hold all
    plot3(allTrajs_r(:,1),allTrajs_r(:,2),allTrajs_r(:,3),'b')
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
end
% ========================================================================
%%% Loading Low-Impact-Angle Data
% ========================================================================
% -------------------------------------------------
% Load Data
% -------------------------------------------------
%%% Load low impact data
lowImpactMat = dlmread(lowImpactAngleFile,',',1,0);
nLowImpactTraj = size(lowImpactMat,1);


if size(lowImpactMat) ~= [1,1] % if there are any low-angle impacts:
    
    % -------------------------------------------------
    % Propagate low-impact-angle trajectories
    % -------------------------------------------------
    if plot_fullLowTrajectories == 1
    %%% Preallocating
    lowTrajs{nLowImpactTraj} = [];

    %%% Setting time vector
    t_i = 0; % sec
    t_f = 4*pi; % Long bc events are watching for impact or escape
    dt = t_f/10000;
    time0_n = t_i:dt:t_f;

    %%% Choosing ode tolerance
    tol = 1e-12;

    %%% Setting integrator options
    options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

    %%% Setting necessary parameters for integration
    prms.u = secondary.MR;
    prms.R2_n = secondary.R_n;
    prms.L1x = L123(1,1);
    prms.L2x = L123(2,1);
    if on_J21 == 1
        prms.R1_n = primary.R/rNorm;
        prms.J21  = primary.J2;
    end

    for kk = 1:nLowImpactTraj
        %%% Propagating trajectory
        if on_J21 == 0
            [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
                time0_n, lowImpactMat(kk,2:7)', options_ImpactEscape, prms);
        elseif on_J21 == 1
            [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn_ZH,...
                time0_n, lowImpactMat(kk,2:7)', options_ImpactEscape, prms);
        end

        lowTrajs{kk} = [X_BCR_n, time_n];
    end
    
    %%% Plotting full low-landing-angle trajectories
    if plot_fullLowTrajectories == 1
        figure; hold all
        plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
        PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
        axis equal
        for kk = 1:length(lowTrajs)
            plot3(lowTrajs{kk}(:,1),lowTrajs{kk}(:,2),lowTrajs{kk}(:,3),'b')
        end
    end
    end % plot_fullLowTrajectories == 1 
    %%% Plotting initial positions of low-landing-angle trajectories

    
end % size(lowImpactMat) ~= [1,1]
% ========================================================================
%%% Determining necessary data for plotting
% ========================================================================

% -------------------------------------------------
% Find Y-Z contour points and making neck section bins
% -------------------------------------------------
%%% Create grid of starting locations based on y-z neck
ys = linspace(-y_neck_upper*2, y_neck_upper*2, 1000);
zs = linspace(-z_neck_upper*2, z_neck_upper*2, 1000);
[Y_yz,Z_yz] = meshgrid(ys,zs);

%%% Calculating JCs across y-z grid
JCs_yz_Lpoint = zeros(size(Y_yz));
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Y_yz,2)
        %%% Zero-Velocity Curve
        if on_J21 == 0
            zv = JacobiConstantCalculator(secondary.MR,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        elseif on_J21 == 1
            zv = JacobiConstantCalculator_J2(secondary.MR,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)],[0,0,0], R1_n, R2_n, J21, J22);
        end
        
        JCs_yz_Lpoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzContourPoints4 ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC_scInitial );

yzContourPoints3 = yzContourPoints4.*bins_neckSectionScalars(3);
yzContourPoints2 = yzContourPoints4.*bins_neckSectionScalars(2);
yzContourPoints1 = yzContourPoints4.*bins_neckSectionScalars(1);





% ========================================================================
%%% Plotting Data
% ========================================================================
% -------------------------------------------------
% Initial conditions in y-z plane of low-angle impactors along with bounding contours
% -------------------------------------------------
if plot_lowTrajectories_r0 == 1
    figure; hold all
%     plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
    axis equal
    if isequal(lowImpactMat,0) == 0
        lowTrajX0s = lowImpactMat(:,2:7);
        plot3(lowTrajX0s(:,1),lowTrajX0s(:,2),lowTrajX0s(:,3),'b.','markersize',6)
        plot3(lowTrajX0s(:,1),lowTrajX0s(:,2),-lowTrajX0s(:,3),'b.','markersize',6)
    end
%     yz4_indices = find(yzContourPoints4(2,:) > 0);
%     plot3(ones(size(yzContourPoints4,2),1).*L123(Lpoint,1), yzContourPoints4(1,yz4_indices), yzContourPoints4(2,yz4_indices),'k','linewidth',3)
    plot3(ones(size(yzContourPoints4,2),1).*L123(Lpoint,1), yzContourPoints4(1,:),yzContourPoints4(2,:),'k','linewidth',3)
    scale = 1.1;
    ylim([-scale*max(abs(yzContourPoints4(1,:))) scale*max(abs(yzContourPoints4(1,:)))]);
    zlim([-scale*max(abs(yzContourPoints4(2,:))) scale*max(abs(yzContourPoints4(2,:)))]);
    view(90,0)
end

% -------------------------------------------------
% Initial conditions in y-z plane along with bounding contours
% -------------------------------------------------
if plot_allTrajectories_r0 == 1
%     %%% Creating ZV contours
%     % Contour bounds
%     xCont_min = L123(1,1)-2*secondary.R_n;
%     xCont_max = L123(2,1)+2*secondary.R_n;
%     yCont_min = -15*secondary.R_n;
%     yCont_max = 15*secondary.R_n;
% 
%     % Creating x-y grid
%     xs = linspace(xCont_min,xCont_max,600);
%     ys = linspace(yCont_min,yCont_max,200);
%     [X_xy, Y_xy] = meshgrid(xs,ys);
%     clear xs ys zs
% 
%     % Calculating JCs across x-y grid
%     JCs_xy = zeros(size(X_xy));
%     for xk = 1:size(X_xy,1)
%         for yk = 1:size(X_xy,2)
%             %%% Zero-Velocity Curve
%             if on_J21 == 0
%                 zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0], [0, 0, 0]);
%             elseif on_J21 == 1
%                 zv = JacobiConstantCalculator_J2(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0], [0, 0, 0], R1_n, R2_n, J21, J22);
%             end
% 
%             JCs_xy(xk,yk) = zv;
%         end
%     end

    figure; hold all
%     yz4_indices = find(yzContourPoints4(2,:) > 0);
%     plot3(ones(size(yzContourPoints4,2),1).*L123(Lpoint,1), yzContourPoints4(1,yz4_indices), yzContourPoints4(2,yz4_indices),'k','linewidth',3)
    plot3(ones(size(yzContourPoints4,2),1).*L123(Lpoint,1), yzContourPoints4(1,:),yzContourPoints4(2,:),'k','linewidth',3)
%     plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
%     [xyContourPoints,href] = contour(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],...
%         'color',colors.std.black,'linewidth',3);
    plot3(ones(size(impactData,1),1).*L123(2,1),impactData(:,c_full_y0_n),impactData(:,c_full_z0_n),'b.','markersize',6)
    plot3(ones(size(impactData,1),1).*L123(2,1),impactData(:,c_full_y0_n),-impactData(:,c_full_z0_n),'b.','markersize',6)
    PlotBoi3('x$_n$','y$_n$','z$_n$',16,'LaTex')
%     view(-35,24)
%     camva(9)
    axis equal
    scale = 1.1;
    ylim([-scale*max(abs(yzContourPoints4(1,:))) scale*max(abs(yzContourPoints4(1,:)))]);
    zlim([-scale*max(abs(yzContourPoints4(2,:))) scale*max(abs(yzContourPoints4(2,:)))]);
    view(90,0)
end

% -------------------------------------------------
% ZV contour of initial energy level across full CR3BP system
% -------------------------------------------------
if plot_fullSystemContour == 1
    %%% Creating ZV contours
    % Contour bounds
    xCont_min = L123(3,1)*1.2;
    xCont_max = L123(2,1)*1.2;
    yCont_min = -1.1;
    yCont_max = 1.1;

    % Creating x-y grid
    xs = linspace(xCont_min,xCont_max,750);
    ys = linspace(yCont_min,yCont_max,750);
    [X_xy, Y_xy] = meshgrid(xs,ys);
    clear xs ys

    % Calculating JCs across x-y grid
    JCs_xy = zeros(size(X_xy));
    for xk = 1:size(X_xy,1)
        for yk = 1:size(X_xy,2)
            %%% Zero-Velocity Curve
            if on_J21 == 0
                zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0], [0, 0, 0]);
            elseif on_J21 == 1
                zv = JacobiConstantCalculator_J2(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0], [0, 0, 0], R1_n, R2_n, J21, J22);
            end
            
            JCs_xy(xk,yk) = zv;
        end
    end

    figure; hold all
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
    plotBodyTexture3(primary.R./rNorm, [-secondary.MR, 0, 0], primary.img)
    [xyContourPoints,href] = contourf(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],...
        'color',colors.std.black,'linewidth',1.5);
    colormap(colors.sch.r6(6,:))
    PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
    view(0,90)
    camva(9)
    axis equal
end

% -------------------------------------------------
% Lat/Lon plot colored by initial neck section (radial)
% -------------------------------------------------
if plot_binNeckSections_rad == 1
    %%% Choosing color for bins of neck sections 
%     binColors_NeckSections = [colors.sch.d4_1(2,:);colors.sch.d4_1(1,:);...
%                              colors.sch.d4_1(4,:);colors.sch.d4_1(3,:)];
                         
    binColors_NeckSections = colorScale([colors.std.mag;colors.std.cyan],4);
                     
    figure; hold all
    for kk = 1:4
        if isempty(binData_neckSections(kk).latLons) == 0
            plot(binData_neckSections(kk).latLons(:,2),binData_neckSections(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_NeckSections(kk,:))
            if on_symmetricImpactResults == 1
                plot(binData_neckSections(kk).latLons(:,2),-binData_neckSections(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_NeckSections(kk,:))
            end
        end
    end
    PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
    xlim([-180 180])
    ylim([-90 90])

    %%% Colorbar for whole figure
    % a = gcf;
    % cbar = colorbar('Position',a.Position);
    cbar1 = colorbar;
    caxis([0 size(binColors_NeckSections,1)]);
    cbar1.FontName     = 'Arial';
    cbar1.FontSize     = 10;
    cbar1.Ticks        = [0.5, 1.5, 2.5, 3.5];
    % cbar1.TickLabels = num2cell(bins_neckSectionScalars(1:end-1).*rNorm);
    cbar1.TickLabels = num2cell([1 2 3 4]);
    % cbar1.Label.String = {sprintf('|r_0| from L_%1.0f, km',Lpoint)};
    cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar1.Label.Position = [.2, 4.2, 0];
    colormap(binColors_NeckSections)
    
    % -------------------------------------------------
    % Plot showing the initial neck section color coding
    % -------------------------------------------------
    if plot_sectionColors == 1
%         yz4_indices = find(yzContourPoints4(2,:) > 0);
%         yz3_indices = find(yzContourPoints3(2,:) > 0);
%         yz2_indices = find(yzContourPoints2(2,:) > 0);
%         yz1_indices = find(yzContourPoints1(2,:) > 0);

        figure; hold all
        plot(yzContourPoints4(1,:), yzContourPoints4(2,:),'k','linewidth',3)
        fill(yzContourPoints4(1,:),yzContourPoints4(2,:),binColors_NeckSections(4,:))
        fill(yzContourPoints3(1,:),yzContourPoints3(2,:),binColors_NeckSections(3,:))
        fill(yzContourPoints2(1,:),yzContourPoints2(2,:),binColors_NeckSections(2,:))
        fill(yzContourPoints1(1,:),yzContourPoints1(2,:),binColors_NeckSections(1,:))
        plot([-max(abs(yzContourPoints4(1,:))) max(abs(yzContourPoints4(1,:)))],[0 0],'k','linewidth',1)

        PlotBoi2('$y_n$','$z_n$',16,'LaTex')
        axis equal

    end
end
% -------------------------------------------------
% Lat/Lon plot colored by initial neck section (vertical)
% -------------------------------------------------

if plot_binNeckSections_vert == 1
    % -------------------
    % Binning data
    % -------------------
    %%% Number of desired vertical bins
    n_vertBins = 4; 
    
    %%% Choosing color for bins of neck sections 
    binColors_vertNeckSections = colorScale([colors.std.mag;colors.std.cyan],n_vertBins);
    
    %%% Find min and max bounds of bins
%     vert_min = min(abs(impactData(:,c_z0_n)));
%     vert_max = max(abs(impactData(:,c_z0_n)));
    vert_min = min(abs(yzContourPoints4(2,:)));
    vert_max = max(abs(yzContourPoints4(2,:)));
    
    %%% Assign bin values
    bins_vertNeckSections = linspace(0, vert_max, n_vertBins+1);
    
    %%% Preallocate bin data structures
    binData_vertNeckSections(n_vertBins).latLons = [];

    %%% Loop through each bin and store appropriate data
    for kk = 1:n_vertBins
        %%% Indices of all impact trajectories that belong to the current bin
        currentBinIndices = find( impactData(:,c_full_z0_n) >= bins_vertNeckSections(kk) &  impactData(:,c_full_z0_n) < bins_vertNeckSections(kk+1));
        
        %%% Grabbing those trajectories and storing their impact lat/lon
        binData_vertNeckSections(kk).latLons = impactData(currentBinIndices,[c_full_impactLatitude,c_full_impactLongitude]);
    end
    
    % -------------------
    % Plotting impact map
    % -------------------
    figure; hold all
    for kk = 1:n_vertBins
        if isempty(binData_vertNeckSections(kk).latLons) == 0
            plot(binData_vertNeckSections(kk).latLons(:,2),binData_vertNeckSections(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_vertNeckSections(kk,:))
            if on_symmetricImpactResults == 1
                plot(binData_vertNeckSections(kk).latLons(:,2),-binData_vertNeckSections(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_vertNeckSections(kk,:))
            end
        end
    end
    PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16, 'LaTex')
    xlim([-180 180])
    ylim([-90 90])

    %%% Colorbar for whole figure
    % a = gcf;
    % cbar = colorbar('Position',a.Position);
    cbar1 = colorbar;
    caxis([0 size(binColors_vertNeckSections,1)]);
    cbar1.FontName     = 'Arial';
    cbar1.FontSize     = 10;
    cbar1.Ticks        = [0.5, 1.5, 2.5, 3.5,4.5];
    % cbar1.TickLabels = num2cell(bins_neckSectionScalars(1:end-1).*rNorm);
    cbar1.TickLabels = num2cell([1 2 3 4 5]);
    % cbar1.Label.String = {sprintf('|r_0| from L_%1.0f, km',Lpoint)};
    cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar1.Label.Position = [.2, 4.2, 0];
    colormap(binColors_vertNeckSections)
    
    % -------------------
    % Plotting vertical neck sections
    % -------------------
    if plot_sectionColors == 1
        figure; hold all
        vertIndices{n_vertBins} = [];
%         yz4_indices = find(yzContourPoints4(2,:) > 0);
        plot(yzContourPoints4(1,:), yzContourPoints4(2,:),'k','linewidth',3)
        for kk = 1:n_vertBins
            %%% Indices of yzContour points in the given vertical range
            vertIndices{kk} = find(yzContourPoints4(2,:) >= bins_vertNeckSections(kk) & yzContourPoints4(2,:) <= bins_vertNeckSections(kk+1)+1e-4);
            %%% Filling vertical layer
            yFill = yzContourPoints4(1,vertIndices{kk});
            zFill_pos = yzContourPoints4(2,vertIndices{kk});
            zFill_neg = -yzContourPoints4(2,vertIndices{kk});
            
            fill(yFill, zFill_pos, binColors_vertNeckSections(kk,:))
            fill(yFill, zFill_neg, binColors_vertNeckSections(kk,:))
            plot([-y_neck_upper y_neck_upper],[0,0],'k','linewidth',1)
        end
        PlotBoi2('$y_n$','$z_n$',16,'LaTex')
        axis equal
    end
end

% -------------------------------------------------
% Lat/Lon plot colored by initial neck section (horizontal)
% -------------------------------------------------
if plot_binNeckSections_horz == 1
    % -------------------
    % Binning data
    % -------------------
    %%% Number of desired horizontal bins
    n_horzBins = 4; 
    
    %%% Choosing color for bins of neck sections 
    binColors_horzNeckSections = colorScale([colors.std.mag;colors.std.cyan],n_horzBins);
    
    %%% Find min and max bounds of bins
%     horz_min = min(impactData(:,c_y0_n));
%     horz_max = max(impactData(:,c_y0_n));
    horz_min = min(yzContourPoints4(1,:));
    horz_max = max(yzContourPoints4(1,:));
    
    %%% Assign bin values
    bins_horzNeckSections = linspace(horz_min, horz_max, n_horzBins+1);
    
    %%% Preallocate bin data structures
    binData_horzNeckSections(n_horzBins).latLons = [];

    %%% Loop through each bin and store appropriate data
    for kk = 1:n_horzBins
        %%% Indices of all impact trajectories that belong to the current bin
        currentBinIndices = find( impactData(:,c_full_y0_n) >= bins_horzNeckSections(kk) &  impactData(:,c_full_y0_n) < bins_horzNeckSections(kk+1));
        
        %%% Grabbing those trajectories and storing their impact lat/lon
        binData_horzNeckSections(kk).latLons = impactData(currentBinIndices,[c_full_impactLatitude,c_full_impactLongitude]);
    end
%     
    % -------------------
    % Plotting
    % -------------------
    figure; hold all
    for kk = 1:n_horzBins
        if isempty(binData_horzNeckSections(kk).latLons) == 0
            plot(binData_horzNeckSections(kk).latLons(:,2),binData_horzNeckSections(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_horzNeckSections(kk,:))
            if on_symmetricImpactResults == 1
                plot(binData_horzNeckSections(kk).latLons(:,2),-binData_horzNeckSections(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_horzNeckSections(kk,:))
            end
        end
    end
    PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
    xlim([-180 180])
    ylim([-90 90])

    %%% Colorbar for whole figure
    % a = gcf;
    % cbar = colorbar('Position',a.Position);
    cbar1 = colorbar;
    caxis([0 size(binColors_horzNeckSections,1)]);
    cbar1.FontName     = 'Arial';
    cbar1.FontSize     = 10;
    cbar1.Ticks        = [0.5, 1.5, 2.5, 3.5,4.5];
    % cbar1.TickLabels = num2cell(bins_neckSectionScalars(1:end-1).*rNorm);
    cbar1.TickLabels = num2cell([1 2 3 4 5]);
    % cbar1.Label.String = {sprintf('|r_0| from L_%1.0f, km',Lpoint)};
    cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar1.Label.Position = [.2, 4.2, 0];
    colormap(binColors_horzNeckSections)
    
%     % -------------------
%     % Plotting the layers in two halves
%     % -------------------
%     if rem(n_horzBins,2) == 0
%         figure; hold all
%         for kk = (n_horzBins/2):-1:1
%             if isempty(binData_horzNeckSections(kk).latLons) == 0
%                 plot(binData_horzNeckSections(kk).latLons(:,2),binData_horzNeckSections(kk).latLons(:,1),...
%                     '.','markersize',15,'color',binColors_horzNeckSections(kk,:))
%                 if on_symmetricImpactResults == 1
%                     plot(binData_horzNeckSections(kk).latLons(:,2),-binData_horzNeckSections(kk).latLons(:,1),...
%                     '.','markersize',15,'color',binColors_horzNeckSections(kk,:))
%                 end
%             end
%         end
%         PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
%         xlim([-180 180])
%         ylim([-90 90])
% 
%         %%% Colorbar for whole figure
%         % a = gcf;
%         % cbar = colorbar('Position',a.Position);
%         cbar1 = colorbar;
%         caxis([0 size(binColors_horzNeckSections,1)]);
%         cbar1.FontName     = 'Arial';
%         cbar1.FontSize     = 10;
%         cbar1.Ticks        = [0.5, 1.5, 2.5, 3.5,4.5];
%         % cbar1.TickLabels = num2cell(bins_neckSectionScalars(1:end-1).*rNorm);
%         cbar1.TickLabels = num2cell([1 2 3 4 5]);
%         % cbar1.Label.String = {sprintf('|r_0| from L_%1.0f, km',Lpoint)};
%         cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
%         cbar1.Label.Position = [.2, 4.2, 0];
%         colormap(binColors_horzNeckSections)
%         
%         figure; hold all
%         for kk = (n_horzBins/2)+1:n_horzBins
%             if isempty(binData_horzNeckSections(kk).latLons) == 0
%                 plot(binData_horzNeckSections(kk).latLons(:,2),binData_horzNeckSections(kk).latLons(:,1),...
%                     '.','markersize',15,'color',binColors_horzNeckSections(kk,:))
%                 if on_symmetricImpactResults == 1
%                     plot(binData_horzNeckSections(kk).latLons(:,2),-binData_horzNeckSections(kk).latLons(:,1),...
%                     '.','markersize',15,'color',binColors_horzNeckSections(kk,:))
%                 end
%             end
%         end
%         PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
%         xlim([-180 180])
%         ylim([-90 90])
% 
%         %%% Colorbar for whole figure
%         % a = gcf;
%         % cbar = colorbar('Position',a.Position);
%         cbar1 = colorbar;
%         caxis([0 size(binColors_horzNeckSections,1)]);
%         cbar1.FontName     = 'Arial';
%         cbar1.FontSize     = 10;
%         cbar1.Ticks        = [0.5, 1.5, 2.5, 3.5,4.5, 5.5, 6.5, 7.5];
%         % cbar1.TickLabels = num2cell(bins_neckSectionScalars(1:end-1).*rNorm);
%         cbar1.TickLabels = num2cell([1 2 3 4 5 6 7 8]);
%         % cbar1.Label.String = {sprintf('|r_0| from L_%1.0f, km',Lpoint)};
%         cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
%         cbar1.Label.Position = [.2, 4.2, 0];
%         colormap(binColors_horzNeckSections)
%     end
    
    % -------------------
    % Plotting horizontal neck sections
    % -------------------
    if plot_sectionColors == 1
        figure; hold all
        horzIndices{n_horzBins} = [];
%         yz4_indices = find(yzContourPoints4(2,:) > 0);
        plot(yzContourPoints4(1,:), yzContourPoints4(2,:),'k','linewidth',3)
        for kk = 1:n_horzBins
            %%% Indices of yzContour points in the given horizontal range
            horzIndices{kk} = find(yzContourPoints4(1,:) >= bins_horzNeckSections(kk) & yzContourPoints4(1,:) < bins_horzNeckSections(kk+1)+1e-4 & yzContourPoints4(2,:) > 0);

            %%% Filling horizontal layer
            fill([yzContourPoints4(1,horzIndices{kk}),bins_horzNeckSections(kk),bins_horzNeckSections(kk+1)],[yzContourPoints4(2,horzIndices{kk}),0,0],binColors_horzNeckSections(kk,:))
            fill([yzContourPoints4(1,horzIndices{kk}),bins_horzNeckSections(kk),bins_horzNeckSections(kk+1)],-[yzContourPoints4(2,horzIndices{kk}),0,0],binColors_horzNeckSections(kk,:))
            plot([-y_neck_upper y_neck_upper],[0,0],'k','linewidth',1)
        end
        PlotBoi2('$y_n$','$z_n$',16,'LaTex')
        axis equal
    end
end

% -------------------------------------------------
% Lat/Lon plot colored by impact angle
% -------------------------------------------------
if plot_binImpactAngles == 1
    %%% Choosing color for bins of impact angles 
    binColors_ImpactAngles = [colors.sch.s5_1(5,:);colors.sch.s5_1(4,:);...
                             colors.sch.s5_1(3,:);colors.sch.s5_1(2,:);...
                             colors.sch.s5_1(1,:)];
                     
    figure; hold all
    for kk = binCount_ImpactAngles:-1:1
%         if isempty(binData_impactAngles(kk).latLons) == 0
%             plot(binData_impactAngles(kk).latLons(:,2),binData_impactAngles(kk).latLons(:,1),...
%             '.','markersize',15,'color',binColors_ImpactAngles(kk,:))
%             if on_symmetricImpactResults == 1
%                 plot(binData_impactAngles(kk).latLons(:,2),-binData_impactAngles(kk).latLons(:,1),...
%                 '.','markersize',15,'color',binColors_ImpactAngles(kk,:))
%             end
%         end
        if isempty(binData_impactAngles(kk).latLons) == 0
            plot(binData_impactAngles(kk).latLons(:,2),binData_impactAngles(kk).latLons(:,1),...
            '.','markersize',15,'color',colors.std.black)
            if on_symmetricImpactResults == 1
                plot(binData_impactAngles(kk).latLons(:,2),-binData_impactAngles(kk).latLons(:,1),...
                '.','markersize',15,'color',colors.std.black)
            end
        end
    end
    PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
    xlim([-180 180])
    ylim([-90 90])
    %%% Colorbar for whole figure
    cbar2 = colorbar;
    caxis([0, 90]);
    cbar2.FontName     = 'Arial';
    cbar2.FontSize     = 10;
    cbar2.Ticks        = [bins_impactAngles(1:end-1), 90];
    cbar2.TickLabels = num2cell([bins_impactAngles(1:end-1), 90]);
    cbar2.Label.String = {'Impact Angle °'};
    cbar2.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar2.Label.Position = [.2, 95, 0];
    binColors_ImpactAngles_spaced = [...
        repmat(binColors_ImpactAngles(1,:),1,1);...
        repmat(binColors_ImpactAngles(2,:),3,1);...
        repmat(binColors_ImpactAngles(3,:),4,1);...
        repmat(binColors_ImpactAngles(4,:),4,1);...
        repmat(binColors_ImpactAngles(5,:),6,1)];
    colormap(binColors_ImpactAngles_spaced)
end


% -------------------------------------------------
% Metadata
% -------------------------------------------------
if plot_metaData == 1
    if on_J21 == 0
        color1 = colors.std.red;
        color2 = colors.std.black;
    elseif on_J21 == 1
        color1 = colors.std.red;
        color2 = colors.std.black;
    end
    figure(20); hold all
    plot(dvLp_mps,nLowImpactTraj,'x','linewidth',2,'markersize',10,'color',color1)
    PlotBoi2('L2 flyover speed, $mps$','Number of low-impact-angle landings',16,'LaTex')
    if on_J21 == 0
        legend('nominal')
    elseif on_J21 == 1
        legend('J21')
    end

    figure(21); hold all
    plot(dvLp_mps,nLowImpactTraj/nImpactTraj*100,'x','linewidth',2,'markersize',10,'color',color1)
    PlotBoi2('L2 flyover speed, $mps$','\% of all impacts that were low-angle',16,'LaTex')
    if on_J21 == 0
        legend('nominal')
    elseif on_J21 == 1
        legend('J21')
    end

    figure(22); hold all
    plot(dvLp_mps,nImpactTraj/nTraj*100,'x','linewidth',2,'markersize',10,'color',color1)
    PlotBoi2('L2 flyover speed, $mps$','\% of all trajectories that impacted',16,'LaTex')
    if on_J21 == 0
        legend('nominal')
    elseif on_J21 == 1
        legend('J21')
    end

    figure(23); hold all
    p1 = plot(dvLp_mps,maxLat,'x','linewidth',2,'markersize',10,'color',color1);
    p2 = plot(dvLp_mps,maxLowLat,'x','linewidth',2,'markersize',10,'color',color2);
    % legend([p1 p2],'MaxLat','MaxLowLat')
    PlotBoi2('L2 flyover speed, $mps$','Max Latitudes',16,'LaTex')
    if on_J21 == 0
        legend([p1 p2],'Max Overall Latitude','Max Low-Angle Latitude')
    elseif on_J21 == 1
        legend([p1 p2],'Max Overall Latitude (J2)','Max Low-Angle Latitude (J2)')
    end

end % plot_metaData

% ========================================================================
%%% Correlation between impact latitude and angle between v0 and y-z plane
% ========================================================================
if plot_latitudeV0AngleCorr == 1
    %%% preallocate vector for angels between v0s and YZ plane
    v0_YZ_angles_deg = zeros(nLowImpactTraj,1);
    
    %%% Vector perpendicular to YZ plane
    nz_vec = [-1, 0, 0];
    
    %%% Column indices
    col_v0  = [5, 6, 7];
    col_lat = 9;
    
    %%% Loop through low-angle impact trajectories and store v0-YZ angles
    for kk = 1:nLowImpactTraj
        v0_YZ_angles_deg(kk) = 90 - acosd(dot(nz_vec,lowImpactMat(kk,col_v0))/norm(lowImpactMat(kk,col_v0)));
    end
    
    corrcoef(v0_YZ_angles_deg,abs(lowImpactMat(:,col_lat)))
    
    %%% Plotting
    figure; hold all
    plot(v0_YZ_angles_deg,abs(lowImpactMat(:,col_lat)),'.','markersize',10,'color',colors.std.blue)
    PlotBoi2('Angle between v0 and impact latitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')

end


% ========================================================================
%%% Impact map binned by time of flight
% ========================================================================
if plot_binTOF == 1
%     %%% set column specifiers
%     c_impactBin = 2;
%     c_neckBin   = 3;
%     c_latitude  = 4;
%     c_longitude = 5;
%     c_impactAng = 6;
%     c_endTime   = 7;
%     c_y0_n      = 8;
%     c_z0_n      = 9;
%     c_v0_az     = 10;
%     c_v0_el     = 11;

%%% Finding bounds on impact time and creating bins
% minTime = min(impactData(:,c_endTime));
% maxTime = max(impactData(:,c_endTime));
% % minTime = 0;
% % maxTime = 2*pi;
% % binCount_endTime = 4;
% % bins_endTimes = linspace(minTime, maxTime, binCount_endTime + 1);
bins_endTimes = [0 pi/2 pi 3*pi/2 2*pi 4*pi];
binCount_endTime = length(bins_endTimes)-1;

binnedEndTimes = discretize(impactData(:,c_full_endTime), bins_endTimes);

binData_endTimes(binCount_endTime).latLons = [];

for kk = 1:binCount_endTime
    currentBinIndices = find(binnedEndTimes == kk);
    
    binData_endTimes(kk).latLons = impactData(currentBinIndices,[c_full_impactLatitude,c_full_impactLongitude]);
end

binColors_endTimes = colorScale([colors.std.mag;colors.std.cyan],binCount_endTime);

    figure; hold all
    for kk = binCount_endTime:-1:1
        if isempty(binData_endTimes(kk).latLons) == 0
            plot(binData_endTimes(kk).latLons(:,2),binData_endTimes(kk).latLons(:,1),...
            '.','markersize',15,'color',binColors_endTimes(kk,:))
            if on_symmetricImpactResults == 1
                plot(binData_endTimes(kk).latLons(:,2),-binData_endTimes(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_endTimes(kk,:))
            end
        end
    end
    PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
    xlim([-180 180])
    ylim([-90 90])    
    
    cbar2 = colorbar;
    caxis([0, 4*pi]);
    cbar2.FontName     = 'Arial';
    cbar2.FontSize     = 10;
    cbar2.Ticks        = [0 pi/2 pi 3*pi/2 2*pi 4*pi];
    cbar2.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi','4\pi'};
    cbar2.Label.String = {'Impact Time'};
    cbar2.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar2.Label.Position = [.2, 95, 0];
    binColors_endTimes_spaced = [...
        binColors_endTimes(1,:);...
        binColors_endTimes(2,:);...
        binColors_endTimes(3,:);...
        binColors_endTimes(4,:);...
        repmat(binColors_endTimes(5,:),4,1)];
    colormap(binColors_endTimes_spaced)
    
%         cbar2 = colorbar;
%     caxis([0, 90]);
%     cbar2.FontName     = 'Arial';
%     cbar2.FontSize     = 10;
%     cbar2.Ticks        = [bins_impactAngles(1:end-1), 90];
%     cbar2.TickLabels = num2cell([bins_impactAngles(1:end-1), 90]);
%     cbar2.Label.String = {'Impact Angle °'};
%     cbar2.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
%     cbar2.Label.Position = [.2, 95, 0];
%     binColors_ImpactAngles_spaced = [...
%         repmat(binColors_ImpactAngles(1,:),1,1);...
%         repmat(binColors_ImpactAngles(2,:),3,1);...
%         repmat(binColors_ImpactAngles(3,:),4,1);...
%         repmat(binColors_ImpactAngles(4,:),4,1);...
%         repmat(binColors_ImpactAngles(5,:),6,1)];
%     colormap(binColors_ImpactAngles_spaced)
    
    figure; hold all
    histogram(impactData(:,c_full_endTime),50,'facecolor','b','facealpha',0.7)
    PlotBoi2('Impact Time (normalized)','Frequency',16,'LaTex')
    xticks([0, pi, 2*pi, 3*pi, 4*pi])
    xticklabels({'0','\pi','2\pi','3\pi','4\pi'})
    xlim([0 4*pi])
end

if calc_maxLatitudeAnalysis == 1
%     figure
%     hist([impactData(:,c_full_impactLatitude); -impactData(:,c_full_impactLatitude)],20)
%     title(sprintf('%1.0f mps',dvLp_mps))
    
    
    
    latitudeDifferences = abs(impactData(:,c_full_impactLatitude)) - abs(impactData(:,c_full_maxLat));
    maxImpactMinusMaxLat = max(abs(impactData(:,c_full_impactLatitude))) - max(abs(impactData(:,c_full_maxLat)))
    max(latitudeDifferences);
    min(latitudeDifferences);
    
    figure; hold all
    histogram(latitudeDifferences,50,'facecolor','b','facealpha',0.7)
    PlotBoi2('$\Delta$ Latitude, $^\circ$','Frequency',16,'LaTex')
end



end


































