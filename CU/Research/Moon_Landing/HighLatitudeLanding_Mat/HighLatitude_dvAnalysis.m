% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
ticWhole = tic;

dataPath_test = '~/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/TestDatabases';
dataPath = '/Volumes/LB_External_Drive/Research/High_Latitude_Landing_Study/Data/Full_Results/';
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
saveData = false;
plotData = false;
makeAndPlotDVBins = false;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Load in datasets
% -------------------------------------------------
%%% High latitude trajectories
% data_high = dlmread([dataPath_test, '/', 'highLatitudeTrajectories.txt'],',',1,0);
load([dataPath, '/', 'highTrajs_50mps_70lat_1pi.mat']) % data_high_1pi;

n_highPoints = size(data_high_1pi,1);

data_high_1pi = data_high_1pi(1:2:n_highPoints, :);
% 989
% data_high_1pi = data_high_1pi(1:100000, :);
n_highPoints = size(data_high_1pi,1);

%%% Low latitude trajectories
% data_low = dlmread([dataPath_test, '/', 'lowLatitudeTrajectory.txt'],',',1,0);

load([dataPath, '/', 'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi.mat']); % data_L2_4pi
n_lowPoints = size(data_L2_4pi,1);

data_L2_4pi = data_L2_4pi(1:4:n_lowPoints, :);
n_lowPoints = size(data_L2_4pi,1);

%%% Column identifiers
c_trajID = 1;   
c_x      = 2;
c_y      = 3;
c_z      = 4;
c_xd     = 5;
c_yd     = 6;
c_zd     = 7;
c_t      = 8;
% -------------------------------------------------
%%% System setup
% -------------------------------------------------
%%% System
primary = bodies.jupiter;   secondary = bodies.europa;

%%% Normalizing factors
rNorm = secondary.a;         % n <-> km

tNorm_1 = sqrt((secondary.a^3)/(primary.u + secondary.u));
tNorm_2 = 1 / secondary.meanMot;

vNorm_1 = rNorm / tNorm_1;
vNorm_2 = rNorm / tNorm_2;

% tNorm = 1/secondary.meanMot; % n <-> sec
% vNorm = rNorm / tNorm;       % n <-> km/sec

%%% prms structure
prms.u = secondary.MR;
prms.n = 1;

%%% Collinear equilibrium points
rLP = EquilibriumPoints(secondary.MR, prms.n, 1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% x-position of secondary body
x_body = 1-prms.u;

% -------------------------------------------------
%%% Options
% -------------------------------------------------
%%% "Low latitude" states must be at least this close to the secondary body
%%% to be considered
maximumRadiusFromSecondary_n = 8*secondary.R_n;

%%% Maximum allowable distance between "low" and "high" states to be
%%% considered co-located
maximumRadiusBetweenPoints_n = 100/rNorm; 

% ========================================================================
%%% Create octants in High-Latitude dataset and get the corresponding
%%% indices
% ========================================================================
%%% Set bound for highest and lowest points of the low trajectory
zMax_low = max(data_L2_4pi(:,c_z)) + 1/rNorm;
zMin_low = min(data_L2_4pi(:,c_z)) - 1/rNorm;

highPointsOctants = cell(n_highPoints,1);

parfor point_k = 1:n_highPoints
    highPointsOctants{point_k}.sec1 = NaN(6,1);
    highPointsOctants{point_k}.sec2 = NaN(6,1);
    highPointsOctants{point_k}.sec3 = NaN(6,1);
    highPointsOctants{point_k}.sec4 = NaN(6,1);
    highPointsOctants{point_k}.sec5 = NaN(6,1);
    highPointsOctants{point_k}.sec6 = NaN(6,1);
    highPointsOctants{point_k}.sec7 = NaN(6,1);
    highPointsOctants{point_k}.sec8 = NaN(6,1);
    
    if isnan(data_high_1pi(point_k,c_x))
        continue
    end
    
    %%% If this point is outside of the z range, skip it
    if (data_high_1pi(point_k,c_z) > zMax_low) || (data_high_1pi(point_k,c_z) < zMin_low)
        continue
    end
    
    %%% Assign octant number of current point
    currentoctantNumber = getoctantNumber(data_high_1pi(point_k,c_x:c_z), x_body);
    
    %%% Store current point in the proper field based on octant number
    if currentoctantNumber == 1
        highPointsOctants{point_k}.sec1 = data_high_1pi(point_k,c_x:c_zd)';
        
    elseif currentoctantNumber == 2
        highPointsOctants{point_k}.sec2 = data_high_1pi(point_k,c_x:c_zd)';
        
    elseif currentoctantNumber == 3
        highPointsOctants{point_k}.sec3 = data_high_1pi(point_k,c_x:c_zd)';
        
    elseif currentoctantNumber == 4
        highPointsOctants{point_k}.sec4 = data_high_1pi(point_k,c_x:c_zd)';
        
    elseif currentoctantNumber == 5
        highPointsOctants{point_k}.sec5 = data_high_1pi(point_k,c_x:c_zd)';
        
    elseif currentoctantNumber == 6
        highPointsOctants{point_k}.sec6 = data_high_1pi(point_k,c_x:c_zd)';
        
    elseif currentoctantNumber == 7
        highPointsOctants{point_k}.sec7 = data_high_1pi(point_k,c_x:c_zd)';
        
    elseif currentoctantNumber == 8
        highPointsOctants{point_k}.sec8 = data_high_1pi(point_k,c_x:c_zd)';
    end
    
end

clear data_high_1pi

%%% Making octant data more accessible
highPointsOctants_structArray = [highPointsOctants{:}];

clear highPointsOctants

highPoints_sec1 = [highPointsOctants_structArray(:).sec1];
highPoints_sec1 = highPoints_sec1(:,~isnan(highPoints_sec1(1,:)))';

highPoints_sec2 = [highPointsOctants_structArray(:).sec2];
highPoints_sec2 = highPoints_sec2(:,~isnan(highPoints_sec2(1,:)))';

highPoints_sec3 = [highPointsOctants_structArray(:).sec3];
highPoints_sec3 = highPoints_sec3(:,~isnan(highPoints_sec3(1,:)))';

highPoints_sec4 = [highPointsOctants_structArray(:).sec4];
highPoints_sec4 = highPoints_sec4(:,~isnan(highPoints_sec4(1,:)))';

highPoints_sec5 = [highPointsOctants_structArray(:).sec5];
highPoints_sec5 = highPoints_sec5(:,~isnan(highPoints_sec5(1,:)))';

highPoints_sec6 = [highPointsOctants_structArray(:).sec6];
highPoints_sec6 = highPoints_sec6(:,~isnan(highPoints_sec6(1,:)))';

highPoints_sec7 = [highPointsOctants_structArray(:).sec7];
highPoints_sec7 = highPoints_sec7(:,~isnan(highPoints_sec7(1,:)))';

highPoints_sec8 = [highPointsOctants_structArray(:).sec8];
highPoints_sec8 = highPoints_sec8(:,~isnan(highPoints_sec8(1,:)))';

clear highPointsOctants_structArray

fprintf('Octants created: %1.4f seconds\n',toc(ticWhole))
% ========================================================================
%%% For every point on the "low" trajectory within some radius of the
%%% secondary body, find the nearest "high" state, if that's within some
%%% distance, find the DV (assuming they were at the same location), and
%%% store all the DV information
% ========================================================================
positionsDistancesOctantsOffsetsAndDVs_cellArray = cell(1, n_lowPoints);

%%
for lowPoint_k = 1:n_lowPoints
    %%% Skip if NaN
    if isnan(data_L2_4pi(lowPoint_k,c_x))
        continue
    end
    
    %%% Check if this point is within the maximum allowable distance of the
    %%% secondary body
    r2 = sqrt((data_L2_4pi(lowPoint_k,c_x)-(1-prms.u))^2 + data_L2_4pi(lowPoint_k,c_y)^2 + data_L2_4pi(lowPoint_k,c_z)^2);
    if (r2 > maximumRadiusFromSecondary_n) || (r2 < secondary.R_n)
        continue
    end
    
    %%% Assign octant number of current point
    currentoctantNumber = getoctantNumber(data_L2_4pi(lowPoint_k,c_x:c_z), x_body);
    
    %%% Determine current octant of high points
    highPointsInOctant = [];
    if currentoctantNumber == 1
        highPointsInOctant = highPoints_sec1;
        
    elseif currentoctantNumber == 2
        highPointsInOctant = highPoints_sec2;
        
    elseif currentoctantNumber == 3
        highPointsInOctant = highPoints_sec3;
        
    elseif currentoctantNumber == 4
        highPointsInOctant = highPoints_sec4;
        
    elseif currentoctantNumber == 5
        highPointsInOctant = highPoints_sec5;
        
    elseif currentoctantNumber == 6
        highPointsInOctant = highPoints_sec6;
        
    elseif currentoctantNumber == 7
        highPointsInOctant = highPoints_sec7;
        
    elseif currentoctantNumber == 8
        highPointsInOctant = highPoints_sec8;
    end
    
    %%% Loop through red points and find the nearest one
    distances_n = ones(size(highPointsInOctant,1),1);
    for highPoint_k = 1:size(highPointsInOctant,1)
        distances_n(highPoint_k) = sqrt((data_L2_4pi(lowPoint_k,c_x) - highPointsInOctant(highPoint_k,1))^2 + (data_L2_4pi(lowPoint_k,c_y) - highPointsInOctant(highPoint_k,2))^2 + (data_L2_4pi(lowPoint_k,c_z) - highPointsInOctant(highPoint_k,3))^2);
    end
    
    nearestPointDistance_n = min(distances_n);
    
    %%% Check if this point is within the maximum allowable point-to-point
    %%% distance
    if nearestPointDistance_n > maximumRadiusBetweenPoints_n
        continue
    end
    
    %%% Store results
    nearestPointIndex = find(distances_n == nearestPointDistance_n);
    DV_n = sqrt((data_L2_4pi(lowPoint_k,c_xd) - highPointsInOctant(highPoint_k,4))^2 + (data_L2_4pi(lowPoint_k,c_yd) - highPointsInOctant(highPoint_k,5))^2 + (data_L2_4pi(lowPoint_k,c_zd) - highPointsInOctant(highPoint_k,6))^2);
    positionsDistancesOctantsOffsetsAndDVs_cellArray{lowPoint_k} = [data_L2_4pi(lowPoint_k,c_x:c_z), r2, currentoctantNumber, nearestPointDistance_n, DV_n]';
    
end

clear highPoints_sec1 highPoints_sec2 highPoints_sec3 highPoints_sec4 highPoints_sec5 highPoints_sec6 highPoints_sec7 highPoints_sec8


fprintf('Results Generated: %1.4f seconds\n',toc(ticWhole))
positionsDistancesOctantsOffsetsAndDVs = [positionsDistancesOctantsOffsetsAndDVs_cellArray{:}]';


% --------------------------
% Column specifiers for data matrix
% --------------------------
c_x                  = 1;
c_y                  = 2;
c_z                  = 3;
c_r2                 = 4;
c_octant             = 5;
c_nearestPointOffset = 6;
c_DV                 = 7;


if saveData
    %%% Create savepath
    savePath = dataPath;
    
    % -------------------------------------------------
    %%% Preparing Save File for new family (southern or northern)
    % -------------------------------------------------
    fileAlreadyExists   = 1;
    fileVersion         = 0;
    newFileName_new = ['positionsDistancesOctantsOffsetsAndDVs'];
    while fileAlreadyExists == 1
        if isfile([savePath,newFileName_new,'.txt']) == 1
            fileVersion = fileVersion + 1;
            newFileName_new = sprintf('%s_%1d',newFileName_new,fileVersion);
        else
            fileAlreadyExists = 0;
        end
    end
    newFileName = [savePath,newFileName_new,'.txt'];
    
    %%% Open File
    newfile = fopen(newFileName,'wt');
   
    % -------------------------------------------------
    %%% Writing data
    % -------------------------------------------------
    %%% Write header
    fprintf(newfile,['x_n,y_n,z_n,r2_n,octant,offsetBetweenPoints_n,DV_n\n']);
    
    %%% Write data
    dataStart = 1;
    dataStop = size(positionsDistancesOctantsOffsetsAndDVs,1);
    
    
    for kk = dataStart:dataStop
        fprintf(newfile,'%1.10f, %1.10f, %1.10f, %1.9f, %1d, %1.9f, %1.5f\n',...
            positionsDistancesOctantsOffsetsAndDVs(kk,1), positionsDistancesOctantsOffsetsAndDVs(kk,2),...
            positionsDistancesOctantsOffsetsAndDVs(kk,3), positionsDistancesOctantsOffsetsAndDVs(kk,4),...
            positionsDistancesOctantsOffsetsAndDVs(kk,5), positionsDistancesOctantsOffsetsAndDVs(kk,6),...
            positionsDistancesOctantsOffsetsAndDVs(kk,7));
            
    end

    %%% Close file
    fclose(newfile);

    save([savePath,'positionsDistancesOctantsOffsetsAndDVs','.mat'], 'positionsDistancesOctantsOffsetsAndDVs','-v7.3')
end


% ========================================================================
%%% Plot results
% ========================================================================
if plotData
    if ~exist('positionsDistancesOctantsOffsetsAndDVs')
        load([dataPath,'positionsDistancesOctantsOffsetsAndDVs_50mps_4piL2.mat']); % data_L1_4pi
    end
    
    % --------------------------
    % Positions of all solution points
    % --------------------------
    figure; hold all
    plot3(positionsDistancesOctantsOffsetsAndDVs(:,c_x),positionsDistancesOctantsOffsetsAndDVs(:,c_y),positionsDistancesOctantsOffsetsAndDVs(:,c_z),'b.')
    PlotBoi3_CR3Bn(20)
    axis equal
    plotSecondary(secondary)
    
    
    % --------------------------
    % X-Y projection of solution points, DV on the z-axis
    % --------------------------
    figure; hold all
    plot3(positionsDistancesOctantsOffsetsAndDVs(:,c_x),positionsDistancesOctantsOffsetsAndDVs(:,c_y),positionsDistancesOctantsOffsetsAndDVs(:,c_DV),'b.')
    PlotBoi3('$x_n$','$y_n$','$\Delta V$',20,'LaTex')
    
    % --------------------------
    % Box plot of DV values per section (octant)
    % --------------------------
    figure; hold all
    boxplot(positionsDistancesOctantsOffsetsAndDVs(:,c_DV).*vNorm_1, positionsDistancesOctantsOffsetsAndDVs(:,c_octant), 'outliersize',2, 'jitter',.4,'whisker',5) % whisker default 1.5
    PlotBoi2('Octant','$\Delta V$, $km/s$',20, 'LaTex')
    % https://www.mathworks.com/help/stats/boxplot.html#bu3pip4
    % "The default value for 'Whisker' corresponds to approximately +/?2.7?
    % and 99.3 percent coverage if the data are normally distributed. The plotted whisker extends to the adjacent value, which is the most extreme data value that is not an outlier."
    
    
    if makeAndPlotDVBins

        %%% DV bins
        n_DV_bins_plus1 = 5;
        DV_bin_values = linspace(min(positionsDistancesOctantsOffsetsAndDVs(:,c_DV)-1e-5), max(positionsDistancesOctantsOffsetsAndDVs(:,c_DV)+1e-5), n_DV_bins_plus1);
        DV_bin_values_kps = DV_bin_values .* vNorm_1;
        %%% Bin colors
%         binColors = colorScale([colors.blue2; colors.red3],n_DV_bins_plus1-1);
%         binColors = [colors.cyan;...
%                      colors.blue2;...
%                      colors.mag;...
%                      colors.red3];
        binColors = colors.sch.d4_1;
        
        
        DVbins = zeros(size(positionsDistancesOctantsOffsetsAndDVs,1),1);
        for kk = 1:size(positionsDistancesOctantsOffsetsAndDVs,1)
            DVbins(kk) = discretize(positionsDistancesOctantsOffsetsAndDVs(kk,c_DV),DV_bin_values);
        end
        
        DVBinData = cell(n_DV_bins_plus1-1,1);
        for kk = 1:(n_DV_bins_plus1-1)
            DVBinData{kk} = positionsDistancesOctantsOffsetsAndDVs(DVbins==kk,:);
        end
        
        
        
        
%         figure('position',[440 315 695 483]); hold all
        figure; hold all
        for kk = 1:(n_DV_bins_plus1-1)
%         for kk = (n_DV_bins_plus1-1):-1:1
%             plot3(DVBinData{kk}(:,c_x),DVBinData{kk}(:,c_y),DVBinData{kk}(:,c_DV),'.','markersize',8,'color',binColors(kk,:))
            plot(DVBinData{kk}(:,c_x),DVBinData{kk}(:,c_y),'.','markersize',8,'color',binColors(kk,:))
        end
        axis equal
        PlotBoi3_CR3Bn(20)
        
        plotSecondary(secondary)
        view(0,90)
        
%         cbar2 = colorbar;
%         caxis([min(positionsDistancesOctantsOffsetsAndDVs(:,c_DV))-1e-5, max(positionsDistancesOctantsOffsetsAndDVs(:,c_DV))+1e-5].*vNorm_1);
%         cbar2.FontName       = 'Arial';
%         cbar2.FontSize       = 10;
%         cbar2.Ticks          = DV_bin_values_kps;
%         cbar2.TickLabels     = num2cell(DV_bin_values_kps);
%         cbar2.Label.String   = {'\DeltaV, km/s'};
%         cbar2.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
%         cbar2.Label.Position = [0.6, 3.67, 0];
%         colormap(binColors)
        
        
        figure('position',[713 268 604 456]); hold all
%         for kk = 1:(n_DV_bins_plus1-1)
        for kk = (n_DV_bins_plus1-1):-1:1
%             plot3(DVBinData{kk}(:,c_x),DVBinData{kk}(:,c_y),DVBinData{kk}(:,c_DV),'.','markersize',8,'color',binColors(kk,:))
            plot(DVBinData{kk}(:,c_x),DVBinData{kk}(:,c_y),'.','markersize',8,'color',binColors(kk,:))
        end
        axis equal
        PlotBoi3_CR3Bn(20)
        
        
        cbar2 = colorbar;
        caxis([min(positionsDistancesOctantsOffsetsAndDVs(:,c_DV))-1e-5, max(positionsDistancesOctantsOffsetsAndDVs(:,c_DV))+1e-5].*vNorm_1);
        cbar2.FontName       = 'Arial';
        cbar2.FontSize       = 12; % 10
        cbar2.Ticks          = DV_bin_values_kps;
        cbar2.TickLabels     = num2cell(DV_bin_values_kps);
        cbar2.Label.String   = {'\DeltaV, km/s'};
        cbar2.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
        cbar2.Label.Position = [0.6, 3.67, 0];
        colormap(binColors)
        
        plotSecondary(secondary)
        view(0,90)
        
        
        
        
        
        
        
        

    
    end
    
end

% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)










function [octant_number] = getoctantNumber(r, x_body)
    if r(3) >= 0
        if r(1) >= x_body
            if r(2) >= 0
                octant_number = 1;
            elseif r(2) < 0
                octant_number = 4;
            end
        elseif r(1) < x_body
            if r(2) >= 0
                octant_number = 2;
            elseif r(2) < 0
                octant_number = 3;
            end
        end
    elseif r(3) < 0
        if r(1) >= x_body
            if r(2) >= 0
                octant_number = 5;
            elseif r(2) < 0
                octant_number = 8;
            end
        elseif r(1) < x_body
            if r(2) >= 0
                octant_number = 6;
            elseif r(2) < 0
                octant_number = 7;
            end
        end
    end
end







