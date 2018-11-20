clear
clc
close all
tic
% 200, 250, 300, 350


% -------------------------------------------------
% Choosing data files
% -------------------------------------------------
%%% Outputs path
MatlabOutputsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/';

for kk = 1:7
    if kk == 1
        %%% 
        logFile            = 'F.iGS_eurL2_SR0.01_100mps_50km_149v0s_log.txt';
        impactFile         = 'F.iGS_eurL2_SR0.01_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = 'F.iGS_eurL2_SR0.01_100mps_50km_149v0s_land.txt';
    end
    if kk == 2
        %%% 
        logFile            = 'F.iGS_eurL2_SR0.04_100mps_50km_149v0s_log.txt';
        impactFile         = 'F.iGS_eurL2_SR0.04_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = 'F.iGS_eurL2_SR0.04_100mps_50km_149v0s_land.txt';
    end
    if kk == 3
        %%% 
        logFile            = 'F.iGS_eurL2_SR0.07_100mps_50km_149v0s_log.txt';
        impactFile         = 'F.iGS_eurL2_SR0.07_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = 'F.iGS_eurL2_SR0.07_100mps_50km_149v0s_land.txt';
    end
    if kk == 4
        %%% 
        logFile            = 'F.iGS_eurL2_SR0.10_100mps_50km_149v0s_log.txt';
        impactFile         = 'F.iGS_eurL2_SR0.10_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = 'F.iGS_eurL2_SR0.10_100mps_50km_149v0s_land.txt';
    end
    if kk == 5
        %%% 
        logFile            = 'F.iGS_eurL2_SR0.13_100mps_50km_149v0s_log.txt';
        impactFile         = 'F.iGS_eurL2_SR0.13_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = 'F.iGS_eurL2_SR0.13_100mps_50km_149v0s_land.txt';
    end
    if kk == 6
        %%% 
        logFile            = 'F.iGS_eurL2_SR0.16_100mps_50km_149v0s_log.txt';
        impactFile         = 'F.iGS_eurL2_SR0.16_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = 'F.iGS_eurL2_SR0.16_100mps_50km_149v0s_land.txt';
    end
    if kk == 7
        %%% 
        logFile            = 'F.iGS_eurL2_SR0.19_100mps_50km_149v0s_log.txt';
        impactFile         = 'F.iGS_eurL2_SR0.19_100mps_50km_149v0s_data.txt';
        lowImpactAngleFile = 'F.iGS_eurL2_SR0.19_100mps_50km_149v0s_land.txt';
    end
    
%     if kk == 8
%         %%% 
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
%     end
%     
%     if kk == 9
%         %%% 
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
%     end
%     
%     if kk == 10
%         %%% 
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
%     end
%     
%     if kk == 11
%         %%% 
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
%     end
%     
%     if kk == 12
%         %%% 
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
%     end
%     
%     if kk == 13
%         %%% 
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
%     end
%     
%     if kk == 14
%         %%% 
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
%     end
    
    %%% Test case
    if kk == 989
%         logFile            = '';
%         impactFile         = '';
%         lowImpactAngleFile = '';
    end

    




% -------------------------------------------------
% Running
% -------------------------------------------------
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
bodies = getBodyData(mbinPath);

iGS_LoadPlot2([MatlabOutputsPath,logFile],[MatlabOutputsPath,impactFile],[MatlabOutputsPath,lowImpactAngleFile],mbinPath)


end
toc





% ========================================================================
%%% Load/Plot Function
% ========================================================================
function iGS_LoadPlot2(logFile,impactFile,lowImpactAngleFile,mbinPath)
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
plot_fullLowTrajectories = 0;
plot_LowTrajectories_r0  = 0;
plot_initialConditions   = 0;
plot_fullSystemContour   = 0;
plot_sectionColors       = 0;

plot_binImpactAngles     = 0;
plot_binNeckSections     = 0;
plot_metaData            = 0;
plot_latitudeV0AngleCorr = 0;

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

if isequal(logData{1}{kk},'SizeRatio') == 1
    SR = str2num(logData{2}{kk});
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
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Correcting R2_n to current size ratio
R2_n = SR * (L123(2,1) - L123(1,1));
% ========================================================================
%%% Loading General Impact Data
% ========================================================================
%%% Load impact data
impactMat = dlmread(impactFile,',',1,0);
nImpactTraj = size(impactMat,1);
if impactMat(end,:) == 0
    nImpactTraj = 0;
end
%%% set column specifiers
c_impactBin = 1;
c_neckBin   = 2;
c_latitude  = 3;
c_longitude = 4;
c_impactAng = 5;

%%% Bin counts
binCount_ImpactAngles = length(bins_impactAngles) - 1;
binCount_NeckSections = length(bins_neckSectionScalars);

%%% Preallocate
binData_impactAngles(binCount_ImpactAngles).latLons = [];
binData_neckSections(length(bins_neckSectionScalars)).latLons = [];

%%% Store impact angle bins
for kk = 1:binCount_ImpactAngles
    %%% Indices of all impact trajectories that belong to the current bin
    currentBinIndices = find(impactMat(:,c_impactBin)==kk);
    
    %%% Grabbing those trajectories and storing their impact lat/lon
    binData_impactAngles(kk).latLons = impactMat(currentBinIndices,[c_latitude,c_longitude]);
end

%%% Store neck section bins
for kk = 1:binCount_NeckSections
    %%% Indices of all impact trajectories that belong to the current bin
    currentBinIndices = find(impactMat(:,c_neckBin)==kk);
    
    %%% Grabbing those trajectories and storing their impact lat/lon
    binData_neckSections(kk).latLons = impactMat(currentBinIndices,[c_latitude,c_longitude]);
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
if lowImpactMat(end,:) == 0
    nLowImpactTraj = 0;
end

nLowImpactTraj/nImpactTraj

if size(lowImpactMat) ~= [1,1] % if there are any low-angle impacts:
    
    % -------------------------------------------------
    % Propagate low-impact-angle trajectories
    % -------------------------------------------------
    if plot_fullLowTrajectories == 1 || plot_LowTrajectories_r0 == 1
    %%% Preallocating
    lowTrajs{nLowImpactTraj} = [];

    %%% Setting time vector
    t_i = 0; % sec
    t_f = 4*pi; % Long bc events are watching for impact or escape
    dt = t_f/10000;
    time0_n = t_i:dt:t_f;

    %%% Choosing ode tolerance
    tol = 2.22045e-14;

    %%% Setting integrator options
    options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

    %%% Setting necessary parameters for integration
    prms.u = secondary.MR;
    prms.R2_n = R2_n;
    prms.L1x = L123(1,1);
    prms.L2x = L123(2,1);

    for kk = 1:size(lowImpactMat,1)
        %%% Propagating trajectory
        [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, lowImpactMat(kk,1:6)', options_ImpactEscape, prms);

        lowTrajs{kk} = [X_BCR_n, time_n];
    end
    
    %%% Plotting full low-landing-angle trajectories
    if plot_fullLowTrajectories == 1
        figure; hold all
        plotBodyTexture3(R2_n,[1-secondary.MR,0,0],secondary.img)
        PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
        axis equal
        for kk = 1:length(lowTrajs)
            plot3(lowTrajs{kk}(:,1),lowTrajs{kk}(:,2),lowTrajs{kk}(:,3),'b')
        end
    end
    
    %%% Plotting initial positions of low-landing-angle trajectories
    if plot_LowTrajectories_r0 == 1
        figure; hold all
        plotBodyTexture3(R2_n,[1-secondary.MR,0,0],secondary.img)
        PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
        axis equal
        for kk = 1:length(lowTrajs)
            plot3(lowTrajs{kk}(1,1),lowTrajs{kk}(1,2),lowTrajs{kk}(1,3),'b.','markersize',16)
        end
    end
    end % plot_fullLowTrajectories == 1 || plot_LowTrajectories_r0 == 1
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
        zv = JacobiConstantCalculator(secondary.MR,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        
        JCs_yz_Lpoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzContourPoints4 ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC_scInitial );

yzContourPoints3 = yzContourPoints4.*bins_neckSectionScalars(3);
yzContourPoints2 = yzContourPoints4.*bins_neckSectionScalars(2);
yzContourPoints1 = yzContourPoints4.*bins_neckSectionScalars(1);



%%% Choosing color for bins of impact angles 
binColors_ImpactAngles = [colors.sch.s5_1(5,:);colors.sch.s5_1(4,:);...
                         colors.sch.s5_1(3,:);colors.sch.s5_1(2,:);...
                         colors.sch.s5_1(1,:)];
                     
%%% Choosing color for bins of neck sections 
binColors_NeckSections = [colors.sch.d4_1(2,:);colors.sch.d4_1(1,:);...
                         colors.sch.d4_1(4,:);colors.sch.d4_1(3,:)];

% ========================================================================
%%% Plotting Data
% ========================================================================

% -------------------------------------------------
% Initial conditions in y-z plane along with bounding contours
% -------------------------------------------------
if plot_initialConditions == 1
    %%% Creating ZV contours
    % Contour bounds
    xCont_min = L123(1,1)-2*secondary.R_n;
    xCont_max = L123(2,1)+2*secondary.R_n;
    yCont_min = -9*secondary.R_n;
    yCont_max = 9*secondary.R_n;

    % Creating x-y grid
    xs = linspace(xCont_min,xCont_max,600);
    ys = linspace(yCont_min,yCont_max,200);
    [X_xy, Y_xy] = meshgrid(xs,ys);
    clear xs ys zs

    % Calculating JCs across x-y grid
    JCs_xy = zeros(size(X_xy));
    for xk = 1:size(X_xy,1)
        for yk = 1:size(X_xy,2)
            %%% Zero-Velocity Curve
            zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0], [0, 0, 0]);

            JCs_xy(xk,yk) = zv;
        end
    end

    figure; hold all
    plot3(ones(size(yzContourPoints4,2),1).*L123(Lpoint,1), yzContourPoints4(1,:),yzContourPoints4(2,:),'k','linewidth',3)
    plotBodyTexture3(R2_n, [1-secondary.MR, 0, 0], secondary.img)
    [xyContourPoints,href] = contour(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],...
        'color',colors.std.black,'linewidth',3);
    PlotBoi3('x$_n$','y$_n$','z$_n$',14,'LaTex')
    view(-35,24)
    camva(9)
    axis equal
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
            zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0], [0, 0, 0]);
            
            JCs_xy(xk,yk) = zv;
        end
    end

    figure; hold all
    plotBodyTexture3(R2_n, [1-secondary.MR, 0, 0], secondary.img)
    plotBodyTexture3(primary.R./rNorm, [-secondary.MR, 0, 0], primary.img)
    [xyContourPoints,href] = contourf(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],...
        'color',colors.std.black,'linewidth',1.5);
    colormap(colors.sch.r6(6,:))
    PlotBoi3('X','Y','Z',14)
    view(0,90)
    camva(9)
    axis equal
end

% -------------------------------------------------
% Plot showing the initial neck section color coding
% -------------------------------------------------
if plot_sectionColors == 1
    figure; hold all
    plot(yzContourPoints4(1,:), yzContourPoints4(2,:),'k','linewidth',3)
    fill(yzContourPoints4(1,:),yzContourPoints4(2,:),binColors_NeckSections(4,:))
    fill(yzContourPoints3(1,:),yzContourPoints3(2,:),binColors_NeckSections(3,:))
    fill(yzContourPoints2(1,:),yzContourPoints2(2,:),binColors_NeckSections(2,:))
    fill(yzContourPoints1(1,:),yzContourPoints1(2,:),binColors_NeckSections(1,:))
    

%     for kk = binCount_NeckSections-1:-1:1
%         plotBody2(bins_neckSections(kk+1),[0,0,0], binColors_NeckSections(kk,:), binColors_NeckSections(kk,:), 1)
%     end

    axis equal
    PlotBoi3('Y','Z','X',14)
end

% -------------------------------------------------
% Lat/Lon plot colored by initial neck section
% -------------------------------------------------
if plot_binNeckSections == 1
    figure; hold all
    for kk = 1:4
        if isempty(binData_neckSections(kk).latLons) == 0
            plot(binData_neckSections(kk).latLons(:,2),binData_neckSections(kk).latLons(:,1),...
                '.','markersize',15,'color',binColors_NeckSections(kk,:))
        end
    end
    title(sprintf('Size Ratio: %1.6f',SR));
    PlotBoi2('Longitude, deg','Latitude, deg',14)
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
end

% -------------------------------------------------
% Lat/Lon plot colored by impact angle
% -------------------------------------------------
if plot_binImpactAngles == 1
    figure; hold all
    for kk = binCount_ImpactAngles:-1:1
    if isempty(binData_impactAngles(kk).latLons) == 0
    plot(binData_impactAngles(kk).latLons(:,2),binData_impactAngles(kk).latLons(:,1),...
    '.','markersize',15,'color',binColors_ImpactAngles(kk,:))
    end
    end
    title(sprintf('Size Ratio: %1.6f',SR));
    PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',14,'LaTex')
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
    color1 = colors.std.red;
    color2 = colors.std.black;

    figure(20); hold all
    plot(dvLp_mps,nLowImpactTraj,'x','linewidth',2,'markersize',10,'color',color1)
    PlotBoi2('L2 flyover speed, $mps$','Number of low-impact-angle landings',14,'LaTex')

    figure(21); hold all
    plot(dvLp_mps,nLowImpactTraj/nImpactTraj*100,'x','linewidth',2,'markersize',10,'color',color1)
    PlotBoi2('L2 flyover speed, $mps$','\% of all impacts that were low-angle',14,'LaTex')

    figure(22); hold all
    plot(dvLp_mps,nImpactTraj/nTraj*100,'x','linewidth',2,'markersize',10,'color',color1)
    PlotBoi2('L2 flyover speed, $mps$','\% of all trajectories that impacted',14,'LaTex')

    figure(23); hold all
    p1 = plot(dvLp_mps,maxLat,'x','linewidth',2,'markersize',10,'color',color1);
    p2 = plot(dvLp_mps,maxLowLat,'x','linewidth',2,'markersize',10,'color',color2);
    % legend([p1 p2],'MaxLat','MaxLowLat')
    PlotBoi2('L2 flyover speed, $mps$','Max Latitudes',14,'LaTex')

end % plot_metaData

% ========================================================================
%%% Correlation between impact latitude and angle between v0 and y-z plane
% ========================================================================
if plot_latitudeV0AngleCorr == 1
    %%% preallocate vector for angels between v0s and YZ plane
    v0_YZ_angles_deg = zeros(size(lowImpactMat,1),1);
    
    %%% Vector perpendicular to YZ plane
    nz_vec = [-1, 0, 0];
    
    %%% Column indices
    col_v0  = [4, 5, 6];
    col_lat = 8;
    
    %%% Loop through low-angle impact trajectories and store v0-YZ angles
    for kk = 1:size(lowImpactMat,1)
        v0_YZ_angles_deg(kk) = 90 - acosd(dot(nz_vec,lowImpactMat(kk,col_v0))/norm(lowImpactMat(kk,col_v0)));
    end
    
    corrcoef(v0_YZ_angles_deg,abs(lowImpactMat(:,col_lat)))
    
    %%% Plotting
    figure; hold all
    plot(v0_YZ_angles_deg,abs(lowImpactMat(:,col_lat)),'.','markersize',10,'color',colors.std.blue)
    PlotBoi2('Angle between v0 and impact latitude, $^\circ$','Latitude, $^\circ$',14,'LaTex')

end


end


































