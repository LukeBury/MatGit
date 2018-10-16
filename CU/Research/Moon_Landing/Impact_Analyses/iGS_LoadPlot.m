clear
clc
close all

logFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/testFile_log.txt';
impactFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/testFile_data.txt';
lowImpactAngleFile = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/testFile_land.txt';

mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
bodies = getBodyData(mbinPath);
primary = bodies.jupiter;
secondary = bodies.europa;

bins_impactAngles = [0, 5, 20, 40, 60, 91];
bins_neckSections = [0, 0.002555110294761, 0.005110220589522, 0.007665330884284, 0.009204876638791];
JC_scInitial = 3.003397404246083;
Lpoint = 2;
y_neck_upper = 0.009204876638791;
z_neck_upper = 0.007665330884284;



iGS_LoadPlot2(logFile,impactFile,lowImpactAngleFile,primary,secondary,bins_impactAngles,bins_neckSections,JC_scInitial,Lpoint,y_neck_upper,z_neck_upper)

function iGS_LoadPlot2(logFile,impactFile,lowImpactAngleFile,primary,secondary,bins_impactAngles,bins_neckSections,JC_scInitial,Lpoint,y_neck_upper,z_neck_upper)
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
% File names
% -------------------------------------------------
% logFile    = 'F.iGS_EurL2_200mps_50km_156v0s_log.txt';
% impactFile = 'F.iGS_EurL2_200mps_50km_156v0s_data.txt';
% lowImpactAngleFile = 'F.iGS_EurL2_200mps_50km_156v0s_land.txt';


% -------------------------------------------------
% Run Switches
% -------------------------------------------------
plot_initialConditions = 0;
plot_fullSystemContour = 0;
plot_sectionColors     = 0;

% ========================================================================
%%% Setting data
% ========================================================================
% %%% Load log data
% logData = dlmread(logFile,':',[5 2 inf 2]);

%%% Color options/schemes
colors = get_colors();

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Find equilibrium points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% ========================================================================
%%% Loading General Impact Data
% ========================================================================
%%% Load impact data
impactData = dlmread(impactFile,',',1,0);

%%% set column specifiers
c_impactBin = 1;
c_neckBin   = 2;
c_latitude  = 3;
c_longitude = 4;
c_impactAng = 5;

%%% Bin counts
binCount_NeckSections = length(bins_neckSections) - 1;
binCount_ImpactAngles = length(bins_impactAngles) - 1;

%%% Preallocate
binData_impactAngles(binCount_ImpactAngles).latLons = [];
binData_neckSections(length(bins_neckSections)).latLons = [];

%%% Store impact angle bins
for kk = 1:binCount_ImpactAngles
    currentBinIndices = find(impactData(:,c_impactBin)==kk);
    binData_impactAngles(kk).latLons = impactData(currentBinIndices,[c_latitude,c_longitude]);
end

for kk = 1:binCount_NeckSections
    currentBinIndices = find(impactData(:,c_neckBin)==kk);
    binData_neckSections(kk).latLons = impactData(currentBinIndices,[c_latitude,c_longitude]);
end


% ========================================================================
%%% Loading Low-Impact-Angle Data
% ========================================================================
% -------------------------------------------------
% Load Data
% -------------------------------------------------
%%% Load low impact data
lowImpactMat = dlmread(lowImpactAngleFile,',',1,0);

if size(lowImpactMat) ~= [1,1] % if there are any low-angle impacts:
    
    % -------------------------------------------------
    % Propagate low-impact-angle trajectories
    % -------------------------------------------------
    %%% Preallocating
    lowTrajs{size(lowImpactMat,1)} = [];

    %%% Setting time vector
    t_i = 0; % sec
    t_f = 6*pi; % Long bc events are watching for impact or escape
    dt = t_f/10000;
    time0_n = t_i:dt:t_f;

    %%% Choosing ode tolerance
    tol = 1e-10;

    %%% Setting integrator options
    options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

    %%% Setting necessary parameters for integration
    prms.u = secondary.MR;
    prms.R2_n = secondary.R_n;
    prms.L1x = L123(1,1);
    prms.L2x = L123(2,1);


    for kk = 1:size(lowImpactMat,1)
        %%% Propagating trajectory
        [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, lowImpactMat(kk,1:6)', options_ImpactEscape, prms);

        lowTrajs{kk} = [X_BCR_n, time_n];
    end

    figure; hold all
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
    axis equal
    for kk = 1:length(lowTrajs)
        plot3(lowTrajs{kk}(:,1),lowTrajs{kk}(:,2),lowTrajs{kk}(:,3))
    end
end % size(lowImpactMat) ~= [1,1]
% ========================================================================
%%% Determining necessary data for plotting
% ========================================================================

% -------------------------------------------------
% Find Y-Z contour points
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
[ yzContourPoints ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC_scInitial );





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
            zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
            JCs_xy(xk,yk) = zv;
        end
    end

    figure; hold all
%     plot3(r0s(:,1),r0s(:,2),r0s(:,3),'.')
    plot3(ones(size(yzContourPoints,2),1).*L123(Lpoint,1), yzContourPoints(1,:),yzContourPoints(2,:),'k','linewidth',3)
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
%     quiver3(X0s(:,1),X0s(:,2),X0s(:,3),X0s(:,4),X0s(:,5),X0s(:,6))
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
            zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
            JCs_xy(xk,yk) = zv;
        end
    end

    figure; hold all
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
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
    fill(yzContourPoints(1,:),yzContourPoints(2,:),binColors_NeckSections(4,:))

    plot(yzContourPoints(1,:), yzContourPoints(2,:),'k','linewidth',3)

    for kk = binCount_NeckSections-1:-1:1
        plotBody2(bins_neckSections(kk+1),[0,0,0], binColors_NeckSections(kk,:), binColors_NeckSections(kk,:), 1)
    end

    axis equal
    PlotBoi3('Y','Z','X',14)
end

% -------------------------------------------------
% Lat/Lon plot colored by initial neck section
% -------------------------------------------------
figure; hold all
for kk = 1:4
    if isempty(binData_neckSections(kk).latLons) == 0
        plot(binData_neckSections(kk).latLons(:,2),binData_neckSections(kk).latLons(:,1),...
            '.','markersize',15,'color',binColors_NeckSections(kk,:))
    end
end
title('Neck Sections')
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
cbar1.Ticks        = sort(0:binCount_NeckSections-1);
cbar1.TickLabels = num2cell(bins_neckSections(1:end-1).*rNorm);
cbar1.Label.String = {sprintf('|r_0| from L_%1.0f, km',Lpoint)};
cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
cbar1.Label.Position = [.2, 4.2, 0];
colormap(binColors_NeckSections)


% -------------------------------------------------
% Lat/Lon plot colored by impact angle
% -------------------------------------------------
figure; hold all
for kk = binCount_ImpactAngles:-1:1
if isempty(binData_impactAngles(kk).latLons) == 0
plot(binData_impactAngles(kk).latLons(:,2),binData_impactAngles(kk).latLons(:,1),...
'.','markersize',15,'color',binColors_ImpactAngles(kk,:))
end
end
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


































