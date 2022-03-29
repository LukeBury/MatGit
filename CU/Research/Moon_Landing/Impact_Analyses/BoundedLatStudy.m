%%% INFORMATION
% When grid searches are performed over intial conditions in the L1 or L2
% necks at low energies, we find that ballistic impact trajectories only 
% reach certain latitudes until a threshold energy is reached which allow
% impacts of all latitudes. Curiously, for sufficiently low energies, the
% bounds set on reachable latitudes don't come from the energy domain -
% these trajectories have a sufficient jacobi constant to reach any point
% on the surface. So what part of phase space describes this boundary? In
% this script, I plan to investigate angular velocity as a boundary, and
% will also study the maximum latitude reached before impact.
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
tic

% ========================================================================
%%% Run Switches
% ========================================================================
plot_angularMomentum = 0;
plot_fullLatLon      = 0;
plot_trajectories    = 1;

% ========================================================================
% Data sets
% ========================================================================
%%% Low-angle Enceladus
dataSets_enc = {'F.iGS_encL2_13mps_4km_149v0s';...
                'F.iGS_encL2_26mps_4km_149v0s';...
                'F.iGS_encL2_39mps_4km_149v0s';...
                'F.iGS_encL2_52mps_4km_149v0s';...
                'F.iGS_encL2_65mps_4km_149v0s';...
                'F.iGS_encL2_78mps_4km_149v0s';...
                'F.iGS_encL2_91mps_4km_149v0s'};
            
%%% Low-angle Europa
dataSets_eur = {'F.iGS_eurL2_50mps_50km_149v0s';...
                'F.iGS_eurL2_100mps_50km_149v0s';...
                'F.iGS_eurL2_150mps_50km_149v0s';...
                'F.iGS_eurL2_200mps_50km_149v0s';...
                'F.iGS_eurL2_250mps_50km_149v0s';...
                'F.iGS_eurL2_300mps_50km_149v0s';...
                'F.iGS_eurL2_350mps_50km_149v0s'};

%%% Full set Europa
dataSets_X0s = {'impactX0s_iGS_eurL2_50mps_250km_22v0s';...
                'impactX0s_iGS_eurL2_100mps_250km_22v0s';...
                'impactX0s_iGS_eurL2_150mps_250km_22v0s';...
                'impactX0s_iGS_eurL2_200mps_250km_22v0s';...
                'impactX0s_iGS_eurL2_250mps_250km_22v0s';...
                'impactX0s_iGS_eurL2_300mps_250km_22v0s';...
                'impactX0s_iGS_eurL2_350mps_250km_22v0s'};

% ========================================================================
% Choose data set and run
% ========================================================================            
%%% Choose data set
dataSets = dataSets_X0s;  
n_dataSets = length(dataSets);
% 989
% n_dataSets = 1;

%%% File location
MatlabOutputsPath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs/';

%%% Preallocating structure for multiple-data-set data
dataSetData{n_dataSets}.AzIs_deg    = [];
dataSetData{n_dataSets}.latLons_deg = [];
dataSetData{n_dataSets}.rBCRs       = [];

% for dataSet_i = 1:n_dataSets
for dataSet_i = [1]
dataSet = dataSets{dataSet_i};

if contains(dataSet,'impactX0s') == 0
    %%% Low-angle landing file
    dataFile = [MatlabOutputsPath, dataSet,'_land.txt'];
elseif contains(dataSet,'impactX0s') == 1
    %%% GeneralX0 impact file
    dataFile = [MatlabOutputsPath, dataSet,'_data.txt'];
end

% ========================================================================
%%% Importing Data and Setting up System
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Bodies
primary = bodies.jupiter;
secondary = bodies.europa;

if contains(secondary.name,'eur') == 1 && contains(dataFile,'eur') == 0
    warning('Body discrepancy');
    return
elseif contains(secondary.name,'enc') == 1 && contains(dataFile,'enc') == 0
    warning('Body discrepancy');
    return
end

%%% Colinear Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec


% ========================================================================
%%% Studying Angular Momentum
% ========================================================================
% ---------------------------------------
% Grabbing initial conditions
% ---------------------------------------
%%% Load low impact data
ImpactX0sMat = dlmread(dataFile,',',1,0);

%%% If this data set is empty, skip to next data set
if isequal(ImpactX0sMat,0)
    continue
end

% ---------------------------------------
% Setting up for integration
% ---------------------------------------
%%% Time
t_i = 0; % sec
t_f = 2*pi;
n_dt = 1000;
time0_n = linspace(t_i,t_f,n_dt);

%%% Choosing ode45 tolerance
tol = 2.22045e-14;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Setting extra parameters
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);

%%% How many trajectories to loop through?
nTraj = size(ImpactX0sMat,1);
% nTraj = 50;

%%% Preallocating for storage of all trajectories
trajs(nTraj).azIs    = [];
trajs(nTraj).latLons = [];
trajs(nTraj).rBCR    = [];


%%% Looping through trajectories
for ii = 1:nTraj
%%% Clearing / Preallocating
    X_BCR_n = [];
    
    %%% Initial conditions
    if contains(dataSet,'impactX0s') == 0
        X0_n = ImpactX0sMat(ii,2:7)';
    elseif contains(dataSet,'impactX0s') == 1
        X0_n = [L123(2,1),ImpactX0sMat(ii,6:10)]';
    end

    % ---------------------------------------
    % Propagating trajectory
    % ---------------------------------------
    %%% Propagating trajectory
    [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
                time0_n, X0_n, options_ImpactEscape, prms);
    
    if plot_trajectories == 1
        %%% Storing trajectory
        trajs(ii).rBCR = [X_BCR_n(:,1:3); nan(1,3)];
    end
    
    % ---------------------------------------
    % Angular momentum and Lat/Lon
    % ---------------------------------------
    %%% Preallocating
    trajs(ii).azIs = zeros(size(X_BCR_n,1),2);
    trajs(ii).latLons = zeros(size(X_BCR_n,1),2);
    trajs(ii).maxLat       = 0;
    trajs(ii).impactLat    = 0;

    for kk = 1:size(X_BCR_n,1)
        if plot_angularMomentum == 1
            %%% Secondary-centric position coordinates
            rSCR = X_BCR_n(kk,1:3) - [1-secondary.MR, 0, 0];
    %         angMom = cross(rSCR, X_BCR_n(kk,4:6));

            %%% Angular momentum
            [ r_BCI ] = r_BCR2BCI( X_BCR_n(kk,1:3), time_n(kk), 1 );
            [ r_BCI_secondary ] = r_BCR2BCI( [1-secondary.MR,0,0], time_n(kk), 1 );
            [ v_BCI ] = v_BCR2BCI( X_BCR_n(kk,4:6), r_BCI, time_n(kk), 1 );
            [ v_BCI_secondary ] = v_BCR2BCI( [0, 0, 0], r_BCI_secondary, time_n(kk), 1 );

            r_SCI = r_BCI - r_BCI_secondary;
            v_SCI = v_BCI - v_BCI_secondary;

            %%% This one's super boring
%             angMom = cross(r_BCI, v_BCI);

            angMom = cross(r_SCI, v_SCI);



            %%% Determine azimuth and elevation at current time, conver to
            %%% inclination
            [Az_rad,El_rad,Rad] = cart2sph(angMom(1), angMom(2), angMom(3));
            i_rad = pi/2 - El_rad;

            %%% Store azimuth and inclination of angular momentum vector
            trajs(ii).azIs(kk,:) = [Az_rad, i_rad];
        end % plot_angularMomentum
        
        if plot_fullLatLon == 1
            %%% Getting latitude/longitude at current time
%             r_SCR_kk = X_BCR_n(kk,1:3) - [1-secondary.MR, 0, 0];
            [lat_deg, lon_deg] = BCR2latlon(X_BCR_n(kk,1:3), 'secondary', secondary.MR);
%             [lat_deg, lon_deg] = surfCR3BP2latlon(r_SCR_kk, 'secondary', secondary.MR);

            %%% Storing lat/lon
            trajs(ii).latLons(kk,:) = [lat_deg, lon_deg];
            
            
        end
    end
    
    if plot_fullLatLon == 1
        %%% Storing maximum lat
        trajs(ii).maxLat = max(abs(trajs(ii).latLons(:,1)));
        
        %%% Storing impact lat
        trajs(ii).impactLat = trajs(ii).latLons(end,1);
    end
end % ii = 1:size(lowImpactMat,1)

%%% Rearranging data and storing to mega-structure
if plot_angularMomentum == 1
    AzIs_deg = vertcat(trajs.azIs).*(180/pi);
    dataSetData{dataSet_i}.AzIs_deg = AzIs_deg;
end
if plot_fullLatLon == 1
    latLons_deg = vertcat(trajs.latLons);
    dataSetData{dataSet_i}.latLons_deg = latLons_deg;
    
    %%% Finding overall max lat
    overallMaxLat = 0;
    maxImpactLat  = 0;
    for ii = 1:length(trajs)
        if trajs(ii).maxLat > overallMaxLat
            overallMaxLat = trajs(ii).maxLat;
        end
        if trajs(ii).impactLat > maxImpactLat
            maxImpactLat = trajs(ii).impactLat;
        end
    end
end
if plot_trajectories == 1
    rBCRs = vertcat(trajs.rBCR);
    dataSetData{dataSet_i}.rBCRs = rBCRs;
end

% ---------------------------------------
% Plotting angular momentum
% ---------------------------------------
if plot_angularMomentum == 1
    figure; hold all
    plot(dataSetData{dataSet_i}.AzIs_deg(:,1),dataSetData{dataSet_i}.AzIs_deg(:,2),'b.')
    PlotBoi2('Azimuth, $^\circ$','Inclination, $^\circ$',16,'LaTex')
    xlim([-180 180])
    ylim([0 180])
end
% ---------------------------------------
% Plotting lat/lons
% ---------------------------------------
if plot_fullLatLon == 1
    figure; hold all
    plot(dataSetData{dataSet_i}.latLons_deg(:,2),dataSetData{dataSet_i}.latLons_deg(:,1),'b.')
    PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
    xlim([-180 180])
    ylim([-90 90])
    
    fprintf('Maximum Latitude Reached During Flight: %1.2f Deg\n',overallMaxLat)
    fprintf('Maximum Latitude at Impact:             %1.2f Deg\n',maxImpactLat)
end
% ---------------------------------------
% Plotting system
% ---------------------------------------
if plot_trajectories == 1
    figure; hold all
    plot3(dataSetData{dataSet_i}.rBCRs(:,1),dataSetData{dataSet_i}.rBCRs(:,2),dataSetData{dataSet_i}.rBCRs(:,3),'b','linewidth',1)
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
    axis equal
end



% figure; hold all
% dataColors = colorScale([colors.std.cyan; colors.std.mag],n_dataSets);
% for kk = n_dataSets:-1:1
%     plot(dataSetData{kk}.AzIs_deg(:,1),dataSetData{kk}.AzIs_deg(:,2),'.','color',dataColors(kk,:));
%     angPlotHandles{kk} = plot(nan,nan,'-','color',dataColors(kk,:),'linewidth',2);
% end
% % angPlotHandles = fliplr(angPlotHandles);
% legend([angPlotHandles{:}],'50 mps','100 mps','150 mps','200 mps','250 mps','300 mps','350 mps')
% PlotBoi2('Azimuth, $^\circ$','Inclination, $^\circ$',16,'LaTex')
% xlim([-180 180])
% ylim([0 180])


end % dataSet_i = 1:n_dataSets


toc





