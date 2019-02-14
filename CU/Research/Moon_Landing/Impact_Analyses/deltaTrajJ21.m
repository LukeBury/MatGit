%%% Information
% In the primary grid search, "trajectory IDs" are store for the results so
% that the same trajectory can be compared with and without J2. This
% comparison is the purpose of this script
% ========================================================================

clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
ticStart = tic;

% ========================================================================
%%% Run Switches
% ========================================================================
plot_angularMomentum = 1;
plot_fullLatLon      = 1;
plot_trajectories    = 1;

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

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% ========================================================================
% Data sets
% ========================================================================
% ---------------------------------------
% Data sets
% ---------------------------------------
dataSets_eur = {'F.iGS_eurL2_50mps_50km_149v0s';...
                'F.iGS_eurL2_100mps_50km_149v0s';...
                'F.iGS_eurL2_150mps_50km_149v0s';...
                'F.iGS_eurL2_200mps_50km_149v0s';...
                'F.iGS_eurL2_250mps_50km_149v0s';...
                'F.iGS_eurL2_300mps_50km_149v0s';...
                'F.iGS_eurL2_350mps_50km_149v0s'};
            
dataSets_eur_J2 = {'F.iGS_eurL2_J21Comp_50mps_50km_149v0s';...
                   'F.iGS_eurL2_J21Comp_100mps_50km_149v0s';...
                   'F.iGS_eurL2_J21Comp_150mps_50km_149v0s';...
                   'F.iGS_eurL2_J21Comp_200mps_50km_149v0s';...
                   'F.iGS_eurL2_J21Comp_250mps_50km_149v0s';...
                   'F.iGS_eurL2_J21Comp_300mps_50km_149v0s';...
                   'F.iGS_eurL2_J21Comp_350mps_50km_149v0s'};            

% ---------------------------------------
% Data info
% ---------------------------------------
%%% File location
MatlabOutputsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/';

%%% Number of data sets
n_dataSets = length(dataSets_eur);
% 989
% n_dataSets = 3;

% ---------------------------------------
% Grabbing full filenames
% ---------------------------------------
%%% Preallocating - All impact data
fileNames_eur_data{n_dataSets}    = '';
fileNames_eur_J2_data{n_dataSets} = '';

%%% Preallocating - Low-angle impactors
fileNames_eur_low{n_dataSets}    = '';
fileNames_eur_J2_low{n_dataSets} = '';

%%% Looping through data sets and storing full filenames
for currentIndex = 1:n_dataSets
    %%% All impact data
    fileNames_eur_data{currentIndex}    = [MatlabOutputsPath,dataSets_eur{currentIndex},'_data.txt'];
    fileNames_eur_J2_data{currentIndex} = [MatlabOutputsPath,dataSets_eur_J2{currentIndex},'_data.txt'];
    
    %%% Low-angle impactors
    fileNames_eur_low{currentIndex}     = [MatlabOutputsPath,dataSets_eur{currentIndex},'_land.txt'];
    fileNames_eur_J2_low{currentIndex}  = [MatlabOutputsPath,dataSets_eur_J2{currentIndex},'_land.txt'];
end


% ========================================================================
% Looping through data sets to study differences caused by J2
% ========================================================================
% 989
% for dataset_i = 1:2
for dataset_i = 1:n_dataSets
% ---------------------------------------
% Turning data into matrices
% ---------------------------------------
%%% Load impact data
dataMat_full    = dlmread(fileNames_eur_data{dataset_i},',',1,0);
dataMat_J2_full = dlmread(fileNames_eur_J2_data{dataset_i},',',1,0);

% dataMat_low1    = dlmread(fileNames_eur_low{1},',',1,0);
% dataMat_J2_low1 = dlmread(fileNames_eur_J2_low{1},',',1,0);

% ---------------------------------------
% Column Identifiers
% ---------------------------------------
%%% Full data sets
c_full_ID             = 1;
c_full_binImpactAngle = 2;
c_full_binNeckSection = 3;
c_full_lat            = 4;
c_full_lon            = 5;
c_full_IA             = 6;
c_full_tf             = 7;
c_full_y0             = 8;
c_full_z0             = 9;
c_full_v0Az           = 10;
c_full_v0El           = 11;

%%% Low-angle data sets
%trajID,x0_n,y0_n,z0_n,dx0_n,dy0_n,dz0_n,bin_neckSection,latitude,longitude,endTime
c_low_ID             = 1;
c_low_x0             = 2;
c_low_y0             = 3;
c_low_z0             = 4;
c_low_xd0            = 5;
c_low_yd0            = 6;
c_low_zd0            = 7;
c_low_binNeckSection = 8;
c_low_lat            = 9;
c_low_lon            = 10;
c_low_tf             = 11;

% ---------------------------------------
% Looping through IDs to find difference
% ---------------------------------------
%%% Preallocating for differences
d_y0s            = nan(size(dataMat_full,1),1);
d_z0s            = nan(size(dataMat_full,1),1);
d_lats           = nan(size(dataMat_full,1),1);
d_lons           = nan(size(dataMat_full,1),1);
d_Az_latLons_deg = nan(size(dataMat_full,1),1);
distanceMags_km  = nan(size(dataMat_full,1),1);

%%% Preallocating for other purposes
lats = nan(size(dataMat_full,1),1);
lons = nan(size(dataMat_full,1),1);
impactQuarter = nan(size(dataMat_full,1),1);



%%% Finding unique and shared IDs with indices
[IDs_J2Has_normalDoesnt, Indices_J2_unique]     = setdiff(dataMat_J2_full(:,c_full_ID),dataMat_full(:,c_full_ID));
[IDs_normalHas_J2Doesnt, Indices_normal_unique] = setdiff(dataMat_full(:,c_full_ID),   dataMat_J2_full(:,c_full_ID));
[IDs_bothHave, Indices_normal_shared]           = intersect(dataMat_full(:,c_full_ID), dataMat_J2_full(:,c_full_ID));

%%% Transposing because Matlab only loops through row vectors
IDs_J2Has_normalDoesnt = IDs_J2Has_normalDoesnt'; Indices_J2_unique     = Indices_J2_unique';
IDs_normalHas_J2Doesnt = IDs_normalHas_J2Doesnt'; Indices_normal_unique = Indices_normal_unique';
IDs_bothHave = IDs_bothHave';                     Indices_normal_shared = Indices_normal_shared';



%%% Looping through IDs shared by each data set
for currentIndex = Indices_normal_shared
    %%% Current ID
    ID = dataMat_full(currentIndex,c_full_ID);
    
    %%% Search J2 data and look for current ID
    matchingJ2Index = find(dataMat_J2_full(:,c_full_ID) == ID);
    
    %%% If current ID not present in J2 data, skip to next ID
    if isempty(matchingJ2Index) == 1
        continue
    end
    
    %%% Storing current lat/lon
    lats(currentIndex) = dataMat_full(currentIndex,c_full_lat);
    lons(currentIndex) = dataMat_full(currentIndex,c_full_lon);
    
    %%% Difference in starting position
    d_y0s(currentIndex) = dataMat_J2_full(matchingJ2Index,c_full_y0) - dataMat_full(currentIndex,c_full_y0);
    d_z0s(currentIndex) = dataMat_J2_full(matchingJ2Index,c_full_z0) - dataMat_full(currentIndex,c_full_z0);
    
    %%% Difference in impact latitude/longitude
    d_lats(currentIndex) = dataMat_J2_full(matchingJ2Index,c_full_lat) - dataMat_full(currentIndex,c_full_lat);
    d_lons(currentIndex) = dataMat_J2_full(matchingJ2Index,c_full_lon) - dataMat_full(currentIndex,c_full_lon);
    
    %%% Determining quarter (bc results mostly symmetric across lat = 0)
    impactQuarter(currentIndex) = determineQuarter(lats(currentIndex),lons(currentIndex));
    
    
    %%% Azimuth of delta lat/lon
    d_Az_latLons_deg(currentIndex) = azimuth(0, 0, d_lats(currentIndex),d_lons(currentIndex));
    
    %%% Magnitude of distance impact site moved by J2
    [rImpact_BCR_nominal] = latlon2BCR('secondary', dataMat_full(currentIndex,c_full_lat), dataMat_full(currentIndex,c_full_lon), secondary.MR, secondary.R_n);
    [rImpact_BCR_J2]      = latlon2BCR('secondary', dataMat_J2_full(matchingJ2Index,c_full_lat), dataMat_J2_full(matchingJ2Index,c_full_lon), secondary.MR, secondary.R_n);
    distanceMags_km(currentIndex) = norm(rImpact_BCR_J2 - rImpact_BCR_nominal)*rNorm;
    
end

% ========================================================================
% Plotting
% ========================================================================
% % ---------------------------------------
% % Plotting full lat/lon mapping with deltas
% % ---------------------------------------
% figure; hold all
% plot(lons,lats,'b.','markersize',7)
% plot([lons, lons + d_lons]',[lats, lats + d_lats]','r','linewidth',0.5)
% PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',16,'LaTex')
% xlim([-180 180])
% ylim([-90 90])

% ---------------------------------------
% Plotting longitudinal-quarter vs azimuth of delta-lat/lon
% ---------------------------------------
figure; hold all
plot(impactQuarter, d_Az_latLons_deg,'b.','markersize',5)
PlotBoi2('Impact Quarter','Azimuth of Lat/Lon Change, $^\circ$',16,'LaTex')
xlim([0.5 4.5])
xticks([1 2 3 4])
ylim([0 360])

% ---------------------------------------
% Plotting histogram of azimuths of delta-lat/lon
% ---------------------------------------
figure; hold all
histogram(d_Az_latLons_deg,250,'facecolor','b','facealpha',0.7)
PlotBoi2('Azimuth of Lat/Lon Change, $^\circ$','Frequency',16,'LaTex')

% figure
% [p,x] = hist(d_Az_latLons_deg); plot(x,p/sum(p)); %PDF
% figure
% [f,x] = ecdf(d_Az_latLons_deg); plot(x,f);
% ---------------------------------------
% Plotting histogram of magnitudes of impact location changes
% ---------------------------------------
figure; hold all
histogram(distanceMags_km,250,'facecolor','b','facealpha',0.7)
PlotBoi2('$\Delta$ Impact Location, $km$','Frequency',16,'LaTex')

% ---------------------------------------
% Bar graph of impact numbers for each dynamic model
% ---------------------------------------
figure
c = categorical({'Unique to Standard CR3BP','Unique to J2-CR3BP','Both'});
bar(c,[length(IDs_normalHas_J2Doesnt), length(IDs_J2Has_normalDoesnt), length(IDs_bothHave)])
PlotBoi2('','Number of Impacts',16,'LaTex')

end % dataset_i = 1:n_dataSets







toc(ticStart)


% ========================================================================
% ========================================================================
% Functions
% ========================================================================
% ========================================================================
function [octant] = determineOctant(lat,lon)
%%% Inputs
% lat - latitude (deg)
% lon - longitude (deg)
%
%%% Outputs
% octant - which octant on the lat/lon map the point is in
%
% 1 | 2 | 3 | 4
% -------------
% 5 | 6 | 7 | 8
%
% ========================================================================
if lat >= 0
    if (lon >= -180) && (lon < -90)
        octant = 1;
    elseif (lon >= -90) && (lon < 0)
        octant = 2;
    elseif (lon >= 0) && (lon < 90)
        octant = 3;
    elseif (lon >= 90) && (lon <= 180)
        octant = 4;
    end
elseif lat < 0
    if (lon >= -180) && (lon < -90)
        octant = 5;
    elseif (lon >= -90) && (lon < 0)
        octant = 6;
    elseif (lon >= 0) && (lon < 90)
        octant = 7;
    elseif (lon >= 90) && (lon <= 180)
        octant = 8;
    end
end

end

function [quarter] = determineQuarter(lat,lon)
%%% Inputs
% lat - latitude (deg)
% lon - longitude (deg)
%
%%% Outputs
% quarter - which quarter on the lat/lon map the point is in
%
% 1 | 2 | 3 | 4
% ========================================================================
if (lon >= -180) && (lon < -90)
    quarter = 1;
elseif (lon >= -90) && (lon < 0)
    quarter = 2;
elseif (lon >= 0) && (lon < 90)
    quarter = 3;
elseif (lon >= 90) && (lon <= 180)
    quarter = 4;
end

end











