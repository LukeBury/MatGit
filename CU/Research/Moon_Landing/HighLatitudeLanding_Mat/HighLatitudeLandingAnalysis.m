% ========================================================================
%%% Description
% ========================================================================
% For analyzing the high-latitude-landing-problem data generated in Julia

% Created: 03/20/20
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
plot_allTrajs       = false;
plot_periapseMap    = true;
plot_apoapseMap     = true;
plot_RPeriapse      = true;
plot_RApoapse       = true;
plot_xyCrossing     = true;
plot_xyCrossing_3D  = true;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Load the data generated in Julia
% -------------------------------------------------
%%% Data Path
dataPath = '/Users/lukebury/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/GoodFullData/';

%%% Data files
file_lowImpactTrajData          = [dataPath, 'lowImpactTrajs_150mps_4pi.txt'];
file_lowImpactEvents_periapses  = [dataPath, 'lowImpactTrajsPeriapses_150mps_4pi.txt'];
file_lowImpactEvents_apoapses   = [dataPath, 'lowImpactTrajsApoapses_150mps_4pi.txt'];
file_lowImpactEvents_zEquals0   = [dataPath, 'lowImpactTrajsZEquals0_150mps_4pi.txt'];
file_highImpactTrajData         = [dataPath, 'highImpactTrajs_150mps_4pi_stopOnImpact.txt'];
file_highImpactEvents_periapses = [dataPath, 'highImpactTrajsPeriapses_150mps_4pi_stopOnImpact.txt'];
file_highImpactEvents_apoapses  = [dataPath, 'highImpactTrajsApoapses_150mps_4pi_stopOnImpact.txt'];
file_highImpactEvents_zEquals0  = [dataPath, 'highImpactTrajsZEquals0_150mps_4pi_stopOnImpact.txt'];

%%% Load data
lowImpactTrajData_BaCR          = dlmread(file_lowImpactTrajData,',',1,0);
lowImpactEvents_BaCR_periapses  = dlmread(file_lowImpactEvents_periapses,',',1,0);
lowImpactEvents_BaCR_apoapses   = dlmread(file_lowImpactEvents_apoapses,',',1,0);
lowImpactEvents_BaCR_zEquals0   = dlmread(file_lowImpactEvents_zEquals0,',',1,0);
highImpactTrajData_BaCR         = dlmread(file_highImpactTrajData,',',1,0);
highImpactEvents_BaCR_periapses = dlmread(file_highImpactEvents_periapses,',',1,0);
highImpactEvents_BaCR_apoapses  = dlmread(file_highImpactEvents_apoapses,',',1,0);
highImpactEvents_BaCR_zEquals0  = dlmread(file_highImpactEvents_zEquals0,',',1,0);

%%% Number of unique trajectories
n_lowTrajs  = sum(isnan(lowImpactTrajData_BaCR(:,1)));
n_highTrajs = sum(isnan(highImpactTrajData_BaCR(:,1)));

%%% Columb specifiers
c_trajID = 1;
c_x      = 2;
c_y      = 3;
c_z      = 4;
c_xd     = 5;
c_yd     = 6;
c_zd     = 7;
c_t      = 8;
c_JC     = 9;
c_R      = 10;
c_a      = 11;
c_e      = 12;
c_i      = 13;
c_E      = 14;
c_H      = 15;
c_hx     = 16;
c_hy     = 17;
c_hz     = 18;
c_lat    = 19;
c_lon    = 20;

% -------------------------------------------------
%%% System
% -------------------------------------------------
primary = bodies.jupiter;   secondary = bodies.europa;

% -------------------------------------------------
%%% Normalizing factors and equillibrium points
% -------------------------------------------------
%%% Normalizing factors
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Collinear equilibrium points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% -------------------------------------------------
%%% Store trajectory data in separate cell arrays
% -------------------------------------------------
%%% Initialize cell arrays
CA_lowImpactTrajData_BaCR          = cell(n_lowTrajs,1);
CA_lowImpactEvents_BaCR_periapses  = cell(n_lowTrajs,1);
CA_lowImpactEvents_BaCR_apoapses  = cell(n_lowTrajs,1);
CA_lowImpactEvents_BaCR_zEquals0   = cell(n_lowTrajs,1);
CA_highImpactTrajData_BaCR         = cell(n_highTrajs,1);
CA_highImpactEvents_BaCR_periapses = cell(n_highTrajs,1);
CA_highImpactEvents_BaCR_apoapses = cell(n_highTrajs,1);
CA_highImpactEvents_BaCR_zEquals0  = cell(n_highTrajs,1);

%%% Loop through trajectories and store the data in cell arrays
for trajID = 1:n_lowTrajs
    CA_lowImpactTrajData_BaCR{trajID}         = lowImpactTrajData_BaCR(lowImpactTrajData_BaCR(:,c_trajID) == trajID,:);
    CA_lowImpactEvents_BaCR_periapses{trajID} = lowImpactEvents_BaCR_periapses(lowImpactEvents_BaCR_periapses(:,c_trajID) == trajID,:);
    CA_lowImpactEvents_BaCR_apoapses{trajID}  = lowImpactEvents_BaCR_apoapses(lowImpactEvents_BaCR_apoapses(:,c_trajID) == trajID,:);
    CA_lowImpactEvents_BaCR_zEquals0{trajID}  = lowImpactEvents_BaCR_zEquals0(lowImpactEvents_BaCR_zEquals0(:,c_trajID) == trajID,:);
end

for trajID = 1:n_highTrajs
    CA_highImpactTrajData_BaCR{trajID}         = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,c_trajID) == trajID,:);
    CA_highImpactEvents_BaCR_periapses{trajID} = highImpactEvents_BaCR_periapses(highImpactEvents_BaCR_periapses(:,c_trajID) == trajID,:);
    CA_highImpactEvents_BaCR_apoapses{trajID}  = highImpactEvents_BaCR_apoapses(highImpactEvents_BaCR_apoapses(:,c_trajID) == trajID,:);
    CA_highImpactEvents_BaCR_zEquals0{trajID}  = highImpactEvents_BaCR_zEquals0(highImpactEvents_BaCR_zEquals0(:,c_trajID) == trajID,:);
end

% -------------------------------------------------
%%% Combining data in large matrices for mass-plotting
% -------------------------------------------------
all_lowImpactTrajData_BaCR         = cell2mat(CA_lowImpactTrajData_BaCR);
all_lowImpactEvents_BaCR_periapses = cell2mat(CA_lowImpactEvents_BaCR_periapses);
all_lowImpactEvents_BaCR_apoapses  = cell2mat(CA_lowImpactEvents_BaCR_apoapses);
all_lowImpactEvents_BaCR_zEquals0  = cell2mat(CA_lowImpactEvents_BaCR_zEquals0);

all_highImpactTrajData_BaCR         = cell2mat(CA_highImpactTrajData_BaCR);
all_highImpactEvents_BaCR_periapses = cell2mat(CA_highImpactEvents_BaCR_periapses);
all_highImpactEvents_BaCR_apoapses  = cell2mat(CA_highImpactEvents_BaCR_apoapses);
all_highImpactEvents_BaCR_zEquals0  = cell2mat(CA_highImpactEvents_BaCR_zEquals0);


% -------------------------------------------------
%%% Calculating position-norms wrt secondary body at events
% -------------------------------------------------
RS_lowIimpactEvents_periapsis = rowNorm(all_lowImpactEvents_BaCR_periapses(:,c_x:c_z) - [1-secondary.MR,0,0]);
RS_lowIimpactEvents_apoapsis  = rowNorm(all_lowImpactEvents_BaCR_apoapses(:,c_x:c_z) - [1-secondary.MR,0,0]);
RS_lowIimpactEvents_zEquals0  = rowNorm(all_lowImpactEvents_BaCR_zEquals0(:,c_x:c_z) - [1-secondary.MR,0,0]);


RS_highIimpactEvents_periapsis = rowNorm(all_highImpactEvents_BaCR_periapses(:,c_x:c_z) - [1-secondary.MR,0,0]);
RS_highIimpactEvents_apoapsis  = rowNorm(all_highImpactEvents_BaCR_apoapses(:,c_x:c_z) - [1-secondary.MR,0,0]);
RS_highIimpactEvents_zEquals0  = rowNorm(all_highImpactEvents_BaCR_zEquals0(:,c_x:c_z) - [1-secondary.MR,0,0]);

% ========================================================================
%%% Plotting isolated trajectories
% ========================================================================
% -------------------------------------------------
%%% Discrete low trajs
% -------------------------------------------------
figure; hold all
title('Low Trajs')
PlotBoi3_CR3Bn(20)
axis equal
plotSecondary(secondary)
for kk = 1:1
    plot3(CA_lowImpactTrajData_BaCR{kk}(:,c_x),CA_lowImpactTrajData_BaCR{kk}(:,c_y),CA_lowImpactTrajData_BaCR{kk}(:,c_z),'color',colors.blue,'linewidth',2)
    plot3(CA_lowImpactEvents_BaCR_periapses{kk}(:,c_x),CA_lowImpactEvents_BaCR_periapses{kk}(:,c_y),CA_lowImpactEvents_BaCR_periapses{kk}(:,c_z),'k.','markersize',18)
    plot3(CA_lowImpactEvents_BaCR_apoapses{kk}(:,c_x),CA_lowImpactEvents_BaCR_apoapses{kk}(:,c_y),CA_lowImpactEvents_BaCR_apoapses{kk}(:,c_z),'k.','markersize',18)
%     plot3(CA_lowImpactEvents_zEquals0{kk}(:,c_x),CA_lowImpactEvents_zEquals0{kk}(:,c_y),CA_lowImpactEvents_zEquals0{kk}(:,c_z),'k.','markersize',18)
end

% -------------------------------------------------
%%% Discrete high trajs
% -------------------------------------------------
figure; hold all
title('High Trajs')
PlotBoi3_CR3Bn(20)
axis equal
% plotSecondary(secondary)
% for kk = 6622
for kk = 6918
    plot3(CA_highImpactTrajData_BaCR{kk}(:,c_x),CA_highImpactTrajData_BaCR{kk}(:,c_y),CA_highImpactTrajData_BaCR{kk}(:,c_z),'color',colors.red,'linewidth',2)
    plot3(CA_highImpactEvents_BaCR_periapses{kk}(:,c_x),CA_highImpactEvents_BaCR_periapses{kk}(:,c_y),CA_highImpactEvents_BaCR_periapses{kk}(:,c_z),'b.','markersize',18)
    plot3(CA_highImpactEvents_BaCR_apoapses{kk}(:,c_x),CA_highImpactEvents_BaCR_apoapses{kk}(:,c_y),CA_highImpactEvents_BaCR_apoapses{kk}(:,c_z),'g.','markersize',18)
%     plot3(CA_highImpactEvents_zEquals0{kk}(:,c_x),CA_highImpactEvents_zEquals0{kk}(:,c_y),CA_highImpactEvents_zEquals0{kk}(:,c_z),'k.','markersize',18)
end


% -------------------------------------------------
%%% All trajs
% -------------------------------------------------
if plot_allTrajs
    figure; hold all
    title('All Trajs')
    PlotBoi3_CR3Bn(20)
    axis equal
    plotSecondary(secondary)
    for kk = 1:n_lowTrajs
        plot3(CA_lowImpactTrajData_BaCR{kk}(:,c_x),CA_lowImpactTrajData_BaCR{kk}(:,c_y),CA_lowImpactTrajData_BaCR{kk}(:,c_z),'color',colors.blue,'linewidth',2)
    end
    for kk = 1:n_highTrajs
        plot3(CA_highImpactTrajData_BaCR{kk}(:,c_x),CA_highImpactTrajData_BaCR{kk}(:,c_y),CA_highImpactTrajData_BaCR{kk}(:,c_z),'color',colors.red,'linewidth',2)
    end

    xlim([0.97, 1.03])
    ylim([-1, 1].*0.02)
end % plot_allTrajs

% ========================================================================
%%% Plotting full sets
% ========================================================================
% -------------------------------------------------
%%% Periapse Mapping
% -------------------------------------------------
if plot_periapseMap
    figure; hold all
    title('Periapse Mapping')
    PlotBoi3_CR3Bn(20)
    axis equal
    p_low  = plot3(all_lowImpactEvents_BaCR_periapses(:,c_x),all_lowImpactEvents_BaCR_periapses(:,c_y),all_lowImpactEvents_BaCR_periapses(:,c_z),'b.','markersize',4);
    p_high = plot3(all_highImpactEvents_BaCR_periapses(:,c_x),all_highImpactEvents_BaCR_periapses(:,c_y),all_highImpactEvents_BaCR_periapses(:,c_z),'r.','markersize',20);
    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
end % plot_periapseMap
% -------------------------------------------------
%%% Apoapse Mapping
% -------------------------------------------------
if plot_apoapseMap
    figure; hold all
    title('Apoapse Mapping')
    PlotBoi3_CR3Bn(20)
    axis equal
    p_low  = plot3(all_lowImpactEvents_BaCR_apoapses(:,c_x),all_lowImpactEvents_BaCR_apoapses(:,c_y),all_lowImpactEvents_BaCR_apoapses(:,c_z),'b.','markersize',4);
    p_high = plot3(all_highImpactEvents_BaCR_apoapses(:,c_x),all_highImpactEvents_BaCR_apoapses(:,c_y),all_highImpactEvents_BaCR_apoapses(:,c_z),'r.','markersize',4);
    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
end % plot_apoapseMap
% -------------------------------------------------
%%% |R| at Periapsis Mapping
% -------------------------------------------------
if plot_RPeriapse
    figure; hold all
    title('|R| of Periapses')
    PlotBoi2('RS$_{per}$','RS$_{per}$',20,'LaTex')
    axis equal
    p_low  = plot(RS_lowIimpactEvents_periapsis,RS_lowIimpactEvents_periapsis,'b.','markersize',8);
    p_high = plot(RS_highIimpactEvents_periapsis,RS_highIimpactEvents_periapsis,'r.','markersize',8);

    p_rad  = plot([0, secondary.R_n],[secondary.R_n, secondary.R_n],'k','linewidth',1);
    plot([secondary.R_n, secondary.R_n],[0, secondary.R_n],'k','linewidth',1)

    xlim([0, max([max(RS_highIimpactEvents_periapsis), max(RS_lowIimpactEvents_periapsis)])])
    ylim([0, max([max(RS_highIimpactEvents_periapsis), max(RS_lowIimpactEvents_periapsis)])])

    [legh,objh] = legend([p_low, p_high, p_rad],'Low','High', 'Europa Radius');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end
% -------------------------------------------------
%%% |R| at Apoapsis Mapping
% -------------------------------------------------
if plot_RApoapse
    figure; hold all
    title('|R| of Apoapses')
    PlotBoi2('RS$_{apo}$','RS$_{apo}$',20,'LaTex')
    axis equal
    p_low  = plot(RS_lowIimpactEvents_apoapsis,RS_lowIimpactEvents_apoapsis,'b.','markersize',8);
    p_high = plot(RS_highIimpactEvents_apoapsis,RS_highIimpactEvents_apoapsis,'r.','markersize',8);

    p_rad  = plot([0, secondary.R_n],[secondary.R_n, secondary.R_n],'k','linewidth',1);
    plot([secondary.R_n, secondary.R_n],[0, secondary.R_n],'k','linewidth',1)

    xlim([0, max([max(RS_highIimpactEvents_apoapsis), max(RS_lowIimpactEvents_apoapsis)])])
    ylim([0, max([max(RS_highIimpactEvents_apoapsis), max(RS_lowIimpactEvents_apoapsis)])])


    [legh,objh] = legend([p_low, p_high, p_rad],'Low','High', 'Europa Radius');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end % plot_RApoapse

% -------------------------------------------------
%%% X-Y Plane Crossing
% -------------------------------------------------
if plot_xyCrossing
    figure; hold all
    title('X-Y Plane Crossing')
    PlotBoi3_CR3Bn(20)
    axis equal
    p_low  = plot(all_lowImpactEvents_BaCR_zEquals0(:,c_x),all_lowImpactEvents_BaCR_zEquals0(:,c_y),'b.','markersize',4);
    p_high = plot(all_highImpactEvents_BaCR_zEquals0(:,c_x),all_highImpactEvents_BaCR_zEquals0(:,c_y),'r.','markersize',4);
    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% X-Y Plane Crossing - 3D
% -------------------------------------------------
if plot_xyCrossing_3D

%     all_lowImpactEvents_BaCR_zEquals0_allData = NaN(size(all_lowImpactEvents_BaCR_zEquals0,1),size(all_lowImpactTrajData_BaCR,2));
%     all_highImpactEvents_BaCR_zEquals0_allData = NaN(size(all_highImpactEvents_BaCR_zEquals0,1),size(all_highImpactTrajData_BaCR,2));
%     
%     for kk = 1:size(all_lowImpactEvents_BaCR_zEquals0,1)
%         [z,y]=ismember(all_lowImpactTrajData_BaCR(:,[1,8]),all_lowImpactEvents_BaCR_zEquals0(kk,[1,8]),'rows');
%         if sum(z) == 1
%             all_lowImpactEvents_BaCR_zEquals0_allData(kk,:) = all_lowImpactTrajData_BaCR(z, :);
%         elseif sum(z)>1
%             temp = all_lowImpactTrajData_BaCR(z, :);
%             all_lowImpactEvents_BaCR_zEquals0_allData(kk,:) = temp(1,:);
%         end
%     end
%     
%     for kk = 1:size(all_highImpactEvents_BaCR_zEquals0,1)
%         [z,y]=ismember(all_highImpactTrajData_BaCR(:,[1,8]),all_highImpactEvents_BaCR_zEquals0(kk,[1,8]),'rows');
%         all_highImpactEvents_BaCR_zEquals0_allData(kk,:) = all_highImpactTrajData_BaCR(z, :);
%     end
    
    all_lowImpactEvents_BaCR_zEquals0_allData  = dlmread([dataPath, 'lowImpactTrajsZEquals0_full.txt'],',',1,0);
    all_highImpactEvents_BaCR_zEquals0_allData = dlmread([dataPath, 'highImpactTrajsZEquals0_full.txt'],',',1,0);
    
    
    figure; hold all
    title('X-Y Plane Crossing')
    axis normal
    PlotBoi3('$x_n$','$y_n$','$a$',20,'LaTeX')
    z_prm_low  = all_lowImpactEvents_BaCR_zEquals0_OEs(:,3);
    z_prm_high = all_highImpactEvents_BaCR_zEquals0_OEs(:,3);
    p_low  = plot3(all_lowImpactEvents_BaCR_zEquals0_allData(:,c_x),all_lowImpactEvents_BaCR_zEquals0_allData(:,c_y), all_lowImpactEvents_BaCR_zEquals0_allData(:,c_a),'b.','markersize',4);
    p_high = plot3(all_highImpactEvents_BaCR_zEquals0_allData(:,c_x),all_highImpactEvents_BaCR_zEquals0_allData(:,c_y),all_highImpactEvents_BaCR_zEquals0_allData(:,c_a),'r.','markersize',4);
    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% ========================================================================
%%% Checking end-times of high trajs to see if most were impacting or not
% ========================================================================
high_endTimes = NaN(length(CA_highImpactTrajData_BaCR),1);
for kk = 1:length(CA_highImpactTrajData_BaCR)
%     kk
%     CA_highImpactTrajData{kk}(end,8)
    high_endTimes(kk) = CA_highImpactTrajData_BaCR{kk}(end,8);
end
% average end time of trajectories
max(high_endTimes)
mean(high_endTimes)

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

% % % savepath = '/Users/lukebury/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/GoodFullData';
% % % filename_high = fullfile(savepath, 'highImpactTrajsZEquals0_full.txt');
% % % filename_low  = fullfile(savepath, 'lowImpactTrajsZEquals0_full.txt');
% % %  
% % % f_low  = fopen(filename_low, 'wt');
% % % f_high = fopen(filename_high, 'wt');
% % %  
% % % fprintf(f_low, 'trajID,x,y,z,xd,yd,zd,t,JC,R,a,e,i,E,H,h_x,h_y,h_z,lat,lon\n');
% % % fprintf(f_high, 'trajID,x,y,z,xd,yd,zd,t,JC,R,a,e,i,E,H,h_x,h_y,h_z,lat,lon\n');
% % %  
% % % for kk = 1:size(all_lowImpactEvents_BaCR_zEquals0_allData,1)
% % %      fprintf(f_low, '%1.0d,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.3f,%1.3f\n',...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,1), all_lowImpactEvents_BaCR_zEquals0_allData(kk,2), all_lowImpactEvents_BaCR_zEquals0_allData(kk,3),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,4), all_lowImpactEvents_BaCR_zEquals0_allData(kk,5), all_lowImpactEvents_BaCR_zEquals0_allData(kk,6),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,7), all_lowImpactEvents_BaCR_zEquals0_allData(kk,8), all_lowImpactEvents_BaCR_zEquals0_allData(kk,9),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,10), all_lowImpactEvents_BaCR_zEquals0_allData(kk,11), all_lowImpactEvents_BaCR_zEquals0_allData(kk,12),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,13), all_lowImpactEvents_BaCR_zEquals0_allData(kk,14), all_lowImpactEvents_BaCR_zEquals0_allData(kk,15),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,16), all_lowImpactEvents_BaCR_zEquals0_allData(kk,17), all_lowImpactEvents_BaCR_zEquals0_allData(kk,18),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,19), all_lowImpactEvents_BaCR_zEquals0_allData(kk,20));
% % % end
% % %  
% % % for kk = 1:size(all_highImpactEvents_BaCR_zEquals0_allData,1)
% % %      fprintf(f_high, '%1.0d,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.3f,%1.3f\n',...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,1), all_highImpactEvents_BaCR_zEquals0_allData(kk,2), all_highImpactEvents_BaCR_zEquals0_allData(kk,3),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,4), all_highImpactEvents_BaCR_zEquals0_allData(kk,5), all_highImpactEvents_BaCR_zEquals0_allData(kk,6),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,7), all_highImpactEvents_BaCR_zEquals0_allData(kk,8), all_highImpactEvents_BaCR_zEquals0_allData(kk,9),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,10), all_highImpactEvents_BaCR_zEquals0_allData(kk,11), all_highImpactEvents_BaCR_zEquals0_allData(kk,12),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,13), all_highImpactEvents_BaCR_zEquals0_allData(kk,14), all_highImpactEvents_BaCR_zEquals0_allData(kk,15),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,16), all_highImpactEvents_BaCR_zEquals0_allData(kk,17), all_highImpactEvents_BaCR_zEquals0_allData(kk,18),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,19), all_highImpactEvents_BaCR_zEquals0_allData(kk,20));
% % % end










